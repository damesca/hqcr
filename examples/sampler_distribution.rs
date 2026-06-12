//! Empirical distribution tester for the three samplers in `poly::sampling`.
//!
//! Drives `sample_fixed_weight` (secret x/y), `sample_fixed_weight_mod`
//! (ephemeral r1/r2/e) and `sample_uniform` (public h) over many independent
//! random XOF seeds, accumulates how often each position in `[0, N)` is set,
//! and writes a self-contained interactive HTML report so the per-position
//! distribution — and any bias — can be inspected by eye and via a chi-square
//! goodness-of-fit statistic.
//!
//! For an *unbiased* fixed-weight sampler every position is a Bernoulli trial
//! with success probability `weight / N`, so the per-position frequency curve
//! should be flat at `weight/N` and `chi2/dof ≈ 1`. For `sample_uniform` each
//! bit is Bernoulli(1/2). The `mod` sampler is *expected* to show mild
//! structure at the very low indices (its backward dedup maps collisions to
//! small positions), and this report is precisely how to visualize that.
//!
//! Run (NEVER executed automatically — provided for the user to run):
//!
//!   cargo run --release --example sampler_distribution                 # HQC-128, 20000 trials
//!   cargo run --release --example sampler_distribution -- 192 50000    # HQC-192, 50000 trials
//!   cargo run --release --example sampler_distribution -- 256 20000
//!
//! Output: `sampler_distribution_<param>.html` in the current directory.

use std::time::{SystemTime, UNIX_EPOCH};

use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;

use hqcr::params::{Hqc128, Hqc192, Hqc256, HqcParams};
use hqcr::poly::sampling::{sample_fixed_weight, sample_fixed_weight_mod, sample_uniform};
use hqcr::poly::Poly;

// ── Tunable defaults ────────────────────────────────────────────────────────
const DEFAULT_TRIALS: usize = 20_000;
/// Number of position-bins for the per-position frequency curve. Positions are
/// grouped into this many contiguous bins and the mean frequency per bin is
/// plotted; this both smooths the curve and keeps the embedded JSON small.
const POS_BINS: usize = 512;
/// Number of buckets for the histogram of per-position set-counts.
const COUNT_BUCKETS: usize = 60;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let param = args.get(1).map(|s| s.as_str()).unwrap_or("128");
    let trials: usize = args
        .get(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(DEFAULT_TRIALS);

    // A fresh per-run nonce so successive runs use independent seeds, while a
    // single run is internally reproducible (seed = nonce || role || trial).
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos() as u64;

    eprintln!("HQC-{param}: running {trials} trials per sampler (nonce={nonce})…");

    let report = match param {
        "128" => run::<Hqc128>(trials, nonce),
        "192" => run::<Hqc192>(trials, nonce),
        "256" => run::<Hqc256>(trials, nonce),
        other => {
            eprintln!("unknown param set '{other}' (expected 128 | 192 | 256)");
            std::process::exit(1);
        }
    };

    let out = format!("sampler_distribution_{param}.html");
    let html = render_html(param, trials, nonce, &report);
    std::fs::write(&out, html).expect("failed to write HTML report");

    eprintln!("\nSummary (chi2/dof ≈ 1.0 means consistent with the ideal flat distribution):");
    for s in &report {
        eprintln!(
            "  {:<24} weight={:<4} E[p]={:.6}  observed mean={:.6}  min={:.6} max={:.6}  chi2/dof={:.4}",
            s.name, s.weight, s.expected_p, s.mean_p, s.min_p, s.max_p, s.chi2_over_dof
        );
    }
    eprintln!("\nWrote {out}  — open it in a browser.");
}

// ── Per-sampler accumulated result ──────────────────────────────────────────
struct SamplerStats {
    name: &'static str,
    weight: usize,
    n: usize,
    trials: usize,
    expected_p: f64,
    mean_p: f64,
    min_p: f64,
    max_p: f64,
    stddev_p: f64,
    chi2: f64,
    dof: usize,
    chi2_over_dof: f64,
    /// Mean per-position frequency in each of POS_BINS contiguous position bins.
    pos_curve: Vec<f64>,
    /// Histogram of per-position set-counts: (bucket_center_count, num_positions).
    count_hist: Vec<(f64, u64)>,
    /// Mean/variance of the count distribution, plus the ideal binomial overlay.
    count_mean: f64,
    count_std: f64,
    binom_mean: f64,
    binom_std: f64,
}

fn run<P: HqcParams>(trials: usize, nonce: u64) -> Vec<SamplerStats> {
    let n = P::N;

    let mut fw = vec![0u64; n]; // sample_fixed_weight  (secret, weight = OMEGA)
    let mut fm = vec![0u64; n]; // sample_fixed_weight_mod (ephemeral, weight = OMEGA_R)
    let mut un = vec![0u64; n]; // sample_uniform (each bit ~ Bernoulli(1/2))

    let report_every = (trials / 10).max(1);
    for t in 0..trials {
        let mut x1 = seeded_xof(nonce, 0, t);
        accumulate::<P>(&sample_fixed_weight::<P>(&mut x1, P::OMEGA), &mut fw);

        let mut x2 = seeded_xof(nonce, 1, t);
        accumulate::<P>(&sample_fixed_weight_mod::<P>(&mut x2, P::OMEGA_R), &mut fm);

        let mut x3 = seeded_xof(nonce, 2, t);
        accumulate::<P>(&sample_uniform::<P>(&mut x3), &mut un);

        if (t + 1) % report_every == 0 {
            eprintln!("  {}/{} trials", t + 1, trials);
        }
    }

    vec![
        summarize("sample_fixed_weight", P::OMEGA, n, trials, &fw),
        summarize("sample_fixed_weight_mod", P::OMEGA_R, n, trials, &fm),
        // For the uniform sampler the "weight" is the expected number of set
        // bits, N/2; we pass N so expected_p resolves to 1/2 below.
        summarize("sample_uniform", (n + 1) / 2, n, trials, &un),
    ]
}

/// Independent SHAKE256 stream per (run, sampler role, trial index).
fn seeded_xof(nonce: u64, role: u8, t: usize) -> impl XofReader {
    let mut h = Shake256::default();
    h.update(&nonce.to_le_bytes());
    h.update(&[role]);
    h.update(&(t as u64).to_le_bytes());
    h.finalize_xof()
}

/// Add one sample's set bits into the per-position counter.
fn accumulate<P: HqcParams>(p: &Poly<P>, counts: &mut [u64]) {
    for (i, c) in counts.iter_mut().enumerate() {
        *c += p.get_bit(i);
    }
}

fn summarize(
    name: &'static str,
    weight: usize,
    n: usize,
    trials: usize,
    counts: &[u64],
) -> SamplerStats {
    let t = trials as f64;
    let expected_p = weight as f64 / n as f64;
    let exp_count = expected_p * t; // expected set-count per position
    // Variance of a Binomial(T, p) count is T·p·(1−p). Using this as the
    // chi-square denominator gives E[χ²/dof] ≈ 1.0 for all sampler types,
    // including sample_uniform (p=0.5). The old denominator T·p was correct
    // only when p≪1 (fixed-weight samplers); for p=0.5 it halved every term,
    // making χ²/dof≈0.5 even for a perfectly uniform sampler.
    let chi2_variance = (t * expected_p * (1.0 - expected_p)).max(f64::EPSILON);

    // Per-position probability stats.
    let mut sum = 0.0f64;
    let mut sum_sq = 0.0f64;
    let mut min_p = f64::INFINITY;
    let mut max_p = f64::NEG_INFINITY;
    let mut chi2 = 0.0f64;
    let mut min_c = u64::MAX;
    let mut max_c = 0u64;
    for &c in counts {
        let p = c as f64 / t;
        sum += p;
        sum_sq += p * p;
        min_p = min_p.min(p);
        max_p = max_p.max(p);
        {
            let d = c as f64 - exp_count;
            chi2 += d * d / chi2_variance;
        }
        min_c = min_c.min(c);
        max_c = max_c.max(c);
    }
    let nf = n as f64;
    let mean_p = sum / nf;
    let var_p = (sum_sq / nf - mean_p * mean_p).max(0.0);
    let dof = n.saturating_sub(1);
    let chi2_over_dof = if dof > 0 { chi2 / dof as f64 } else { 0.0 };

    // Per-position frequency curve, binned across the index range.
    let bins = POS_BINS.min(n);
    let mut pos_curve = vec![0.0f64; bins];
    for (b, slot) in pos_curve.iter_mut().enumerate() {
        let lo = b * n / bins;
        let hi = ((b + 1) * n / bins).max(lo + 1).min(n);
        let s: u64 = counts[lo..hi].iter().sum();
        *slot = s as f64 / ((hi - lo) as f64 * t);
    }

    // Histogram of per-position set-counts.
    let (count_hist, count_mean, count_std) =
        build_count_hist(counts, min_c, max_c, t);

    // Ideal binomial overlay: each position is Binomial(trials, expected_p).
    let binom_mean = exp_count;
    let binom_std = (t * expected_p * (1.0 - expected_p)).max(0.0).sqrt();

    SamplerStats {
        name,
        weight,
        n,
        trials,
        expected_p,
        mean_p,
        min_p,
        max_p,
        stddev_p: var_p.sqrt(),
        chi2,
        dof,
        chi2_over_dof,
        pos_curve,
        count_hist,
        count_mean,
        count_std,
        binom_mean,
        binom_std,
    }
}

fn build_count_hist(
    counts: &[u64],
    min_c: u64,
    max_c: u64,
    t: f64,
) -> (Vec<(f64, u64)>, f64, f64) {
    let buckets = COUNT_BUCKETS;
    let span = (max_c - min_c).max(1);
    let mut hist = vec![0u64; buckets];
    let mut mean = 0.0f64;
    let mut mean_sq = 0.0f64;
    for &c in counts {
        let idx = (((c - min_c) * buckets as u64) / (span + 1)) as usize;
        hist[idx.min(buckets - 1)] += 1;
        let cf = c as f64;
        mean += cf;
        mean_sq += cf * cf;
    }
    let nf = counts.len() as f64;
    mean /= nf;
    let std = (mean_sq / nf - mean * mean).max(0.0).sqrt();

    let out = hist
        .into_iter()
        .enumerate()
        .map(|(i, freq)| {
            let center =
                min_c as f64 + (i as f64 + 0.5) * (span as f64 + 1.0) / buckets as f64;
            (center, freq)
        })
        .collect();
    // mean/std reported as fractions of a probability are not needed here; keep
    // raw counts. `t` retained for signature symmetry.
    let _ = t;
    (out, mean, std)
}

// ── HTML report ──────────────────────────────────────────────────────────────

fn render_html(param: &str, trials: usize, nonce: u64, report: &[SamplerStats]) -> String {
    let data_json = report
        .iter()
        .map(sampler_json)
        .collect::<Vec<_>>()
        .join(",\n");

    format!(
        r##"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>HQC-{param} sampler distributions</title>
<style>
  :root {{ color-scheme: dark; }}
  body {{ font-family: ui-monospace, Menlo, Consolas, monospace; background:#0f1115; color:#d7dae0; margin:0; padding:24px; }}
  h1 {{ font-size:18px; margin:0 0 4px; }}
  .sub {{ color:#8a90a0; font-size:12px; margin-bottom:24px; }}
  .card {{ background:#161922; border:1px solid #242a38; border-radius:10px; padding:16px 18px; margin-bottom:28px; }}
  .card h2 {{ font-size:15px; margin:0 0 2px; color:#9ecbff; }}
  .stats {{ display:flex; flex-wrap:wrap; gap:6px 22px; font-size:12px; color:#aab1c2; margin:8px 0 14px; }}
  .stats b {{ color:#e8ebf2; font-weight:600; }}
  .verdict {{ font-weight:700; }}
  .ok {{ color:#5fd38d; }}
  .warn {{ color:#f0c674; }}
  .bad {{ color:#ff6b6b; }}
  .charts {{ display:grid; grid-template-columns:1fr 1fr; gap:18px; }}
  @media (max-width:900px) {{ .charts {{ grid-template-columns:1fr; }} }}
  .chart-title {{ font-size:12px; color:#8a90a0; margin-bottom:4px; }}
  canvas {{ width:100%; height:240px; background:#0c0e13; border-radius:6px; display:block; }}
  .legend {{ font-size:11px; color:#8a90a0; margin-top:6px; }}
  .legend span {{ display:inline-block; margin-right:14px; }}
  .swatch {{ display:inline-block; width:10px; height:10px; border-radius:2px; vertical-align:middle; margin-right:5px; }}
</style>
</head>
<body>
<h1>HQC-{param} — empirical sampler distributions</h1>
<div class="sub">{trials} trials per sampler &middot; run nonce {nonce} &middot; ideal &chi;&sup2;/dof &asymp; 1.0</div>
<div id="root"></div>

<script>
const DATA = [
{data_json}
];

function vClass(r) {{ return (r<1.05 && r>0.95) ? "ok" : (r<1.2 && r>0.8) ? "warn" : "bad"; }}
function vText(r)  {{ return (r<1.05 && r>0.95) ? "consistent with flat (unbiased)" : (r<1.2 && r>0.8) ? "mild deviation" : "significant deviation"; }}

// Generic line plotter on a high-DPI canvas.
function lineChart(cv, series, opts) {{
  const dpr = window.devicePixelRatio || 1;
  const W = cv.clientWidth, H = cv.clientHeight;
  cv.width = W*dpr; cv.height = H*dpr;
  const ctx = cv.getContext('2d'); ctx.scale(dpr,dpr);
  const m = {{l:58,r:12,t:10,b:24}};
  const pw = W-m.l-m.r, ph = H-m.t-m.b;
  let ymin=opts.ymin, ymax=opts.ymax;
  if (ymin===undefined||ymax===undefined) {{
    ymin=Infinity; ymax=-Infinity;
    for (const s of series) for (const y of s.y) {{ if(y<ymin)ymin=y; if(y>ymax)ymax=y; }}
  }}
  if (ymax===ymin) {{ ymax=ymin+1; }}
  const pad=(ymax-ymin)*0.08; ymin-=pad; ymax+=pad;
  const n = series[0].y.length;
  const X = i => m.l + (n<=1?0:pw*i/(n-1));
  const Y = v => m.t + ph*(1-(v-ymin)/(ymax-ymin));
  // grid + y labels
  ctx.strokeStyle='#222838'; ctx.fillStyle='#6b7282'; ctx.font='10px monospace'; ctx.lineWidth=1;
  for (let g=0; g<=4; g++) {{
    const v=ymin+(ymax-ymin)*g/4, y=Y(v);
    ctx.beginPath(); ctx.moveTo(m.l,y); ctx.lineTo(W-m.r,y); ctx.stroke();
    ctx.fillText(v.toExponential(2), 4, y+3);
  }}
  // x labels (position range)
  if (opts.xmax!==undefined) {{
    ctx.fillStyle='#6b7282';
    for (let g=0; g<=4; g++) {{
      const x=m.l+pw*g/4;
      ctx.fillText(Math.round(opts.xmax*g/4), x-8, H-8);
    }}
  }}
  for (const s of series) {{
    ctx.strokeStyle=s.color; ctx.lineWidth=s.width||1.5;
    if (s.dash) ctx.setLineDash(s.dash); else ctx.setLineDash([]);
    ctx.beginPath();
    s.y.forEach((v,i)=>{{ const x=X(i),y=Y(v); i?ctx.lineTo(x,y):ctx.moveTo(x,y); }});
    ctx.stroke();
  }}
  ctx.setLineDash([]);
}}

// Bar plotter for the count histogram, with an optional overlay curve.
function barChart(cv, bars, overlay) {{
  const dpr = window.devicePixelRatio || 1;
  const W = cv.clientWidth, H = cv.clientHeight;
  cv.width = W*dpr; cv.height = H*dpr;
  const ctx = cv.getContext('2d'); ctx.scale(dpr,dpr);
  const m = {{l:48,r:12,t:10,b:24}};
  const pw = W-m.l-m.r, ph = H-m.t-m.b;
  const xs = bars.map(b=>b[0]), ys = bars.map(b=>b[1]);
  let xmin=Math.min(...xs), xmax=Math.max(...xs);
  let ymax=Math.max(...ys,1); if(overlay) ymax=Math.max(ymax,...overlay.map(o=>o[1]));
  if (xmax===xmin) xmax=xmin+1;
  const X = v => m.l + pw*(v-xmin)/(xmax-xmin);
  const Y = v => m.t + ph*(1-v/ymax);
  ctx.strokeStyle='#222838'; ctx.fillStyle='#6b7282'; ctx.font='10px monospace';
  for (let g=0; g<=4; g++) {{
    const y=m.t+ph*g/4; ctx.beginPath(); ctx.moveTo(m.l,y); ctx.lineTo(W-m.r,y); ctx.stroke();
    ctx.fillText(Math.round(ymax*(1-g/4)), 4, y+3);
  }}
  for (let g=0; g<=4; g++) {{ const v=xmin+(xmax-xmin)*g/4; ctx.fillText(Math.round(v), m.l+pw*g/4-10, H-8); }}
  const bw = Math.max(1, pw/bars.length - 1);
  ctx.fillStyle='#4a8cff';
  for (const [x,y] of bars) {{ const px=X(x)-bw/2, py=Y(y); ctx.fillRect(px,py,bw,m.t+ph-py); }}
  if (overlay) {{
    ctx.strokeStyle='#ff9f43'; ctx.lineWidth=2; ctx.beginPath();
    overlay.forEach((o,i)=>{{ const px=X(o[0]),py=Y(o[1]); i?ctx.lineTo(px,py):ctx.moveTo(px,py); }});
    ctx.stroke();
  }}
}}

function gaussianOverlay(bars, mean, std, area) {{
  if (std<=0) return null;
  const k = area/(std*Math.sqrt(2*Math.PI));
  return bars.map(b=>{{ const z=(b[0]-mean)/std; return [b[0], k*Math.exp(-0.5*z*z)]; }});
}}

const root = document.getElementById('root');
for (const d of DATA) {{
  const card = document.createElement('div'); card.className='card';
  const vc = vClass(d.chi2_over_dof);
  card.innerHTML = `
    <h2>${{d.name}}</h2>
    <div class="stats">
      <span>weight <b>${{d.weight}}</b></span>
      <span>N <b>${{d.n}}</b></span>
      <span>E[p]=weight/N <b>${{d.expected_p.toExponential(4)}}</b></span>
      <span>observed mean <b>${{d.mean_p.toExponential(4)}}</b></span>
      <span>min <b>${{d.min_p.toExponential(3)}}</b></span>
      <span>max <b>${{d.max_p.toExponential(3)}}</b></span>
      <span>&sigma; <b>${{d.stddev_p.toExponential(3)}}</b></span>
      <span>&chi;&sup2;=${{d.chi2.toFixed(1)}} / dof=${{d.dof}} = <b class="verdict ${{vc}}">${{d.chi2_over_dof.toFixed(4)}}</b></span>
      <span class="verdict ${{vc}}">${{vText(d.chi2_over_dof)}}</span>
    </div>
    <div class="charts">
      <div>
        <div class="chart-title">Per-position selection frequency (binned across [0, N))</div>
        <canvas class="c-curve"></canvas>
        <div class="legend"><span><i class="swatch" style="background:#4a8cff"></i>observed</span><span><i class="swatch" style="background:#5fd38d"></i>ideal weight/N</span></div>
      </div>
      <div>
        <div class="chart-title">Histogram of per-position set-counts vs ideal binomial</div>
        <canvas class="c-hist"></canvas>
        <div class="legend"><span><i class="swatch" style="background:#4a8cff"></i>observed positions</span><span><i class="swatch" style="background:#ff9f43"></i>Binomial(T, p)</span></div>
      </div>
    </div>`;
  root.appendChild(card);

  const curve = card.querySelector('.c-curve');
  lineChart(curve, [
    {{ y: d.pos_curve, color:'#4a8cff', width:1.2 }},
    {{ y: d.pos_curve.map(()=>d.expected_p), color:'#5fd38d', width:1.5, dash:[5,4] }},
  ], {{ xmax: d.n }});

  const hist = card.querySelector('.c-hist');
  const area = d.n * (d.count_hist.length>1 ? (d.count_hist[1][0]-d.count_hist[0][0]) : 1);
  barChart(hist, d.count_hist, gaussianOverlay(d.count_hist, d.binom_mean, d.binom_std, area));
}}
</script>
</body>
</html>
"##
    )
}

fn sampler_json(s: &SamplerStats) -> String {
    let pos_curve = f64_array(&s.pos_curve);
    let count_hist = s
        .count_hist
        .iter()
        .map(|(c, f)| format!("[{:.3},{}]", c, f))
        .collect::<Vec<_>>()
        .join(",");
    format!(
        r#"  {{ "name":"{}", "weight":{}, "n":{}, "trials":{}, "expected_p":{}, "mean_p":{}, "min_p":{}, "max_p":{}, "stddev_p":{}, "chi2":{}, "dof":{}, "chi2_over_dof":{}, "count_mean":{}, "count_std":{}, "binom_mean":{}, "binom_std":{}, "pos_curve":[{}], "count_hist":[{}] }}"#,
        s.name,
        s.weight,
        s.n,
        s.trials,
        s.expected_p,
        s.mean_p,
        s.min_p,
        s.max_p,
        s.stddev_p,
        s.chi2,
        s.dof,
        s.chi2_over_dof,
        s.count_mean,
        s.count_std,
        s.binom_mean,
        s.binom_std,
        pos_curve,
        count_hist,
    )
}

fn f64_array(v: &[f64]) -> String {
    v.iter()
        .map(|x| format!("{:.8}", x))
        .collect::<Vec<_>>()
        .join(",")
}
