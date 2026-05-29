# HQC — Hamming Quasi-Cyclic: Matemáticas, Esquema KEM e Implicaciones de Seguridad

> Notas de referencia para la implementación en Rust nativo de HQC,
> basadas en las especificaciones oficiales (versión 22/08/2025, `pqc-hqc.org`)
> y el paper fundacional [AGM+2018].

---

## Índice

1. [Contexto y referencias](#1-contexto-y-referencias)
2. [Álgebra subyacente: el anillo $\mathcal{R}$ y los códigos quasi-cíclicos](#2-álgebra-subyacente)
3. [El problema de seguridad: QCSD](#3-el-problema-de-seguridad-qcsd)
4. [HQC-PKE: cifrado IND-CPA](#4-hqc-pke)
5. [El código corrector de errores concatenado RMRS](#5-el-código-rmrs)
6. [HQC-KEM: transformación a IND-CCA2](#6-hqc-kem)
7. [Parámetros concretos](#7-parámetros-concretos)
8. [Análisis de seguridad](#8-análisis-de-seguridad)
9. [Implicaciones para la implementación en Rust](#9-implicaciones-para-la-implementación)
10. [Apéndice: referencias y glosario](#10-apéndice)

---

## 1. Contexto y referencias

### 1.1 Estado NIST

HQC fue seleccionado por NIST el **11 de marzo de 2025** como quinto algoritmo
de intercambio de clave post-cuántico a estandarizar, en calidad de alternativa
code-based al ya estandarizado ML-KEM (Kyber). El estándar FIPS definitivo se
espera alrededor de 2027.

### 1.2 Referencias primarias

| Documento | Descripción |
|---|---|
| **[AGM+2016]** Aguilar-Melchor, Blazy, Deneuville, Gaborit, Zémor. *Efficient Encryption from Random Quasi-Cyclic Codes.* IEEE Trans. Inf. Theory 64(5), 2018. arXiv:1612.05572 | Paper fundacional. Propone el framework de cifrado sobre QC codes y define HQC y RQC por primera vez. |
| **[HQC-Spec-2025]** Aguilar-Melchor et al. *Hamming Quasi-Cyclic (HQC).* Especificaciones oficiales v22/08/2025. https://pqc-hqc.org | Especificaciones completas, incluyendo la actualización al transform SFO⊥_m (salted FO), parámetros alineados con FIPS-203, y el análisis de DFR actual. |
| **[HHK17]** Hofheinz, Hövelmanns, Kiltz. *A Modular Analysis of the Fujisaki-Okamoto Transformation.* TCC 2017. | El transform FO que convierte PKE IND-CPA en KEM IND-CCA2. Usado como base del transform HHK. |
| **[Ber10]** Bernstein. *Grover vs. McEliece.* PQCrypto 2010. | Análisis del efecto de Grover sobre ataques ISD al syndrome decoding. |
| **[FS09]** Finiasz, Sendrier. *Security Bounds for the Design of Code-Based Cryptosystems.* ASIACRYPT 2009. | Bound de complejidad para ISD, esencial para el dimensionado de parámetros. |

### 1.3 La idea en una frase

HQC es un cifrado ElGamal-like sobre un anillo de polinomios binarios cíclicos,
donde el "ruido" que permite corregir errores en el descifrado se mantiene en
magnitud controlada gracias a la estructura quasi-cíclica del problema duro
subyacente.

---

## 2. Álgebra subyacente

### 2.1 El anillo base $\mathcal{R}$

Sea $n$ un **número primo primitivo** (es decir, tal que el polinomio
$(X^n - 1)/(X - 1)$ es irreducible en $\mathbb{F}_2[X]$). Define el anillo:

$$\mathcal{R} := \mathbb{F}_2[X]/(X^n - 1)$$

Los elementos de $\mathcal{R}$ son polinomios binarios de grado $< n$, o
equivalentemente vectores de $\mathbb{F}_2^n$. La multiplicación en $\mathcal{R}$
es la convolución cíclica módulo $X^n - 1$:

$$w_k = \sum_{i+j \equiv k \pmod{n}} u_i \cdot v_j \quad \forall k \in \{0, \ldots, n-1\}$$

La condición de que $n$ sea primo primitivo garantiza que $X^n - 1$ se
factoriza en exactamente **dos factores irreducibles** sobre $\mathbb{F}_2$:
$X^n - 1 = (X-1) \cdot f(X)$ donde $f$ tiene grado $n-1$ e irreducible.
Esto impide ataques algebraicos que explotan factorizaciones con muchos
factores pequeños.

### 2.2 La matriz circulante

Para $\mathbf{v} \in \mathbb{F}_2^n$, la **matriz circulante** inducida por
$\mathbf{v}$ es:

$$\mathrm{rot}(\mathbf{v}) = \begin{pmatrix} v_0 & v_{n-1} & \cdots & v_1 \\ v_1 & v_0 & \cdots & v_2 \\ \vdots & & \ddots & \vdots \\ v_{n-1} & v_{n-2} & \cdots & v_0 \end{pmatrix} \in \mathbb{F}_2^{n \times n}$$

La propiedad fundamental es que la multiplicación en $\mathcal{R}$ equivale a
producto matriz-vector:

$$\mathbf{u} \cdot \mathbf{v} = \mathbf{u} \times \mathrm{rot}(\mathbf{v})^\top = \mathbf{v} \times \mathrm{rot}(\mathbf{u})^\top$$

Esto permite que toda la clave pública sea compacta: una matriz quasi-cíclica
de índice 2 queda completamente descrita por un único polinomio de $\mathcal{R}$.

### 2.3 Peso de Hamming y vectores de peso fijo

Para $\mathbf{v} \in \mathbb{F}_2^n$, su **peso de Hamming** $\omega(\mathbf{v})$
es el número de coordenadas no nulas.

Para un entero positivo $\omega$, el conjunto de vectores de peso exactamente
$\omega$ se denota:

$$\mathcal{R}_\omega := \{ \mathbf{v} \in \mathcal{R} \mid \omega(\mathbf{v}) = \omega \}$$

El producto $\mathbf{e} \cdot \mathbf{r}$ con $\mathbf{e}, \mathbf{r} \in \mathcal{R}_\omega$ tiene
peso de Hamming acotado. Esto es la piedra angular del análisis de descifrado
correcto:

$$\omega(\mathbf{e} \cdot \mathbf{r}) \leq \omega(\mathbf{e}) \cdot \omega(\mathbf{r}) = \omega^2$$

(cota superior; la media real es $n/2 \cdot (1 - (1-2\omega/n)^2)$ aproximadamente).

### 2.4 Códigos Quasi-Cíclicos (QC)

Un código **Quasi-Cíclico de índice $s$** de longitud $sn$ es aquel para el que
aplicar un desplazamiento circular simultáneo a todos los $s$ bloques de $n$
posiciones produce otro codeword. En términos polinomiales: si
$\mathbf{c} = (c_0, \ldots, c_{s-1}) \in \mathcal{C}$, entonces
$(X \cdot c_0, \ldots, X \cdot c_{s-1}) \in \mathcal{C}$.

Un código quasi-cíclico **sistemático de índice $s$ y tasa $1/s$** tiene
matriz de paridad:

$$H = \begin{pmatrix} I_n & 0 & \cdots & 0 & A_0 \\ 0 & I_n & \cdots & & A_1 \\ & & \ddots & & \vdots \\ 0 & \cdots & 0 & I_n & A_{s-2} \end{pmatrix}$$

donde cada $A_i$ es una matriz **circulante** $n \times n$ — es decir, está
completamente determinada por un polinomio $\mathbf{h}_i \in \mathcal{R}$.

Para HQC se utilizan códigos de índice $s = 2$ (HQC-128, HQC-192) y $s = 3$
(HQC-256). La clave privada contiene los polinomios $\mathbf{h}_i$ que
generan las matrices $A_i$.

---

## 3. El problema de seguridad: QCSD

### 3.1 Syndrome Decoding (SD): el problema duro raíz

**Definición (SD, versión computacional):** dados $H \in \mathbb{F}_2^{(n-k) \times n}$
aleatoria e $\mathbf{y} \in \mathbb{F}_2^{n-k}$, encontrar $\mathbf{x} \in \mathbb{F}_2^n$
de peso $\omega$ tal que:

$$H \mathbf{x}^\top = \mathbf{y}^\top$$

Este problema es **NP-completo** (Berlekamp-McEliece-van Tilborg, 1978) y
equivalente al problema LPN con un número fijo de muestras.

### 3.2 2-QCSD-P: la variante quasi-cíclica usada en HQC

HQC basa su seguridad en la variante quasi-cíclica con condición de paridad.
Sean $b_1 \in \{0,1\}$ y $b_2 = \omega + b_1 \cdot \omega \pmod 2$. Define:

$$\mathbb{F}_{2,b_1}^n = \{ \mathbf{h} \in \mathbb{F}_2^n \mid \omega(\mathbf{h}) \equiv b_1 \pmod 2 \}$$

**Distribución 2-QCSD-P$(n, \omega, b_1)$:**

1. Samplear $H = (I_n \mid \mathrm{rot}(\mathbf{h})) \in \mathbb{F}_{2,b_1}^{n \times 2n}$
2. Samplear $\mathbf{x} = (\mathbf{x}_1, \mathbf{x}_2)$ con $\omega(\mathbf{x}_1) = \omega(\mathbf{x}_2) = \omega$
3. Computar $\mathbf{y}^\top = H \mathbf{x}^\top$
4. Salida: $(H, \mathbf{y})$

**Problema decisional 2-DQCSD-P:** distinguir $(H, \mathbf{y})$ de esta
distribución de $(H, \mathbf{y})$ uniforme sobre $\mathbb{F}_{2,b_1}^{n \times 2n} \times \mathbb{F}_2^n$.

La **condición de paridad** ($b_1, b_2$) es una precaución técnica para
eliminar distinguidores triviales basados en la paridad del síndrome; tiene
un coste de a lo sumo 1 bit de seguridad.

### 3.3 3-QCSD-PT: variante con truncación (HQC-256)

Para HQC-256 se usa un código de índice 3, donde la longitud natural del
mensaje codificado es $n_1 n_2$ (no primo). HQC trabaja con $n$ el primo
primitivo inmediatamente mayor que $n_1 n_2$ y trunca $\ell = n - n_1 n_2$
bits finales. La variante del problema duro se llama **3-DQCSD-PT** (con
Paridad y Truncación).

La truncación rompe la estructura cíclica para los últimos $\ell$ bits, lo que
paradójicamente **debilita levemente** la ventaja del atacante que intenta
explotar la estructura QC.

---

## 4. HQC-PKE

### 4.1 Estructura matemática

HQC-PKE es un sistema de **cifrado probabilístico tipo ElGamal** en el anillo
$\mathcal{R}$, donde el "logaritmo discreto" es reemplazado por el syndrome
decoding.

Trabajamos en el par de anillos $\mathcal{R}^2 = \mathcal{R} \times \mathcal{R}$.
La clave privada es un par de vectores sparse $(\mathbf{x}, \mathbf{y}) \in \mathcal{R}_\omega^2$,
y la clave pública es el par $(\mathbf{h}, \mathbf{s})$ donde:

$$\mathbf{s} = \mathbf{x} + \mathbf{h} \cdot \mathbf{y} \in \mathcal{R}$$

con $\mathbf{h} \in \mathcal{R}$ una "base pública" aleatoria, y la igualdad
en el sentido de $\mathcal{R}$ (convolución cíclica módulo $X^n - 1$).

La seguridad de esta relación reduce a 2-DQCSD-P: el par $(\mathbf{h}, \mathbf{s})$
es computacionalmente indistinguible de un par $(h, s)$ uniforme en $\mathcal{R}^2$.

### 4.2 Generación de claves — HQC-PKE.Keygen

```
Entradas: ninguna
Salidas: (ek_PKE, dk_PKE)  // clave de cifrado, clave de descifrado

1. Samplear seed ←$ {0,1}^256
2. (seed_dk, seed_ek) = SHA3-512(seed)
3. h ←$ R   usando SHAKE256(seed_ek, n)          // base pública, uniforme
4. (x, y) ←$ R_ω × R_ω  usando SHAKE256(seed_dk) // clave privada, sparse
5. s = x + h·y   (convolución cíclica en R)
6. ek_PKE = (h, s)                               // clave pública
7. dk_PKE = (seed_dk, ek_PKE)                    // clave privada
```

**Observación sobre el almacenamiento:** $\mathbf{h}$ y $\mathbf{s}$ son polinomios en
$\mathcal{R}$ de $n$ bits. La clave pública es exactamente $2n$ bits — sin
expansión respecto al tama~no del problema. Esto es consecuencia directa de
usar códigos quasi-cíclicos: si se usara una matriz $H$ completamente
aleatoria en SD clásico, la clave pública sería $O(n^2)$ bits.

### 4.3 Cifrado — HQC-PKE.Encrypt

```
Entradas: ek_PKE = (h, s), m ∈ F_2^k, θ ∈ {0,1}^256
Salidas: ciphertext c = (u, v)

1. Codificar: m' = C.Encode(m)   (Reed-Muller/Reed-Solomon, longitud n bits)
2. Samplear r_1, r_2 ←$ R_ω, e ←$ R_ω   usando SHAKE256(θ)
3. u = r_1 + h·r_2                          // "componente pública"
4. v = m' + s·r_2 + e                       // "componente privada"
5. Output c = (u, v)
```

La estructura del ciphertext $(u, v)$ es la de un **ciphertext ElGamal**:
$u$ es una "clave de sesión enmascarada" y $v$ es el mensaje enmascarado con esa clave.

### 4.4 Descifrado — HQC-PKE.Decrypt

```
Entradas: dk_PKE = (x, y, h, s), c = (u, v)
Salidas: m ∈ F_2^k

1. Computar: v - u·y
   = (m' + s·r_2 + e) - (r_1 + h·r_2)·y
   = m' + (x + h·y)·r_2 + e - r_1·y - h·r_2·y
   = m' + x·r_2 + h·y·r_2 + e - r_1·y - h·y·r_2
   = m' + x·r_2 - r_1·y + e
   = m' + err                           // donde err = x·r_2 + r_1·y + e

2. Decodificar: m = C.Decode(m' + err)
```

**El análisis de descifrado correcto** depende del peso de Hamming del error:

$$\omega(\mathbf{err}) = \omega(\mathbf{x} \cdot \mathbf{r}_2 + \mathbf{r}_1 \cdot \mathbf{y} + \mathbf{e})$$

Cada término tiene peso acotado: $\omega(\mathbf{x} \cdot \mathbf{r}_2) \leq \omega^2$,
$\omega(\mathbf{r}_1 \cdot \mathbf{y}) \leq \omega^2$, $\omega(\mathbf{e}) = \omega$.
Por cota de unión, el peso total es a lo sumo $2\omega^2 + \omega$, que debe
ser menor que la capacidad de corrección del código $\mathcal{C}$.

Si el error excede la capacidad del decoder, el descifrado falla. La
**Decryption Failure Rate (DFR)** se analiza precisamente y se requiere que sea
$\leq 2^{-\lambda}$ para nivel de seguridad $\lambda$.

---

## 5. El código corrector de errores concatenado RMRS

El componente más complejo de HQC desde el punto de vista algorítmico es el
código que permite recuperar $\mathbf{m}'$ de $\mathbf{m}' + \mathbf{err}$.
La elección del código condiciona directamente la DFR y el tama~no de parámetros.

### 5.1 Historia de los decoders en HQC

Las especificaciones originales (2017) usaban un **código tensor BCH⊗Repetición**.
En 2020 se introdujo el decoder **RMRS** (Reed-Muller / Reed-Solomon concatenado),
que es estrictamente mejor en todos los aspectos y es el único que usa la
especificación actual.

### 5.2 Código Reed-Solomon como código externo

Un código **Reed-Solomon** $[n_2, k_2, d_2]$ sobre $\mathbb{F}_{2^m}$ con
$n_2 = 2^m - 1$, $k_2 = n_2 - 2t$, $d_2 = 2t + 1$ puede corregir hasta $t$
errores.

Los mensajes de $k$ bits de HQC se codifican primero como palabras de
Reed-Solomon sobre $\mathbb{F}_{2^m}$. Para HQC-128: $m = 8$, $n_2 = 255$
símbolos de 8 bits, $t = 57$ (puede corregir 57 errores de símbolo).

### 5.3 Código Reed-Muller como código interno

El código de **Reed-Muller** de primer orden $\mathcal{RM}(1, m)$ tiene parámetros
$[2^m, m+1, 2^{m-1}]$:

- Longitud: $n_1 = 2^m$ bits
- Dimensión: $m + 1$ bits (codifica $m+1$ bits por bloque)
- Distancia mínima: $2^{m-1}$ — puede corregir hasta $2^{m-2} - 1$ errores

La decodificación de RM(1, m) se realiza eficientemente con la **Fast Hadamard
Transform (FHT)** en $O(n_1 \log n_1)$:

$$(\hat{f}_0, \hat{f}_1, \ldots, \hat{f}_{n_1-1}) = \text{WHT}(f_0, f_1, \ldots, f_{n_1-1})$$

$$\hat{f}_k = \sum_{j=0}^{n_1-1} (-1)^{\langle k, j \rangle} f_j$$

donde $\langle \cdot, \cdot \rangle$ es el producto escalar binario de los
índices vistos como vectores de $m$ bits.

### 5.4 Código concatenado RMRS

El código **concatenado** $\mathcal{C} = \mathcal{C}_\text{RS} \circ \mathcal{C}_\text{RM}$:

1. Tomar el mensaje $\mathbf{m} \in \mathbb{F}_2^k$.
2. Codificar con Reed-Solomon: $\mathbf{m}_\text{RS} = \text{RS.Encode}(\mathbf{m}) \in \mathbb{F}_{2^m}^{n_2}$.
3. Interpretar cada símbolo $\mathbb{F}_{2^m}$ como un bloque de $m$ bits.
4. Codificar cada bloque con RM(1, m): $\mathbf{m}' = (\text{RM.Encode}(\text{sym}_1) \| \cdots \| \text{RM.Encode}(\text{sym}_{n_2})) \in \mathbb{F}_2^{n_1 n_2}$.

La longitud del mensaje codificado es $n_1 n_2$ bits, que se embebe en $\mathcal{R}$
de longitud $n$ (con truncación de los últimos $\ell = n - n_1 n_2$ bits).

**Decodificación:**

1. Dividir en $n_2$ bloques de $n_1 = 2^m$ bits.
2. Aplicar FHT a cada bloque para decodificar RM: obtener $n_2$ símbolos de $\mathbb{F}_{2^m}$.
3. Aplicar el decoder algebraico de Reed-Solomon (e.g., Berlekamp-Massey o
   Euclidean): corregir hasta $t$ errores de símbolo.
4. Salida: el mensaje original $\mathbf{m}$.

La concatenación hace que la capacidad de corrección **efectiva** sea mucho
mayor que cualquiera de los dos códigos por separado: RM maneja los errores
"de bit" dentro de cada símbolo, y RS corrige errores "de símbolo" (incluso
si el FHT falla completamente en algunos bloques, RS puede recuperar esos símbolos).

---

## 6. HQC-KEM

### 6.1 Del cifrado IND-CPA al KEM IND-CCA2

HQC-PKE es IND-CPA pero **no** IND-CCA2. Para convertirlo en un KEM IND-CCA2
se aplica la **transformación de Fujisaki-Okamoto (FO)**, específicamente la
variante **SFO⊥_m** (salted FO with implicit rejection for messages), introducida
en la versión 2025 de las especificaciones.

La idea central de FO es hacer que el ciphertext sea **determinístico dado el
mensaje y la clave pública**: el atacante que modifica un ciphertext obtiene un
resultado sin relación con el mensaje original, porque el receptor re-cifra y
compara.

### 6.2 Generación de claves — HQC-KEM.Keygen

```
Entradas: ninguna
Salidas: (ek_KEM, dk_KEM)

1. seed_KEM ←$ {0,1}^256
2. (seed_PKE.dk, seed_PKE.ek) = SHA3-512(seed_KEM)
3. (ek_PKE, dk_PKE) = HQC-PKE.Keygen(seed_PKE.dk, seed_PKE.ek)
4. ek_KEM = ek_PKE = (h, s)
5. dk_KEM = (seed_KEM, ek_KEM)
   // formato alternativo comprimido: dk_KEM = (seed_KEM)
```

### 6.3 Encapsulación — HQC-KEM.Encaps

```
Entradas: ek_KEM, mensaje aleatorio m ←$ {0,1}^256, salt ←$ {0,1}^256
Salidas: (K, c)  // clave compartida K, ciphertext c

1. m ←$ {0,1}^256  (aleatorio)
2. salt ←$ {0,1}^256  (aleatorio)
3. (K, θ) = G(m, ek_KEM, salt)  donde G = SHAKE256
   // K es la clave compartida candidata (256 bits)
   // θ es la semilla de aleatoriedad para cifrar
4. c = HQC-PKE.Encrypt(ek_KEM, m, θ)   con salt embebido
5. Output (K, c)
```

La derivación de $(K, \theta)$ con SHA3/SHAKE256 absorbe tanto $m$ como
$\text{ek}_\text{KEM}$ y el $\text{salt}$ para resistir ataques multi-target:
si el atacante ve múltiples ciphertexts bajo la misma clave pública, no puede
reutilizar trabajo entre distintos desafíos.

### 6.4 Desencapsulación — HQC-KEM.Decaps

```
Entradas: dk_KEM = (seed_KEM, ek_KEM), c = (u, v, salt)
Salidas: K  // clave compartida (o clave aleatoria si c es inválido)

1. Reconstruir (ek_PKE, dk_PKE) desde seed_KEM
2. m' = HQC-PKE.Decrypt(dk_PKE, c)
3. Reconstruir: (K', θ') = G(m', ek_KEM, salt)
4. c' = HQC-PKE.Encrypt(ek_KEM, m', θ')   // re-cifrado
5. Si c' == c:
      Output K'
   Else:
      Output J(seed_KEM, c)   // "implicit rejection": K pseudoaleatoria
```

El paso de **implicit rejection** (paso 5, rama else) es crítico para
IND-CCA2: si el ciphertext es inválido, en lugar de devolver error se devuelve
una clave pseudoaleatoria $J(\text{seed}_\text{KEM}, c)$. El atacante no puede
distinguir si el descifrado tuvo éxito o no — el oráculo de desencapsulación
siempre devuelve algo que parece una clave válida.

La re-cifrado del paso 4 garantiza que cualquier modificación del ciphertext
sea detectada.

### 6.5 El flujo completo en un diagrama

```
Sender (Encaps)                          Receiver (Decaps)
────────────────                         ─────────────────
m ←$ {0,1}^256
salt ←$ {0,1}^256
(K, θ) = SHAKE256(m ∥ ek ∥ salt)
c = PKE.Enc(ek, m, θ)
                          ── c ──▶
                                         m' = PKE.Dec(dk, c)
                                         (K', θ') = SHAKE256(m' ∥ ek ∥ salt)
                                         c'' = PKE.Enc(ek, m', θ')
                                         if c'' == c: output K'
                                         else:        output J(seed_KEM, c)
Output K
```

Si el canal es honesto: $m' = m$, luego $K' = K$. Las dos partes comparten la
misma clave $K$ sin haberla transmitido directamente.

---

## 7. Parámetros concretos

### 7.1 Tabla de parámetros (especificaciones 2025)

| Parámetro | HQC-128 | HQC-192 | HQC-256 |
|---|---|---|---|
| Nivel de seguridad NIST | 1 (128 bits) | 3 (192 bits) | 5 (256 bits) |
| $n$ (longitud bloque QC) | 17669 | 35851 | 57637 |
| $\omega$ (peso de x, y, r₁, r₂, e) | 66 | 100 | 131 |
| $k$ (longitud mensaje en bits) | 256 | 384 | 512 |
| $n_1$ (longitud bloque RM) | 256 ($2^8$) | 256 | 256 |
| $n_2$ (número de símbolos RS) | 90 | 120 | 150 |
| $t$ (capacidad corrección RS) | 57 | 75 | 99 |
| DFR objetivo | $\leq 2^{-128}$ | $\leq 2^{-192}$ | $\leq 2^{-256}$ |
| **Clave pública** $|\text{ek}|$ | **2249 bytes** | **4522 bytes** | **7245 bytes** |
| **Ciphertext** $|c|$ | **4497 bytes** | **9026 bytes** | **14469 bytes** |
| **Clave privada** $|\text{dk}|$ | **40 bytes** (seed) | **40 bytes** | **40 bytes** |

La clave privada es de solo 40 bytes porque se almacena únicamente la semilla
desde la que se regenera $(\mathbf{x}, \mathbf{y})$ cuando se necesita.

### 7.2 Comparativa con ML-KEM (Kyber)

| Métrica | ML-KEM-128 | HQC-128 |
|---|---|---|
| Clave pública | 800 bytes | 2249 bytes |
| Ciphertext | 768 bytes | 4497 bytes |
| Clave privada | 1632 bytes | 40 bytes |
| Seguridad base | LWE sobre retículos | SD sobre códigos QC |
| DFR | $2^{-164}$ aprox. | $2^{-128}$ |

HQC tiene ciphertexts más grandes que Kyber, pero su base de seguridad
es completamente ortogonal (teoría de códigos vs. retículos), lo que
justifica su selección como alternativa.

---

## 8. Análisis de seguridad

### 8.1 Reducción de seguridad formal

La seguridad de HQC-PKE se reduce al problema 2-DQCSD-P mediante una
reducción estándar tipo "distinguidor → decoder":

**Teorema (informal):** Si existe un adversario $\mathcal{A}$ que rompe
IND-CPA de HQC-PKE con ventaja $\epsilon$ y tiempo $t$, entonces existe un
algoritmo $\mathcal{B}$ que resuelve 2-DQCSD-P$(n, \omega)$ con ventaja
$\epsilon/2$ y tiempo similar a $t$.

La reducción funciona construyendo una instancia del problema QCSD como
clave pública y simulando el oráculo de cifrado usando la estructura del
problema.

La seguridad de HQC-KEM IND-CCA2 sigue del teorema de Fujisaki-Okamoto (FO)
aplicado a PKE IND-CPA, específicamente usando el análisis del transform
$\text{SFO}^\perp_m$ en el **ROM** (Random Oracle Model).

### 8.2 Ataques al syndrome decoding: ISD

El mejor ataque clásico conocido contra SD es **Information Set Decoding (ISD)**,
en sus variantes modernas:

- **Prange (1962):** baseline $O(2^{0.12n})$
- **Stern/Dumer (1989):** mejora mediante colisiones meet-in-the-middle
- **BJMM (2012):** Becker-Joux-May-Meurer, actualmente state-of-the-art para
  muchos parámetros
- **MMT (2011):** May-Meurer-Thomae
- **BCJ:** Becker-Coron-Joux (representación sobre enteros)

La complejidad del mejor ataque contra SD$(n, k, \omega)$ se estima con la
fórmula de **Finiasz-Sendrier**:

$$\log_2 W(n, k, \omega) \approx n \cdot h\left(\frac{\omega}{n}\right) - k$$

donde $h(p) = -p\log_2 p - (1-p)\log_2(1-p)$ es la entropía binaria.

Los parámetros de HQC están dimensionados para que el mejor ataque ISD
requiera $\geq 2^\lambda$ operaciones de bit para nivel de seguridad $\lambda$.

### 8.3 Efecto cuántico (Grover + BJMM)

El algoritmo de Grover proporciona una aceleración cuadrática para la búsqueda
no estructurada. Para ataques a SD:

- Un atacante cuántico puede usar **Grover + BJMM** para obtener una mejora
  cuadrática sobre el mejor atacante clásico.
- El efecto concreto sobre ISD es reducir la complejidad en un factor de
  $\approx \sqrt{\cdot}$ en la parte de búsqueda exhaustiva.

HQC compensa esto eligiendo $n$ y $\omega$ tales que incluso con la mejora
de Grover, la complejidad sea $\geq 2^{128}$ (para HQC-128), etc.

### 8.4 Ataques estructurales sobre códigos QC

La estructura quasi-cíclica podría en principio ser explotada. Los ataques
conocidos que intentan aprovecharse de QC:

- **Distinguish attacks:** intentan distinguir $(H, \mathbf{y})$ QC de uniforme.
  Los mejores tienen impacto sublineal en $n$ — despreciable para $n$ grande.
- **Algebraic attacks:** explotan la factorización de $X^n - 1$. Bloqueados
  por la condición de $n$ **primo primitivo** (solo dos factores irreducibles).
- **Parity attacks:** distinguidores basados en la paridad del síndrome.
  Bloqueados por la condición 2-DQCSD-**P** (con paridad).

La conclusión del análisis de los autores y la comunidad NIST es que, en la
práctica, el mejor ataque contra HQC es ISD sobre la instancia QC con un
factor adicional sublineal de ventaja por la estructura — el dimensionado de
parámetros absorbe este factor con margen.

### 8.5 Decryption Failure Rate (DFR) como vector de ataque

La DFR no es solo una métrica de corrección — es también una consideración
de seguridad. Si la DFR es demasiado alta, el atacante puede provocar fallos
controlados para extraer información sobre la clave secreta (**fault attacks**,
**DFR-based CCA**). Las especificaciones requieren que la DFR sea
$\leq 2^{-\lambda}$ para cada nivel $\lambda$.

El análisis de DFR en HQC requiere estudiar la distribución del peso de Hamming
del producto $\mathbf{x} \cdot \mathbf{r}_2 + \mathbf{r}_1 \cdot \mathbf{y} + \mathbf{e}$,
que se modela mediante aproximaciones combinatorias. La dificultad es que
$\mathbf{x} \cdot \mathbf{r}_2$ en $\mathcal{R}$ no tiene coordenadas independientes —
el modelo de Bernoulli independiente es una aproximación que sobreestima el
DFR (es conservadora).

---

## 9. Implicaciones para la implementación en Rust

### 9.1 Módulos principales

Una implementación Rust completa necesita los siguientes componentes:

```
hqc/
├── poly/              # Aritmética en R = F_2[X]/(X^n - 1)
│   ├── mul.rs         # Multiplicación QC (el hot path)
│   └── weight.rs      # Peso de Hamming, sampling
├── codes/
│   ├── reed_muller.rs  # RM(1,m), encoder + FHT decoder
│   └── reed_solomon.rs # RS[n2,k2,d2], encoder + BM decoder
├── pke.rs             # HQC-PKE: Keygen, Encrypt, Decrypt
├── kem.rs             # HQC-KEM: Keygen, Encaps, Decaps
├── hash.rs            # SHAKE256, SHA3-512 (wrappers)
└── params.rs          # HQC-128, HQC-192, HQC-256
```

### 9.2 La multiplicación de polinomios: el cuello de botella

La operación más costosa es la **multiplicación en $\mathcal{R}$**:
$\mathbf{u} \cdot \mathbf{v}$ con $n \approx 17000$–$57000$.

Opciones de implementación, en orden de eficiencia:

| Método | Complejidad | Notas |
|---|---|---|
| **Schoolbook** | $O(n^2)$ | Inaceptable para $n \geq 17000$ |
| **Karatsuba** | $O(n^{1.585})$ | Útil para tamaños intermedios |
| **NTT / FFT binaria** | $O(n \log n)$ | No directo en $\mathbb{F}_2[X]$ — requiere levantamiento |
| **Gao-Mateer (additive FFT)** | $O(n \log n)$ | FFT aditiva sobre $\mathbb{F}_{2^m}$, usada en la ref. C |
| **CLMUL / carry-less multiplication** | $O(n^2/64)$ con SIMD | x86: `pclmulqdq`; ARM: `vmull`. Vectorizable eficientemente |

La referencia C de HQC usa una combinación de **Karatsuba + SIMD word-level**
para la multiplicación de palabras de 64 bits, alcanzando complejidad efectiva
$\sim O(n^2/64)$ con buenas constantes.

Para una implementación Rust portátil (sin unsafe), una primera versión limpia
con Karatsuba word-level ya es competitiva. La ruta hacia SIMD usa
`std::arch` con target feature `avx2` / `aes` (para `pclmulqdq`).

**Nota sobre vectores sparse:** en la clave generación, $\mathbf{x}$ y $\mathbf{y}$
son de peso $\omega \approx 66$–$131$. La multiplicación por un vector sparse
puede implementarse como suma de $\omega$ rotaciones, con complejidad
$O(\omega \cdot n / 64)$ usando operaciones bitwise sobre `u64` arrays — mucho
más rápido que la multiplicación general para claves.

### 9.3 Fast Hadamard Transform

El decoder Reed-Muller requiere la WHT (Walsh-Hadamard Transform) sobre bloques
de $n_1 = 256$ bits. Esta operación es altamente paralelizable:

```
// WHT de longitud 256 con aritmética sobre i16
fn wht_256(input: [i16; 256]) -> [i16; 256] {
    let mut a = input;
    let mut h = 1usize;
    while h < 256 {
        for i in (0..256).step_by(h * 2) {
            for j in i..i+h {
                let x = a[j];
                let y = a[j + h];
                a[j]     = x + y;
                a[j + h] = x - y;
            }
        }
        h *= 2;
    }
    a
}
```

Los 256 símbolos RS se decodifican **independientemente** — paralelización
trivial con Rayon o SIMD sobre los 8 bloques de $n_1/\text{simd\_lanes}$ por pasada.

### 9.4 Reed-Solomon: el decoder algebraico

El decoder RS sobre $\mathbb{F}_{2^8}$ requiere:

1. Calcular los **síndromes** $S_j = \sum_i e_i \alpha^{ij}$ donde $\alpha$ es
   una raíz primitiva en $\mathbb{F}_{2^8}$.
2. Berlekamp-Massey (o Euclidean) para encontrar el **error locator polynomial**.
3. Chien search para encontrar las **posiciones de error**.
4. Forney algorithm para los **valores de error** (en $\mathbb{F}_{2^8}$ binario
   los valores son siempre 1, así que este paso es trivial).

Toda la aritmética es sobre $\mathbb{F}_{2^8}$ que se implementa eficientemente
con tablas de logaritmos discretos (tabla de antilogaritmos en base $\alpha$).

### 9.5 Constant-time: el requisito más delicado

HQC transmite información sensible a través de ramas condicionales si no se
implementa en tiempo constante. Los puntos críticos son:

- **Sampling de vectores de peso fijo:** `SampleFixedWeightVect$` usa rejection
  sampling y puede tener número de iteraciones variable — las iteraciones deben
  estar acotadas y la constante de tiempo mantenida mediante `cmov` / select
  constante.
- **Comparación de ciphertexts** en Decaps: usar `subtle::ct_eq` en lugar de `==`.
- **Implicit rejection:** la rama `if c' == c` debe evaluarse en tiempo constante
  y el `select` entre $K'$ y $J$ no debe filtrar en qué rama se tomó.
- **Reed-Solomon decoder:** el Chien search tiene número de raíces variable;
  la iteración debe recorrer siempre todos los $n_2$ candidatos.

El crate `subtle` de RustCrypto proporciona las primitivas necesarias
(`Choice`, `ConditionallySelectable`, `ConstantTimeEq`).

### 9.6 Estructura de tipos genéricos sugerida

```rust
/// Parámetros de HQC parametrizados por nivel de seguridad
pub trait HqcParams {
    const N:      usize;   // longitud bloque QC
    const OMEGA:  usize;   // peso de los vectores secretos
    const K:      usize;   // longitud mensaje en bits
    const N1:     usize;   // n1 = 2^m, tamaño bloque RM
    const N2:     usize;   // número de símbolos RS
    const T:      usize;   // capacidad corrección RS
    const LAMBDA: usize;   // nivel de seguridad en bits
}

pub struct Hqc128;
impl HqcParams for Hqc128 {
    const N: usize = 17669;
    const OMEGA: usize = 66;
    // ...
}

pub struct EncapsulationKey<P: HqcParams> {
    h: GfPoly<P>,   // polinomio h en R
    s: GfPoly<P>,   // polinomio s = x + h·y
}

pub struct DecapsulationKey<P: HqcParams> {
    seed: [u8; 40], // desde el que se regenera (x, y)
    ek:   EncapsulationKey<P>,
}
```

---

## 10. Apéndice

### 10.1 Glosario

| Término | Definición |
|---|---|
| **IND-CPA** | Indistinguishability under Chosen Plaintext Attack. El adversario no puede distinguir el cifrado de dos mensajes de su elección. |
| **IND-CCA2** | Indistinguishability under adaptive Chosen Ciphertext Attack. El adversario tiene acceso a un oráculo de descifrado para cualquier ciphertext excepto el desafío. |
| **SD / Syndrome Decoding** | Dado $(H, \mathbf{y})$, encontrar $\mathbf{x}$ sparse con $H\mathbf{x}^\top = \mathbf{y}^\top$. NP-completo. |
| **QCSD** | SD en el contexto de códigos quasi-cíclicos (la matriz $H$ tiene estructura circulante). |
| **ISD** | Information Set Decoding. Familia de ataques aleatorios contra SD basados en elegir subconjuntos de coordenadas. |
| **DFR** | Decryption Failure Rate. Probabilidad de que el descifrado produzca un error. |
| **Reed-Muller RM(1,m)** | Código lineal $[2^m, m+1, 2^{m-1}]$, decodificable con FHT. |
| **Reed-Solomon RS** | Código MDS $[n_2, k_2, n_2-k_2+1]$ sobre $\mathbb{F}_{2^m}$. Máxima distancia separable. |
| **FHT / WHT** | Fast (Walsh-)Hadamard Transform. Decodifica Reed-Muller en $O(n \log n)$. |
| **FO / HHK transform** | Fujisaki-Okamoto transform en su versión modular (Hofheinz-Hövelmanns-Kiltz). Convierte PKE IND-CPA en KEM IND-CCA2 en el ROM. |
| **ROM** | Random Oracle Model. Modelo donde las funciones hash se tratan como oráculo aleatorio ideal. |
| **Primitive prime** | Primo $n$ tal que $(X^n-1)/(X-1)$ es irreducible en $\mathbb{F}_2[X]$. |
| **Implicit rejection** | En desencapsulación, devolver un valor pseudoaleatorio en vez de error cuando el ciphertext es inválido. Crítico para IND-CCA2. |

### 10.2 Flujo de dependencias matemáticas

```
Primitive prime n
    │
    ▼
R = F_2[X]/(X^n - 1)          R_ω (peso fijo ω)
    │                              │
    ├──────────────────────────────┤
    ▼                              ▼
Quasi-Cyclic codes          Secret (x,y) ∈ R_ω × R_ω
    │                              │
    ▼                              ▼
2-DQCSD-P problem          Public key: s = x + h·y
    │                              │
    └──────────────────────────────┘
                │
                ▼
           HQC-PKE (IND-CPA)
                │
        RMRS codec (RM + RS)
                │
                ▼
         HHK/SFO⊥ transform
                │
                ▼
           HQC-KEM (IND-CCA2)
```

### 10.3 Referencias completas

1. **[AGM+2016/2018]** C. Aguilar-Melchor, O. Blazy, J.-C. Deneuville, P.
   Gaborit, G. Zémor. *Efficient Encryption from Random Quasi-Cyclic Codes.*
   IEEE Trans. Inf. Theory 64(5):3927–3943, 2018. arXiv:1612.05572.
   → **El paper original que define HQC.**

2. **[HQC-Spec-2025]** Aguilar-Melchor et al. *Hamming Quasi-Cyclic (HQC),
   specifications v22/08/2025.* https://pqc-hqc.org/doc/hqc_specifications_2025_08_22.pdf
   → **Las especificaciones oficiales actuales (NIST submission).**

3. **[HHK17]** D. Hofheinz, K. Hövelmanns, E. Kiltz. *A Modular Analysis of
   the Fujisaki-Okamoto Transformation.* TCC 2017, LNCS 10677.
   → El análisis formal del transform FO que usa HQC.

4. **[BMW+12]** A. Becker, A. Joux, A. May, A. Meurer. *Decoding Random
   Binary Linear Codes in $2^{n/20}$: How $1+1=0$ Improves ISD.*
   EUROCRYPT 2012. → El mejor ataque ISD clásico (BJMM).

5. **[FS09]** M. Finiasz, N. Sendrier. *Security Bounds for the Design of
   Code-Based Cryptosystems.* ASIACRYPT 2009. → Estimación de seguridad para
   SD, usada en el dimensionado de parámetros HQC.

6. **[GM10]** S. Gao, T. Mateer. *Additive Fast Fourier Transforms over Finite
   Fields.* IEEE Trans. Inf. Theory 56(12), 2010. → La FFT aditiva binaria
   usada para la multiplicación de polinomios eficiente.