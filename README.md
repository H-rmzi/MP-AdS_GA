
# MATLAB Codes for the Myers–Perry–AdS Gravitational Atom

## Overview

This repository contains MATLAB codes supporting the analytical investigation of massive scalar quasibound states around small, singly spinning Myers–Perry–anti-de Sitter (MP–AdS) black holes.

The associated manuscript is entitled:

**The Singly Spinning Myers–Perry–AdS Gravitational Atom: Fine Structure, Rotational Hyperfine Splitting, and Superradiant Widths**

The study examines how the black-hole gravitational field modifies the discrete global-AdS scalar spectrum through conservative shifts, rotational frame dragging, and horizon-induced damping or superradiant growth. The *gravitational atom* terminology describes an AdS-confined, atom-like spectral organization; it does not imply an exact hydrogenic spectrum or an exact relativistic harmonic oscillator.

The repository contains three MATLAB scripts reproducing the manuscript's quantitative figures and **one self-contained MATLAB file** reproducing selected analytical coefficients and an independent radial-node consistency check. These programs evaluate analytical approximations; they do not constitute a numerical solver for the full rotating quasibound-state eigenvalue problem.

## Repository Contents

| File | Description |
|---|---|
| `PotentialFigure_MPAdS.m` | Evaluates and plots effective radial potentials for two representative scalar modes, including trapping regions, turning points, and the asymptotic AdS wall. |
| `LevelDiagram_MPAdS.m` | Reproduces analytical fine-structure and leading rotational hyperfine splitting of the scalar frequency levels. |
| `GrowthDiagram_MPAdS.m` | Reproduces the magnitude and sign of horizon-induced damping and superradiant growth rates. |
| `MPAdS_Analytical_Supplement.m` | **Single-file analytical supplement:** evaluates the static resonant coefficient `Q_out`, the mixed rotational–resonant coefficients `A_N` and `B_N` for `N=0,1`, and an optional independent node-free radial check. All helper functions are included in this file. |

The three figure scripts remain separate and unchanged. Only **one MATLAB source file** is required for the additional analytical-coefficient calculations.

## Physical Model and Parameters

The calculations concern a massive test scalar field in a singly spinning MP–AdS background with spacetime dimension \(D=4+n\). All dimensional quantities use geometrized units \(G=c=1\), with normalization conventions defined in the manuscript.

The principal illustrative parameters are:

| Parameter | Value |
|---|---|
| Spacetime dimension \(D\) | 7 |
| AdS curvature radius \(\ell\) | 50 |
| Scalar mass \(\mu_s\) | 0.1 |
| Black-hole mass parameter \(h=r_m^4\) | 0.810513 |
| Dimensionless scalar mass \(\beta=\mu_s\ell\) | 5 |
| Representative figure shell \(\mathcal K=2N+l\) | 2 |

The figure scripts use different rotation parameters:

| Script | Rotation parameter | Purpose |
|---|---|---|
| `PotentialFigure_MPAdS.m` | \(a=0\) | Static effective potentials and trapping regions |
| `LevelDiagram_MPAdS.m` | \(a=0.1\) | Conservative fine structure and leading rotational hyperfine splitting |
| `GrowthDiagram_MPAdS.m` | \(a=0.1\) | Leading horizon-induced damping and superradiant growth |

The representative modes displayed in the figures belong to the \(\mathcal K=2\) shell: \((N,l)=(1,0)\) and \((0,2)\), with \(m=-2,-1,0,1,2\) for the latter family.

The **analytical supplement has a different, deliberately restricted mode selection**: \(D=7\), \(l=2\), and \(N=0,1\). These modes belong to shells \(\mathcal K=2\) and \(\mathcal K=4\), respectively. The standard reference calculation uses \(\beta=5\).

## Analytical Methodology

The MATLAB computations use the exact global-AdS scalar normal-mode functions and frequencies, conservative gravitational perturbations, leading frame-dragging corrections, and matched-asymptotic expressions for the effect of the black-hole horizon.

The additional analytical supplement evaluates the finite-part coefficient \(\mathcal Q_{\mathrm{out}}\) for the nonrotating resonant mode \((N,l)=(0,2)\) and the coefficients controlling the part of the **real frequency that is odd in azimuthal number \(m\)** at linear order in rotation.

For the latter calculation, the conventions are

\[
H=\frac{h}{\ell^4},\qquad
\xi=\frac{ma}{\ell},\qquad
\Omega=\ell\omega^{\mathrm{stat}},
\]

with the formal small-black-hole expansion

\[
\left.\frac{\partial\operatorname{Re}\Omega}{\partial\xi}\right|_{\xi=0}
=H\Lambda_N+H^2\left[A_N\ln\left(\frac{r_+^{(0)}}{\ell}\right)+B_N\right]+\cdots.
\]

Here \(r_+^{(0)}\) is the horizon radius of the nonrotating background. The finite coefficient \(B_N\) is defined **together with the displayed logarithmic convention**; only the combined logarithmic and finite contribution has invariant physical meaning under changes of logarithmic scale.

The optional *node-free* check solves the first-order inhomogeneous radial problem without dividing by the square of the reference radial eigenfunction at its internal node for \(N=1\). It tests the finite-part evaluation used for \(B_1\).

## Requirements

- **MATLAB R2025b:** the author-reported successful execution environment for the single-file analytical supplement, including both `PASS` tests.
- No Symbolic Math Toolbox is needed for `MPAdS_Analytical_Supplement.m`.
- The original figure scripts use MATLAB functionality including `tiledlayout` and `exportgraphics` (MATLAB R2020a or later); the complete four-file collection has not been separately certified on all earlier releases.
- No independent numerical quasibound-state eigensolver is required.

## Usage

Clone or download the repository, open MATLAB, and set the current working directory to the folder containing the `.m` files.

### Reproduce the manuscript figures

Run the existing scripts separately:

```matlab
PotentialFigure_MPAdS
LevelDiagram_MPAdS
GrowthDiagram_MPAdS
```

PDF/PNG exports and CSV data, where implemented, follow the settings in the individual figure scripts. Previously generated output files with the same names may be overwritten.

### Reproduce the analytical coefficients

Run the **single** supplemental MATLAB file:

```matlab
MPAdS_Analytical_Supplement
```

This runs the analytical coefficient calculations **and** the independent node-free check. For a faster run without the node-free check:

```matlab
MPAdS_Analytical_Supplement(false)
```

The reference parameter can also be supplied explicitly:

```matlab
results = MPAdS_Analytical_Supplement(true, 5);
```

The source supports a nonnegative value of \(\beta\), but the bundled reference assertions and reported validation below are specifically for \(\beta=5\).

The supplement writes output files to a `results/` folder beside the MATLAB source file, including `analytical_coefficients.csv`, `run_summary.txt`, and (when the node-free check is enabled) `nodefree_convergence_matlab.csv`. Existing files with the same names may be overwritten. The function also returns a MATLAB structure with the results.

## Reference Results and Validation

For \(D=7\), \(l=2\), \(\beta=5\):

| Quantity | \(N=0\) | \(N=1\) |
|---|---:|---:|
| \(\Omega_0\) | 10.831 | 12.831 |
| First-order static coefficient \(\Omega_1\) | -19.683 | -62.134 |
| Leading rotational coefficient \(\Lambda_N\) | 19.683 | 77.880 |
| AdS radial norm | \(2.386\times10^{-4}\) | \(2.540\times10^{-5}\) |
| Far-zone finite part | 0.3042 | 1.851 |
| Logarithmic coefficient \(A_N\) | -2096 | \(-1.968\times10^4\) |
| Finite coefficient \(B_N\) | 815.6 | 5374 |

The static coefficient is reproduced through both the closed expression and a beta-integral reconstruction:

```text
Q_out (closed)    = 5.811149809582
Q_out (beta form) = 5.811149809582
```

For the independent node-free check at \(\epsilon=0.01\), the MATLAB R2025b output reported:

```text
Node-free I_FP = 1.8507353670
Reference I_FP = 1.8506888890

PASS: independent node-free finite-part check.
PASS: static and mixed-coefficient reference tests.
```

The node-free value differs from the reference by approximately \(4.648\times10^{-5}\) at that finite cutoff, within the specified convergence tolerance. The author-reported console output confirms successful execution of the unified file in MATLAB R2025b; it does not establish that the full rotating QBS equations have been independently solved in MATLAB.

## Reproducibility and Scope

The codes reproduce selected results derived within the **small-black-hole, low-frequency, and slow-rotation** approximations of the manuscript. The mixed rotational–resonant expansion is taken at linear order in \(a\), with the derivative evaluated at \(a=0\) before studying the small-hole limit.

Numerical quadrature, root finding, and plotting in the scripts serve to evaluate analytical formulas or perform consistency checks. They are **not** independent numerical determinations of the full rotating MP–AdS quasibound-state spectrum. The code does not prove a uniform bound on all neglected asymptotic terms or supply exact finite-spin QBS frequencies.

## Citation and Archived Version

The manuscript's existing Zenodo citation refers to the **previously archived** software version:

**Software:** *Analytical MATLAB Codes for Fine Structure, Rotational Hyperfine Splitting, and Superradiant Widths in Myers-Perry-AdS*

**Archived version:** `v1.0.1`

**DOI:** https://doi.org/10.5281/zenodo.23676729

This DOI identifies the prior release. **Do not assume that version `v1.0.1` includes `MPAdS_Analytical_Supplement.m`.** After the updated repository is archived as a new version within the existing Zenodo record, its version-specific citation and, where appropriate, the manuscript's data-availability statement should be updated to identify the version that actually contains the analytical supplement.

If these codes are used in academic research, please cite the relevant archived software version and the associated research article when available.

## Copyright and Usage

Copyright © 2026 Hamed Ramezani.

This repository is made available under the **Creative Commons Attribution 4.0 International (CC BY 4.0)** license.

License: https://creativecommons.org/licenses/by/4.0/

When using these scripts in academic research, please cite the corresponding archived software release and associated article.
