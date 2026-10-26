# MATLAB Codes for the Myers–Perry–AdS Gravitational Atom

## Overview

This repository contains MATLAB scripts supporting the analytical investigation of massive scalar quasibound states around small, singly spinning Myers–Perry–anti-de Sitter (MP–AdS) black holes.

The associated manuscript is entitled:

**The Singly Spinning Myers–Perry–AdS Gravitational Atom: Fine Structure, Rotational Hyperfine Splitting, and Superradiant Widths**

The study examines how a black hole modifies the discrete scalar spectrum of global AdS through conservative gravitational corrections, rotational frame dragging, and horizon-induced damping or superradiant growth.

The MATLAB scripts evaluate analytical expressions derived in the manuscript and reproduce its three quantitative figures. They do not implement an independent numerical solution of the full quasibound-state eigenvalue problem.

## Repository Contents

| File | Description |
|---|---|
| `PotentialFigure_MPAdS.m` | Evaluates and plots the effective radial potentials for two representative scalar modes, including their trapping regions, turning points, and the asymptotic AdS wall. |
| `LevelDiagram_MPAdS.m` | Reproduces the analytical fine-structure and rotational hyperfine splitting of the scalar frequency levels. |
| `GrowthDiagram_MPAdS.m` | Reproduces the magnitudes and signs of the horizon-induced damping and superradiant growth rates. |

## Physical Model and Parameters

The calculations concern a massive test scalar field in a singly spinning MP–AdS background in $D=4+n$ spacetime dimensions.

The principal illustrative example uses:

| Parameter | Value |
|---|---|
| Spacetime dimension $D$ | 7 |
| AdS curvature radius $\ell$ | 50 |
| Scalar mass $\mu_s$ | 0.1 |
| Black-hole mass parameter $r_m^4$ | 0.810513 |
| Oscillator shell $\mathcal K=2N+l$ | 2 |

**The rotation parameter depends on the figure:**

| Script | Rotation parameter | Purpose |
|---|---|---|
| `PotentialFigure_MPAdS.m` | $a=0$ | Static-background effective potentials and trapping regions |
| `LevelDiagram_MPAdS.m` | $a=0.1$ | Conservative fine structure and rotational hyperfine splitting |
| `GrowthDiagram_MPAdS.m` | $a=0.1$ | Horizon-induced decay and superradiant growth |

The representative scalar modes belong to the shell $\mathcal K=2$:

- $(N,l)=(1,0)$
- $(N,l)=(0,2)$, with azimuthal sublevels $m=-2,-1,0,+1,+2$

All calculations use geometrized units, $G=c=1$, with the normalization conventions specified in the manuscript.

## Analytical Methodology

The computational workflow is based on:

- The exact global-AdS scalar normal-mode spectrum.
- Analytical conservative frequency corrections induced by the black-hole gravitational field.
- Leading rotational hyperfine splitting associated with frame dragging.
- Matched-asymptotic expressions for horizon-induced damping and superradiant growth.
- Analytical effective-potential expressions, followed by numerical evaluation of radial profiles and turning points.

The scripts evaluate these expressions for the specified parameter choices and construct the corresponding figures.

The results are valid within the approximation regimes discussed in the manuscript, particularly the small-black-hole and slow-rotation limits.

## Requirements

- MATLAB with support for `tiledlayout` and `exportgraphics` (MATLAB R2020a or later).
- No independent numerical quasibound-state eigensolver is required.

## Usage

Download or clone the repository, open MATLAB, and set the current working directory to the repository folder.

Run the scripts separately:

```matlab
PotentialFigure_MPAdS
LevelDiagram_MPAdS
GrowthDiagram_MPAdS
```

Each script generates its corresponding scientific figure. PDF and PNG exports, together with numerical CSV data where implemented, are generated according to the output settings defined in the script.

The scripts may overwrite previously generated output files with the same names.

## Reproducibility and Scope

The repository provides the MATLAB implementation used to evaluate the analytical formulas and generate the quantitative figures accompanying the manuscript.

The effective-potential figure uses a nonrotating background, whereas the spectral-level and growth/decay figures incorporate the leading rotational corrections.

Numerical root finding, quadrature, and plotting operations are used only to evaluate analytical expressions and construct figures. They should not be interpreted as an independent numerical determination of the full MP–AdS quasibound-state spectrum.

The reported numerical values must be understood within the accuracy and parameter regimes of the underlying analytical approximations.

## Citation

If you refer to this software in academic work, please cite the associated research article and the archived software release once their bibliographic details and DOI become available.

Citation metadata are provided in `CITATION.cff`.

## Copyright and Usage

Copyright © 2026 Hamed Ramezani.

This repository is made available under the Creative Commons Attribution 4.0 International (CC BY 4.0) license.

License: https://creativecommons.org/licenses/by/4.0/

When using these scripts in academic research, please also cite the associated software release and research article.
