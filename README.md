# MATLAB Codes for the Myers–Perry–AdS Gravitational Atom

## Overview

This repository contains MATLAB scripts supporting the analytical study of massive scalar quasibound states around small, singly spinning Myers–Perry–anti-de Sitter (MP–AdS) black holes.

The associated manuscript is entitled:

**The Singly Spinning Myers–Perry–AdS Gravitational Atom: Fine Structure, Rotational Hyperfine Splitting, and Superradiant Widths**

The analysis investigates how a black hole modifies the discrete scalar spectrum of global AdS through conservative gravitational corrections, rotational frame dragging, and horizon-induced damping or superradiant growth.

The MATLAB scripts evaluate analytical expressions derived in the manuscript and reproduce its quantitative figures. They do not implement an independent numerical solution of the full quasibound-state eigenvalue problem.

## Repository Contents

| File | Description |
|---|---|
| `PotentialFigure_MPAdS.m` | Computes and plots the effective radial potentials for the two representative modes, including the trapping regions and turning points. |
| `LevelDiagram_MPAdS.m` | Reproduces the analytical fine-structure and rotational hyperfine level diagram. |
| `GrowthDiagram_MPAdS.m` | Plots the horizon-induced damping and superradiant growth rates for the representative scalar modes. |

## Physical Model

The calculations concern a massive test scalar field propagating on a singly spinning MP–AdS black-hole background in \(D=4+n\) spacetime dimensions.

The principal illustrative example uses:

| Parameter | Value |
|---|---|
| Spacetime dimension \(D\) | 7 |
| AdS curvature radius \(\ell\) | 50 |
| Scalar mass \(\mu_s\) | 0.1 |
| Rotation parameter \(a\) | 0.1 |
| Black-hole mass parameter \(r_m^4\) | 0.810513 |
| Oscillator shell \(\mathcal K=2N+l\) | 2 |

The representative scalar modes are \((N,l)=(1,0)\) and \((N,l)=(0,2)\). The latter includes rotational sublevels labeled by the azimuthal number \(m\).

All calculations use geometrized units, \(G=c=1\), with the normalization conventions specified in the manuscript.

## Methodology

The computational workflow is based on:

- The exact global-AdS scalar normal-mode spectrum.
- Analytical conservative frequency shifts induced by the black-hole geometry.
- Leading rotational splitting associated with frame dragging.
- Matched-asymptotic expressions for horizon-induced damping and superradiant growth.
- Analytical effective-potential expressions and numerical evaluation of their turning points.

The MATLAB scripts evaluate these analytical results for the specified parameter choices and generate the corresponding figures.

The validity of the results is restricted to the approximation regimes discussed in the manuscript, including the small-black-hole and slow-rotation assumptions.

## Requirements

- MATLAB with support for local functions in scripts.
- No independent numerical quasibound-state solver is required.

## Usage

Download or clone the repository, open MATLAB, and set the current folder to the repository directory.

Run the scripts separately:

```matlab
PotentialFigure_MPAdS
LevelDiagram_MPAdS
GrowthDiagram_MPAdS
```

Each script produces the corresponding figure. Output file names and export settings are defined within the scripts.

## Reproducibility

The repository is intended to support reproducibility of the analytical plots and illustrative numerical evaluations presented in the manuscript.

The reported numerical values should be interpreted within the order of the analytical approximations, rather than as independent high-precision numerical quasibound-state frequencies.

## Citation

If you use these scripts, please cite the associated research article and the archived software release once their bibliographic information and DOI become available.

## License

License information will be provided with the first archived software release.
