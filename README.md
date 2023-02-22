# construction

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/jasenfinch/construction/workflows/R-CMD-check/badge.svg)](https://github.com/jasenfinch/construction/actions)
[![Codecov test coverage](https://codecov.io/gh/jasenfinch/construction/branch/master/graph/badge.svg)](https://codecov.io/gh/jasenfinch/construction?branch=master)
[![GitHub release](https://img.shields.io/github/release/jasenfinch/construction.svg)](https://GitHub.com/jasenfinch/construction/releases/)
<!-- badges: end -->

> Consensus Structural Chemical Classifications For Putative Molecular Formula Assignments of *m/z* Features from ESI-MS

## Overview

This R package provides an approach for compiling consensus structural chemical classifications of putative molecular formulas assigned to *m/z* features from electrospray ionisation mass spectrometry (ESI-MS) approaches. 
Putative ionisation products can be searched against the [KEGG compound database](https://www.genome.jp/kegg/compound/) and/or [PubChem](https://pubchem.ncbi.nlm.nih.gov/) to identify candidate compounds.
These candidate compounds are then searched against the [ClassyFire](http://classyfire.wishartlab.com/) database from which a consensus structural chemical classification is calculated.

## Installation

The `construction` package can be installed from GitHub using the
following:

```
remotes::install_github('jasenfinch/construction')
```

## Learn more

The package documentation can be browsed online at
<https://jasenfinch.github.io/construction/>.

If you believe you’ve found a bug in `construction`, please file a bug (and, if possible, a [reproducible example](https://reprex.tidyverse.org)) at
<https://github.com/jasenfinch/construction/issues>.
