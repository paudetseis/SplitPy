# SplitPy: Software for teleseismic shear-wave splitting analysis

![](./splitpy/examples/figures/SplitPy_logo.png)

Seismic anisotropy refers to the property of seismic waves to propagate
at different wavespeeds depending on the direction of propagation. This
property can be related to the coherent alignment of rock-forming minerals,
which is thought to reflect the current dynamics or fossilized structure of Earth
materials due to tectonic deformation. In the upper mantle, seismic anisotropy 
is most easily determined using the distortion of teleseismic body waves with a 
known initial polarity, typically core-refracted shear-waves (SKS, SKKS). 

SplitPy is a teleseismic shear-wave Splitting Toolbox based on the 
Matlab Tool [`SplitLab`](http://splitting.gm.univ-montp2.fr), 
but with modifications from [Wustefeld et al (2008)](#references). 
Additional error surface implementation has been added, however these error 
surfaces have not been fully tested. The code produces output identical to
those in [Audet et al. (2016)](#references)

[![DOI](https://zenodo.org/badge/211722700.svg)](https://zenodo.org/badge/latestdoi/211722700)
[![Build Status](https://travis-ci.org/paudetseis/SplitPy.svg?branch=master)](https://travis-ci.org/paudetseis/SplitPy)

## Citing

If you use `SplitPy` in your work, please cite the Zenodo DOI above.

## Documentation 
Installation, Usage, API documentation and tutorials are described at 
https://paudetseis.github.io/SplitPy/.

## References

- Audet, P., Sole, C., and Schaeffer, A.J. (2016). Control of lithospheric
  inheritance on neotectonic activity in northwestern Canada? Geology,
  44, 807-810, https://doi.org/10.1130/G38118.1

- Wustefeld, A., and Bokelmann, G. (2007). Null detection in shear-wave splitting 
  measurements. Bulletin of the Seismological Society of America, 97, 1204-1211,
  https://doi.org/10.1785/0120060190

- Wustefeld, A., Bokelmann, G., Zeroli, C., and Barruol, G. (2008). SplitLab: 
  A shear-wave splitting environment in Matlab. Computers & Geoscience, 34, 
  515-528, https://doi.org/10.1016/j.cageo.2007.08.002
