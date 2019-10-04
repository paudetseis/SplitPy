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

## Installation

### Dependencies

The current version was developed using **Python3.7** \
Also, the following packages are required:

- [`obspy`](https://github.com/obspy/obspy/wiki)
- [`stdb`](https://github.com/schaefferaj/StDb)
- [`dill`](https://pypi.org/project/dill/)
- [`PyQt5`](https://pypi.org/project/PyQt5/)

#### Conda environment

We recommend creating a custom 
[conda environment](https://conda.io/docs/user-guide/tasks/manage-environments.html)
where `SplitPy` can be installed along with its dependencies.

```bash
   conda create -n split python=3.7 obspy dill -c conda-forge
```

Activate the newly created environment:

```bash
   conda activate split
```

Install remaining dependencies using `pip` inside the `split` environment. 
Note that you need to install `StDb` from source at this time:

```bash
   git clone https://github.com/paudetseis/StDb.git
   cd StDb
   pip install .
   pip install PyQt5
```

### Installing from source

Download or clone the repository:
```bash
git clone https://github.com/paudetseis/SplitPy.git
cd SplitPy
```

Next we recommend following the steps for creating a `conda` environment 
(see [above](#conda-environment)). Then install using `pip`:

```bash
pip install .
``` 

## Usage 

### Documentation & Tutorials

The documentation for all classes and functions in `splitpy` and some tutorials
on how to use the python scripts bundled with the package can be accessed 
from https://paudetseis.github.io/SplitPy/.

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
