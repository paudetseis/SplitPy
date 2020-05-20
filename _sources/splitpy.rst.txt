.. figure:: ../splitpy/examples/figures/SplitPy_logo.png
   :align: center

Licence
=======

Copyright 2019 Pascal Audet & Andrew Schaeffer

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Installation
============

Dependencies
------------

The current version has been tested using **Python > 3.6** \
Also, the following package is required:

- `stdb <https://github.com/paudetseis/StDb>`_

Other required packages (e.g., ``obspy``, ``PyQt5``)
will be automatically installed by ``stdb``.

Conda environment
-----------------

We recommend creating a custom ``conda`` environment
where ``SplitPy`` can be installed along with some of its dependencies.

.. sourcecode:: bash

   conda create -n split python=3.7 obspy -c conda-forge

Activate the newly created environment:

.. sourcecode:: bash

   conda activate split

Install remaining dependencies using ``pip`` inside the ``split`` environment:

.. sourcecode:: bash

   pip install stdb

Installing from Pypi
--------------------

*This option is not available at this time*

Installing from source
----------------------

- Clone the repository:

.. sourcecode:: bash

   git clone https://github.com/paudetseis/SplitPy.git
   cd SplitPy

- Install using pip:

.. sourcecode:: bash

   pip install .

Citing SplitPy
--------------

If you use SplitPy in your work, please cite the Zenodo DOI: https://doi.org/10.5281/zenodo.3564780

"""

.. automodule:: splitpy
   :members:

Basic Usage
===========

