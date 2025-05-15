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

The current version has been tested using **Python = 3.12** \
Also, the following package is required:

- `stdb <https://github.com/paudetseis/StDb>`_

Other required packages (e.g., ``obspy``, ``PyQt5``)
will be automatically installed by ``stdb``.

Conda environment
-----------------

We recommend creating a custom ``conda`` environment
where ``SplitPy`` can be installed along with some of its dependencies.

.. sourcecode:: bash

   conda create -n split -c conda-forge python=3.12 obspy

Activate the newly created environment:

.. sourcecode:: bash

   conda activate split

Install remaining dependencies using ``pip`` inside the ``split`` environment:

.. sourcecode:: bash

   pip install git+https://github.com/schaefferaj/stdb

Installing from GitHub development branch
-----------------------------------------

.. sourcecode:: bash

   pip install git+https://github.com/paudetseis/splitpy

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

Using local data
================

The main script packaged with ``RfPy`` uses FDSN web services through and ``ObsPy`` `Client` to load waveform data. For waveform data locally stored on your hard drive, the scripts can use a `Client` that reads a `SeisComP Data Structure <https://docs.obspy.org/packages/autogen/obspy.clients.filesystem.sds.html>`_ archive containing SAC or miniSEED waveform data. Check out the scripts ``rfpy_calc`` below and the argument ``--local-data`` and ``--dtype`` for more details.

Station Metadata
----------------

If you have data stored locally on your drive, it is likely you also have a station `XML <https://www.fdsn.org/xml/station/>`_ file containing the metadata. The corresponding ObsPy documentation is `here <https://docs.obspy.org/packages/obspy.core.inventory.html>`_. 

To convert the station `XML` file to an input that can be read by ``OrientPy``, you run the command ``gen_stdb station.xml`` (only available on StDb version 0.2.7), which will create the file ``station.pkl``. If you don't have a station `XML` file but you have a dataless SEED file, you can convert it first to `XML` using `this tools <https://seiscode.iris.washington.edu/projects/stationxml-converter>`_.

Waveform Data
-------------

The SDS folder containing the waveform data has the structure:

.. code-block:: python

   archive
     + year
       + network code
         + station code
           + channel code + type
             + one file per day and location, e.g. NET.STA.LOC.CHAN.TYPE.YEAR.DOY


For example:

.. code-block:: python

   SDS/
     2020/
       NY/
         TGTN/
           HHZ.D/ 
             NY.TGTN..HHZ.D.2020.332
             ...


Note, the filename does not include the extension (`.MSEED` or `.SAC`), and the characters `.D` (for type Data) that appear in both the channel code and the filename. Note also the two dots (`..`). If there is a location code, it should appear between those dots (e.g., for a location code `10`, the corresponding filename should be `NY.TGTN.10.HHZ.D.2020.332`). There is no location code for the NY.TGTN data, and this field is simply absent from the filenames. Finally, the day-of-year (DOY) field must be zero-padded to be exactly 3 characters.
