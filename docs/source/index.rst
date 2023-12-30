RetroBioCat 2.0
===============

RetroBioCat 2.0 (rbc2) is a python package for computer-aided synthesis planning
in biocatalysis.  Version 2 builds on the original RetroBioCat package, creating
a more modular and extensible framework allowing for the incorporation of other
synthesis planning approaches which can be combined together.  This allows for hybrid
synthsis planning, for example combining chemistry and biocatalysis.

Unlike the original RetroBioCat package, rbc2 does not feature a web-app element,
and is instead designed to be used as soley a python package.  This allows for more
flexibility in how the package is used.

Installation
------------
It is highly recommended to install this package in its own environment, for example a conda environment.
Python >= 3.9 is required.

::

    pip install rbc2


**A note on M1 Macs:**  There seems to be an issue installing tables on M1 macs, to do with location HDF5.
Assuming you have installed the necessary files using homebrew, inside your environment run the following before installation:

::

    export HDF5_DIR=/opt/homebrew/opt/hdf5
    export BLOSC_DIR=/opt/homebrew/opt/c-blosc


Testing
-------


Usage
-----












