
.. pymet documentation master file, created by
   sphinx-quickstart on Wed Mar 20 22:01:03 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pymet's documentation!
=================================

This is the documentation for pymet. pymet is Python based Meteorological Analysis and Graphics module.

Installing
===========

These are external packages which you will need to install before using pymet.

- NumPy <http://www.numpy.org/>
- SciPy <http://www.scipy.org/>

and module :mod:`pymet.metplt` needs to install Matplotlib and Basemap Matplotlib Toolkit,

- Matplotlib <http://matplotlib.org/>
- Basemap Matplotlib Toolkit <http://matplotlib.org/basemap/>

and module  :mod:`pymet.io`  needs to install python-grads,

- python-grads <http://sourceforge.net/projects/opengrads/files/python-grads/>

Before install python-grads, you need to patch gacore.py in python-grads with patches in patch directory. ::

  $ tar xvf pygrads-1.1.8.tar.gz
  $ cd pygrads-1.1.8/grads
  $ cp <pymet base directory>/patch/pygrads-1.1.8_gacore.py.patch ./
  $ patch -c < pygrads-1.1.8_gacore.py.patch

After install external packages, you will need to compile internal module by f2py, ::

  $ cd pymet
  $ f2py -c -m _internal _internal.f90


References
===========

.. toctree::
   :maxdepth: 1
   :titlesonly:
  
   pymet.constants
   pymet.dynamics
   pymet.grid
   pymet.stats
   pymet.tools
   pymet.io
   pymet.field
   pymet.metplt

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

