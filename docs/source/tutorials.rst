*********
Tutorials
*********

Getting started
===============

For a simple walk-through of using |Triumvirate|, see `'Quick Guide'
<tutorials/QuickGuide.html>`_.

To learn about program parameter configuration,
see `'Parameter Configuration' <tutorials/Configuration.html>`_.

.. toctree::
    :hidden:
    :maxdepth: 1

    tutorials/QuickGuide.ipynb
    tutorials/Configuration.ipynb


Usage guides
============

For advanced usage of the Python package, consult the following tutorials
for various components:

- `Customised Logger <tutorials/Logger.html>`_
- `Parameter Set <tutorials/Parameters.html>`_
- `Binning Scheme <tutorials/Binning.html>`_
- `Particle Catalogue <tutorials/Catalogue.html>`_
- `Clustering Measurements <tutorials/Measurements.html>`_

.. toctree::
    :hidden:
    :maxdepth: 1

    tutorials/Logger.ipynb
    tutorials/Parameters.ipynb
    tutorials/Binning.ipynb
    tutorials/Catalogue.ipynb
    tutorials/Measurements.ipynb


Although the C++ library is intended as the backend for the Python package,
a 'black-box' C++ program is provided to offer users an end-to-end clustering
measurement pipeline. The guide to running the program can be found in
`'Running the C++ program' <tutorials/QuickGuide.html#running-the-c-program>`_.

.. admonition:: OpenMP parallelisation

    For both Python and C++ programs built with OpenMP enabled
    (see ':doc:`installation`'), you can set the environmental variable
    ``OMP_NUM_THREADS`` for multithreaded execution:

    .. code-block:: console

        $ export OMP_NUM_THREADS=<num-threads>


For developers, ':doc:`apidoc_cpp/apidoc_cpp`' contains the full C++
API reference. To reuse C++ routines, ensure the linker finds ``libtrv``,
e.g. by providing the path to the static library (built following
':doc:`installation`') with the flag :code:`-L<path-to-libtrv.a>`
during compilation.


.. |Triumvirate| raw:: html

    <span style="font-variant: small-caps">Triumvirate</span>
