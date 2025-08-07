*********
Tutorials
*********

Getting started
===============

For a simple walk-through of using |Triumvirate|, see `Quick Guide
<tutorials/QuickGuide.html>`_.

To learn about program parameter configuration,
see `Parameter Configuration <tutorials/Configuration.html>`_.

.. toctree::
    :hidden:
    :maxdepth: 1

    tutorials/QuickGuide.ipynb
    tutorials/Configuration.ipynb


Usage guides
============

For advanced usage of the Python package, consult the following sections
for various functionalities:

- `Customised Logger <tutorials/Logger.html>`_
- `Parameter Set <tutorials/Parameters.html>`_
- `Particle Catalogue <tutorials/Catalogue.html>`_
- Clustering Measurements

  - `Binning Scheme <tutorials/measurement/Binning.html>`_
  - `Clustering Statistics <tutorials/measurement/Statistics.html>`_

- Window Convolution

  - `Window Function Multipoles <tutorials/window/WindowFunctionMultipoles.html>`_
  - `Window Convolution Series <tutorials/window/WindowConvolutionSeries.html>`_
  - `Window Convolution Matrix <tutorials/window/WindowConvolutionMatrix.html>`_
  - `Window Convolution Example 1 <tutorials/window/WindowConvolution-QSO_complete.html>`_
  - `Window Convolution Example 2 <tutorials/window/WindowConvolution-LRG_FA.html>`_


.. toctree::
    :hidden:
    :maxdepth: 1

    tutorials/Logger.ipynb
    tutorials/Parameters.ipynb
    tutorials/Catalogue.ipynb
    tutorials_measurement
    tutorials_window


Although the C++ library is intended as the backend for the Python package,
a 'black-box' C++ program is provided to offer users an end-to-end clustering
measurement pipeline. The guide to running the program can be found in
`'Running the C++ program' <tutorials/QuickGuide.html#running-the-c-program>`_.

.. admonition:: OpenMP parallelisation

    For both Python and C++ programs built with OpenMP enabled
    (see :doc:`installation`), you can set the environmental variable
    ``OMP_NUM_THREADS`` for multithreaded execution:

    .. code-block:: sh

        export OMP_NUM_THREADS=<num-threads>

    In general, the following process affinity settings are found to be
    performant:

    .. code-block:: sh

        export OMP_PLACES=threads
        export OMP_PROC_BIND=spread

    However, the optimal settings may vary depending on the system.


For developers, :doc:`apidoc_cpp/apidoc_cpp` contains the full C++
API reference. To reuse C++ routines, ensure the linker finds ``libtrv``,
e.g. by providing the path to the static library (built following
:doc:`installation`) with the flag :code:`-L<path-to-libtrv.a>`
during compilation.


.. |Triumvirate| raw:: html

    <span style="font-variant: small-caps">Triumvirate</span>
