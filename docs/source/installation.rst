************
Installation
************

Python package
==============

.. image:: https://img.shields.io/pypi/v/Triumvirate?logo=PyPI&color=informational
    :target: https://pypi.org/project/Triumvirate
    :alt: PyPI
.. image:: https://img.shields.io/conda/vn/msw/triumvirate?logo=Anaconda&color=informational
    :target: https://anaconda.org/msw/triumvirate
    :alt: conda

|br| |Triumvirate| as a Python package is distributed through `PyPI
<https://pypi.org/project/Triumvirate>`_ and `conda
<https://anaconda.org/msw/triumvirate>`_.

The recommended installation method is `conda` (from channel 'msw'):

.. code-block:: console

    $ conda install -c msw triumvirate

Alternatively, `pip` may be used, ideally in a virtual environment
(e.g. created with ``conda create -n <env>`` and activated with
``conda activate <env>``):

.. code-block:: console

    $ python -m pip install Triumvirate

.. attention::

    Note the capitalisation of the package name in the commands above.

By default, the installed package from public distributions is built with
OpenMP enabled (if available).


C++ library & program
=====================

|Triumvirate| as either a C++ library or a 'black-box' C++ program can be
built using `make`.

First, ensure that the required dependencies are installed (see
'`Dependencies`_' below).

Next, obtain the source code by cloning the GitHub repository
and change to the repository directory:

.. code-block:: console

    $ git clone git@github.com:MikeSWang/Triumvirate.git [--branch <branch-or-release>]
    $ cd Triumvirate

Then to compile, run:

.. code-block:: console

    $ make clean
    $ make cppinstall|cpplibinstall|cppappbuild [useomp=(true|1)]

Here ``cppinstall`` builds both a static library and a binary executable,
``cpplibinstall`` only the former and ``cppappbuild`` only the latter.
To enable OpenMP parallelisation, append ``useomp=true`` or ``useomp=1`` to
the end of the second line as shown above (see also '`OpenMP support`_' below).

By default, the static library is compiled to ``build/lib/libtrv.a`` and the
binary executable is compiled to ``build/bin/triumvirate`` in the repository
directory.

.. hint::

    The default |Makefile|_ (located at the repository diretory root)
    should work in most build environments, but you might need to modify it
    as needed for your own system.


Development mode
================

Both the Python package and the C++ library/program can be set up in
development mode with `make`. As in '`C++ library & program`_' above, first
check the required dependencies are installed (see '`Dependencies`_'
below), then ``git clone`` the GitHub repository and ``git checkout``
the branch/release to be edited:

.. code-block:: console

    $ git clone git@github.com:MikeSWang/Triumvirate.git --branch <branch-or-release>
    $ cd Triumvirate

Then at the repository directory root, run

.. code-block:: console

    $ make clean
    $ make [py|cpp]install [useomp=(true|1)]

where ``install`` builds both and ``pyinstall``/``cppinstall`` is for
Python/C++ build only; you may also replace this with ``cpplibinstall`` or
``cppappbuild`` above to compile the C++ static library or binary executable
only. Again, to enable OpenMP parallelisation, append ``useomp=true`` or
``useomp=1`` to the end of the second line as shown above (see also
'`OpenMP support`_' below).

.. note::

    The latest release is on the |main|_ branch. You might need to modify
    the default |Makefile|_ (located at the repository diretory root)
    as appropriate for your needs.


Dependencies
============

.. warning::

    If you encounter build issues associated with the ARM architecture
    on macOS, please consult the relevant sources for help.

There are two required dependencies when building the C++ program and the
Python package in development mode: the GSL and FFTW3 libraries.

If you are using a `conda` environment, both libraries can be installed with

.. code-block:: console

    $ conda install -c conda-forge gsl fftw

Outside `conda`, they can be installed using e.g. `apt` on Debian-based
Linux distributions such as Ubuntu or `brew` on macOS:

.. tabs::

    .. code-tab:: console Linux (e.g. Debian-based)

        $ [sudo] apt[-get] install [-y] libgsl(-dev|27) libfftw3-(dev|3)

    .. code-tab:: console macOS

        $ brew install gsl fftw

(select the appropriate package release in round brackets above).

To check the compilation flags and linker options needed when using these
libraries, use e.g. `pkg-config`:

.. code-block:: console

    $ pkg-config --cflags --libs gsl fftw3

.. hint::

    These options are automatically configured and used by the default
    |Makefile|_ (located at the repository diretory root) during
    build processes.


OpenMP support
==============

.. attention::

    Building the C++ program, or the Python package in development mode,
    requires a C++ compiler with OpenMP support. On Linux platforms,
    the GNU compiler with libgomp should suffice; by contrast, macOS
    systems may not come with OpenMP-supported compilers.

On Linux platforms, we recommend setting the following environmental variables
for building with OpenMP:

.. code-block:: console
    :caption: Linux

    # Optional: only if `CXX` previously set to a non-GNU compiler.
    $ export CXX=g++
    $ export CFLAGS="${CFLAGS} -fopenmp"
    # Optional: only if `-fopenmp` alone does not link properly.
    $ export LDFLAGS="${LDFLAGS} -lgomp"

On macOS systems, we recommend one first installs either the GNU compiler or
the LLVM compiler plus libomp using Homebrew,

.. code-block:: console
    :caption: macOS

    $ brew install gcc          # GNU compiler; or
    $ brew install llvm libomp  # LLVM compiler with libomp

and then set the relevant environmental variables,

.. tabs::

    .. code-tab:: console macOS with GNU compiler

        # HINT: assuming version 11 for brew formula 'gcc'
        $ export CXX=$(brew --prefix gcc)/bin/g++-11
        $ export CFLAGS="${CFLAGS} -fopenmp"

    .. code-tab:: console macOS with LLVM compiler plus libomp

        $ export CXX=$(brew --prefix llvm)/bin/clang++
        $ export CFLAGS="${CFLAGS} -Xpreprocessor -fopenmp"
        $ export LDFLAGS="${LDFLAGS} -L$(brew --prefix libomp)/lib -lomp"

These are the instructions for general OpenMP compilation; for the OpenMP
build of |Triumvirate|, the default |Makefile|_ (located at the repository
diretory root) assumes the GNU compiler and automatically sets the
appropriate comliation options; if you would like to use the LLVM compiler
on macOS, simply comment out the instructions for the GNU compiler and
uncomment the corresponding lines for the LLVM compiler.


Parallelised building
=====================

Building the C++ program, or the Python package in development mode,
can be parallelised by passing the ``-j[N] -O`` option to `make`,
where the optional parameter ``N`` is the number of concurrent jobs
(see also `'GNU Make Manual' <https://www.gnu.org/software/make/manual/
html_node/Options-Summary.html>`_).

Installing |Triumvirate| in development mode directly with `pip` in
editable mode can also be parallelised by setting the environmental
variable ``PY_BUILD_PARALLEL`` to ``-j[N]`` akin to the above:

.. code-block:: console

    $ export PY_BUILD_PARALLEL=-j   # use all available CPUs; or
    $ export PY_BUILD_PARALLEL=-j4  # use e.g. 4 CPUs


.. |Triumvirate| raw:: html

    <span style="font-variant: small-caps">Triumvirate</span>


.. |main| replace:: ``main``
.. _main: https://github.com/MikeSWang/Triumvirate/tree/main

.. |Makefile| replace:: ``Makefile``
.. _Makefile: _static/Makefile

.. |br| raw:: html

    <br/>
