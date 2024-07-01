************
Installation
************

Dependencies
============

|Triumvirate| as a compiled C++ library/program or as a Python package built
*from a source distribution* has the following dependencies.


GSL & FFTW3
-----------

The GSL and FFTW3 libraries are required dependencies. If they have already
been installed on your system, please ignore this subsection.

In a Conda environment, both libraries can be installed with

.. code-block:: console

    $ conda install -c conda-forge gsl fftw

Outside Conda environments, they can be installed using a package manager:

.. tab-set::

    .. tab-item:: Linux

        .. code-block:: console

            # Debian-based OS (e.g. Ubuntu); or
            $ [sudo] apt[-get] install [-y] libgsl(-dev|27) libfftw3-(dev|3)
            # Red Hat--based OS (e.g. CentOS)
            $ [sudo] yum install [-y] gsl[-devel] fftw3-devel

    .. tab-item:: macOS

        .. code-block:: console

            $ brew install gsl fftw

(select the appropriate option/version in the brackets above).

Another option is to build these libraries from the source, following the
instructions from `'GSL - GNU Scientific Library'
<https://www.gnu.org/software/gsl/>`_ and `'FFTW' <https://www.fftw.org>`_.
In this case, path environmental variables may need to be set/modified
to ensure that they are included and linked during |Triumvirate| compilation.

To check the compilation options for using these libraries, `pkg-config`
may be of help:

.. code-block:: console

    $ pkg-config --cflags --libs gsl fftw3

However, it is not guaranteed that `pkg-config` always detects the
relevant settings, and it is important to ensure the set-up of `pkg-config`
matches that of the dependencies (e.g. both are installed by Conda in the
same Conda environment).

.. hint::

    These library compilation options are automatically configured and set
    by the default |Makefile|_ in the package distribution(s)
    (located at the repository directory root) during build processes.


OpenMP library
--------------

An OpenMP library is an optional dependency for OpenMP-enabled installation.
If a compatible OpenMP library has been distributed with your compiler,
please ignore this subsection.

For GCC compilers, the relevant library is ``libgomp`` (typically distributed
with the compiler). On Linux systems, it can be additionally installed either
in a Conda environment with

.. code-block:: console
    :caption: Linux

    $ conda install -c conda-forge libgomp

or outside Conda environments using a package manager:

.. code-block:: console
    :caption: Linux

    # Debian-based OS (e.g. Ubuntu); or
    $ [sudo] apt[-get] install [-y] libgomp1
    # Red Hat--based OS (e.g. CentOS)
    $ [sudo] yum install [-y] libgomp

For LLVM compilers, the relevant library is ``libomp``. It can be installed
either in a Conda environment with

.. code-block:: console

    $ conda install -c conda-forge llvm-openmp

or outside Conda environments using the Homebrew package manager:

.. code-block:: console

    $ brew install libomp


Python package
==============

.. image:: https://img.shields.io/pypi/v/Triumvirate?logo=PyPI&color=informational
    :target: https://pypi.org/project/Triumvirate
    :alt: PyPI

.. image:: https://img.shields.io/conda/v/msw/triumvirate?logo=Anaconda&color=informational
    :target: https://anaconda.org/msw/triumvirate
    :alt: Conda

|br| |Triumvirate| as a Python package is distributed through |PyPI-repo|
and |conda-repo|. For dependency management, it is recommended that a
virtual environment should be created for installing and using the package
(e.g. a Conda environment created with ``conda create -n <env>`` and
activated with ``conda activate <env>``).

To install from |PyPI-repo|, execute in shell:

.. code-block:: console

    $ python -m pip install triumvirate

To install using |conda-repo|, execute in shell:

.. code-block:: console

    $ conda install -c msw triumvirate

By default, the package is installed with OpenMP enabled if it is supported.

.. hint::

    Conda packages are built with dependencies such as ``numpy`` and
    ``scipy`` sourced from the ``conda-forge`` channel. For consistency
    and avoidance of dependency conflicts, it is recommended that
    ``conda-forge`` should be set as the highest-priority channel,

    .. code-block:: console

        $ conda config --prepend channels conda-forge

    and optionally, it is good practice to use ``strict`` channel priority,

    .. code-block:: console

        $ conda config --set channel_priority strict


C++ library & program
=====================

|Triumvirate| as either a static library or a binary executable can be
built using `make`, provided that dependency requirements are satisfied
(see '`Dependencies`_' above).

First, obtain the source by cloning the GitHub repository and change into
its local directory path:

.. code-block:: console

    $ git clone git@github.com:MikeSWang/Triumvirate.git [--branch <branch-or-release>]
    $ cd Triumvirate

Then, execute in shell:

.. code-block:: console

    $ make clean
    $ make cppinstall|cpplibinstall|cppappbuild [useomp=(true|1)]

Here ``cppinstall`` builds both the static library and the binary executable,
``cpplibinstall`` only the former and ``cppappbuild`` only the latter.
To enable OpenMP parallelisation (see '`OpenMP support`_' below), append
``useomp=true`` or ``useomp=1`` to the end of the second line as shown above.

By default, the static library is compiled to ``build/lib/libtrv.a`` and
the binary executable is compiled to ``build/bin/triumvirate`` in the
repository directory.

.. hint::

    The default |Makefile|_ (located at the repository directory root)
    should work in most build environments, but may need to be modified
    as appropriate for the build environment.


Development mode
================

.. image:: https://img.shields.io/github/v/release/MikeSWang/Triumvirate?display_name=tag&sort=semver&logo=Git
    :target: https://github.com/MikeSWang/Triumvirate/releases/latest
    :alt: Release

|br| Both the Python package and the C++ library/program can be set up in
development mode with `make`, provided that dependency requirements are
satisfied (see '`Dependencies`_' above).

As in '`C++ library & program`_' above, first ``git clone`` the desired
branch/release from the GitHub repository and change into the repository
directory path:

.. code-block:: console

    $ git clone git@github.com:MikeSWang/Triumvirate.git --branch <branch-or-release>
    $ cd Triumvirate

Then, execute in shell:

.. code-block:: console

    $ make clean
    $ make ([py|cpp]install)|(cpp[libinstall|appbuild]) [useomp=(true|1)]

where ``cpplibinstall`` or ``cppappbuild`` respectively builds the C++
static library or binary executable only, ``cppinstall`` builds both,
``pyinstall`` builds the Python package only, and ``install`` builds
all of the above. As before, to enable OpenMP parallelisation (see
'`OpenMP support`_' below), append ``useomp=true`` or ``useomp=1`` to the
end of the second line as shown above.

.. note::

    The latest release is on the |main|_ branch. The default |Makefile|_
    (located at the repository directory root) should work in most build
    environments, but may need to be modified as appropriate.


OpenMP support
==============

When building from a source distribution with OpenMP parallelisation,
the compiler must support OpenMP and the OpenMP library
(see '`OpenMP library`_') must be compatible.

By default, |Makefile|_ in source distributions and `setup.py` in Python
package distributions both assume the GCC compiler and OpenMP library
and configure the OpenMP-enabled compilation options accordingly.


Using `make`
------------

By default, OpenMP is *disabled* for `make`-based installation. To *enable*
OpenMP parallelisation, pass ``useomp=true`` or ``useomp=1`` to `make`.

To override the compilation settings used in the default |Makefile|_, set the
environmental variables as shown in the following examples for macOS:

.. tab-set::

    .. tab-item:: GCC
        :sync: gcc

        .. code-block:: console

            # Set GCC compiler (version 12 assumed here).
            $ export CXX=$(brew --prefix gcc)/bin/g++-12

    .. tab-item:: LLVM
        :sync: llvm

        .. code-block:: console

            # Assume Homebrew LLVM compiler and OpenMP library instead here.
            $ export CXX=$(brew --prefix llvm)/bin/clang++
            # Set OpenMP compilation flags.
            $ export CXXFLAGS_OMP="-I$(brew --prefix libomp)/include -fopenmp"
            # Set OpenMP linker flags.
            $ export LDFLAGS_OMP="-L$(brew --prefix libomp)/lib -lomp"

These commands are also included in the default |Makefile|_ (though commented
out as they are an alternative to the default GCC set-up).


Python setup
------------

By default, OpenMP support is automatically detected for Python installation
(except when building through `make`). To *disable* OpenMP parallelisation
explicitly, set the environmental variable ``PY_NO_OMP`` with
:code:`export PY_NO_OMP=''` (and unset with :code:`unset PY_NO_OMP` to
re-enable it); to *enforce* OpenMP parallelisation explicitly, set the
environmental variable ``PY_OMP`` (to any non-empty value).

To override the compilation settings used by ``setup.py`` (e.g. to use a
different compiler suite to the default GCC), set the environmental variables
as shown in the following examples (for macOS with Homebrew package manager):

.. tab-set::

    .. tab-item:: GCC
        :sync: gcc

        .. code-block:: console

            # Set GCC compiler (version 12 assumed here).
            $ export PY_CXX=$(brew --prefix gcc)/bin/g++-12

    .. tab-item:: LLVM
        :sync: llvm

        .. code-block:: console

            # Set LLVM compiler.
            $ export PY_CXX=$(brew --prefix llvm)/bin/clang++
            # Set OpenMP compilation flags.
            $ export PY_CXXFLAGS_OMP="-I$(brew --prefix libomp)/include -fopenmp"
            # Set OpenMP linker flags.
            $ export PY_LDFLAGS_OMP="-L$(brew --prefix libomp)/lib -lomp"


Parallelised building
=====================

Building the C++ library/program or the Python package from a
source distribution can be parallelised.

When using `make`, pass the ``-j[N]`` option where the optional parameter
``N`` is the number of concurrent jobs (see also `'GNU Make Manual'
<https://www.gnu.org/software/make/manual/html_node/Options-Summary.html>`_).

For the Python setup, set the environmental variable ``PY_BUILD_PARALLEL``
to ``-j[N]`` akin to the above, e.g. :code:`export PY_BUILD_PARALLEL=-j`
to use all available CPUs or :code:`export PY_BUILD_PARALLEL=-j4` to use four.


.. |Triumvirate| raw:: html

    <span style="font-variant: small-caps">Triumvirate</span>


.. |PyPI-repo| replace:: PyPI
.. _PyPI-repo: https://pypi.org/project/Triumvirate

.. |conda-repo| replace:: Conda
.. _conda-repo: https://anaconda.org/msw/triumvirate

.. |main| replace:: ``main``
.. _main: https://github.com/MikeSWang/Triumvirate/tree/main

.. |Makefile| replace:: ``Makefile``
.. _Makefile: _static/Makefile

.. |br| raw:: html

    <br/>
