============
Installation
============

Python package
==============

|Triumvirate| as a Python package is distributed through `PyPI
<https://pypi.org/project/Triumvirate>`_ and `conda-forge
<https://anaconda.org/conda-forge/triumvirate>`_.

By default, we recommend you install using `conda`:

.. code:: console

    $ conda install triumvirate

Alternatively, you may use `pip`, ideally in a virtual environment
(e.g. created with ``conda create -n <env>`` and activated with
``conda activate <env>``):

.. code:: console

    $ python -m pip install Triumvirate

.. attention::

    Please note the capitalisation of the package name in the commands above.


C++ program
===========

The C++ program can be customised and compiled to a binary executable using
`make`.

You need to first obtain the source code by cloning into the GitHub repository
and change to the branch corresponding to the version needed [1]_:

.. code:: console

    $ git clone git@github.com:MikeSWang/Triumvirate.git
    $ cd Triumvirate
    $ git checkout <version-branch>

Then to compile [2]_, run:

.. code:: console

    $ make clean
    $ make cppinstall

To enable OpenMP parallelisation, append ``useomp=true`` or ``useomp=1`` to
the end of the second line above.

By default, the binary executable is compiled to ``build/triumvirate`` in
the repository directory.

.. admonition:: C++ library

    In the future, the C++ code will also be released as a library.

.. [1] The latest release is on the ``main`` branch.

.. [2] The default ``Makefile`` (located in the repository root diretory)
       suits most use cases, but you may modify it as appropriate for your need.


Development mode
================

The Python package can be set up (with optional OpenMP parallelisation) in
development mode using `make` as above for C++, with the replacement of
``cppinstall`` with ``pyinstall``. The C++ program can be compiled
simultaneously if ``cppinstall``/``pyinstall`` is replaced with ``install``.


.. |Triumvirate| raw:: html

    <span style="font-variant: small-caps">Triumvirate</span>
