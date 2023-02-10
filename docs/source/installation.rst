************
Installation
************

Python package
==============

|Triumvirate| as a Python package is distributed through `PyPI
<https://pypi.org/project/Triumvirate>`_ and `conda-forge
<https://anaconda.org/conda-forge/triumvirate>`_.

We recommend you install using `conda`:

.. code:: console

    $ conda install triumvirate

Alternatively, you may use `pip`, ideally in a virtual environment
(e.g. created with ``conda create -n <env>`` and activated with
``conda activate <env>``):

.. code:: console

    $ python -m pip install Triumvirate

.. attention::

    Please note the capitalisation of the package name in the commands above.

By default, the installed package from public distributions is built with
OpenMP enabled.


C++ program
===========

|Triumvirate| as a 'black-box' C++ program is built as a binary executable
using `make`.

You need to first obtain the source code by cloning into the GitHub repository
and change to the repository directory:

.. code:: console

    $ git clone git@github.com:MikeSWang/Triumvirate.git --branch stable
    $ cd Triumvirate

Then to compile, run:

.. code:: console

    $ make clean
    $ make cppinstall [useomp=[true|1]]

To enable OpenMP parallelisation, append ``useomp=true`` or ``useomp=1`` to
the end of the second line as shown above.

By default, the binary executable is compiled to ``build/triumvirate`` in
the repository directory.

.. admonition:: C++ library

    In the future, the C++ code will also be released as a library.


Development mode
================

Both the Python package and the C++ program can be set up in development
mode with `make`. As in `C++ program`_ above, first ``git clone`` this
repository and ``git checkout`` the branch/release you would like to edit

.. code:: console

    $ git clone git@github.com:MikeSWang/Triumvirate.git --branch <branch-or-release>
    $ cd Triumvirate

then at the repository directory root run

.. code:: console

    $ make clean
    $ make [py|cpp]install [useomp=[true|1]]

where ``install`` builds both and ``pyinstall``/``cppinstall`` is for
Python/C++ build only. Again, to enable OpenMP parallelisation, append
``useomp=true`` or ``useomp=1`` to the end of the second line as shown above.

The latest release is on the ``main`` branch. The default ``Makefile``
(located at the repository diretory root) suits most use cases, but you may
modify it as appropriate for your need.


.. |Triumvirate| raw:: html

    <span style="font-variant: small-caps">Triumvirate</span>
