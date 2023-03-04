************
Installation
************

Python package
==============

.. image:: https://img.shields.io/pypi/v/Triumvirate?logo=PyPI&color=informational
    :target: https://pypi.org/project/Triumvirate
    :alt: PyPI
.. image:: https://img.shields.io/conda/vn/conda-forge/triumvirate?logo=conda-forge
    :target: https://anaconda.org/conda-forge/triumvirate
    :alt: conda

|br| |Triumvirate| as a Python package is distributed through `PyPI
<https://pypi.org/project/Triumvirate>`_ and `conda-forge
<https://anaconda.org/conda-forge/triumvirate>`_.

We recommend you install using `conda`:

.. code-block:: console

    $ conda install triumvirate

Alternatively, you may use `pip`, ideally in a virtual environment
(e.g. created with ``conda create -n <env>`` and activated with
``conda activate <env>``):

.. code-block:: console

    $ python -m pip install Triumvirate

.. attention::

    Please note the capitalisation of the package name in the commands above.

By default, the installed package from public distributions is built with
OpenMP enabled (if available).


C++ library & program
=====================

|Triumvirate| as either a C++ library or a 'black-box' C++ program can be
built using `make`.

You need to first obtain the source code by cloning into the GitHub repository
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
the end of the second line as shown above.

By default, the static library is compiled to ``build/lib/libtrv.a`` and the
binary executable is compiled to ``build/bin/triumvirate`` in the repository
directory.


Development mode
================

Both the Python package and the C++ library/program can be set up in
development mode with `make`. As in `C++ library & program`_ above, first
``git clone`` this repository and ``git checkout`` the branch/release you
would like to edit:

.. code-block:: console

    $ git clone git@github.com:MikeSWang/Triumvirate.git --branch <branch-or-release>
    $ cd Triumvirate

Then at the repository directory root run

.. code-block:: console

    $ make clean
    $ make [py|cpp]install [useomp=(true|1)]

where ``install`` builds both and ``pyinstall``/``cppinstall`` is for
Python/C++ build only; you may also replace this with ``cpplibinstall`` or
``cppappbuild`` as above to compile the C++ static library or binary executable
only. Again, to enable OpenMP parallelisation, append ``useomp=true`` or
``useomp=1`` to the end of the second line as shown above.

The latest release is on the ``main`` branch. The default ``Makefile``
(located at the repository diretory root) suits most use cases, but you may
modify it as appropriate for your need.


.. |Triumvirate| raw:: html

    <span style="font-variant: small-caps">Triumvirate</span>


.. |br| raw:: html

    <br/>
