![Triumvirate](https://github.com/MikeSWang/Triumvirate/raw/dev/docs/source/_static/Triumvirate.png)

[![arXiv eprint](
https://img.shields.io/badge/arXiv-yymm.nnnnn-important
)](https://arxiv.org/abs/yymm.nnnnn)
[![GitHub release (latest by date)](
https://img.shields.io/github/v/release/MikeSWang/Triumvirate?label=release
)](https://github.com/MikeSWang/Triumvirate/releases/latest)
[![Documentation status](
https://readthedocs.org/projects/triumvirate/badge/?version=latest
)](https://triumvirate.readthedocs.io/en/latest)
[![Build status](
https://travis-ci.com/MikeSWang/Triumvirate.svg?branch=dev
)](https://travis-ci.com/MikeSWang/Triumvirate)
[![Licence](
https://img.shields.io/badge/licence-GPLv3-informational
)](https://github.com/mikeswang/Triumvirate/tree/dev/LICENCE)


# Three-Point Clustering Correlators in Large-Scale Structure Analysis

<span style="font-variant: small-caps">Triumvirate</span> is a Python/C++
software package for the measurement and model comparison of
three-point clustering correlators in large-scale structure analysis.

## Pre-release guide

<mark style="text-color: orange">If you can see this text, please
note that the repository is currently in pre-release mode, with specific
instructions detailed below overriding any other instructions in the rest
of this guide.</mark>

After `git clone`-ing this repository, at the root of the directory
run the following commands to compile the code:
```
make clean
make all
```
If you only wish to use the Python (or C++) package, replace `make all`
with `make pyinstall` (or `make cppinstall`).  Depending on your
machine environment, certain parts of the [``Makefile``](
https://github.com/MikeSWang/Triumvirate/tree/dev/Makefile)
] may need replacing (e.g. the C++ compiler).

See Python scripts/notebooks in
[``examples/``](https://github.com/MikeSWang/Triumvirate/tree/dev/examples)
for simple examples of how to perform basic tasks.  For API documentation,
you may need to build it yourself using ``sphinx-autodoc``.

TODO: Add instructions for generating API docs here.

Please note that both testing and documentation are incomplete and there
may be unknown bugs.  Any feedback is welcome, and more importantly,
if you would like to contribute to testing and documentation please get
in touch.


## Installation


## Documentation

Comprehensive documentation can be found at
[mikeswang.github.io/Triumvirate](https://mikeswang.github.io/Triumvirate).
Tutorials and scripts can also be found in
[``examples/``](https://github.com/MikeSWang/Triumvirate/tree/dev/examples).


## Attribution

The underlying C++ code was originally developed by Naonori Sugiyama.
Please refer to/cite Sugiyama et al. (2018) [[doi:10.1093/mnras/sty3249](
https://doi.org/10.1093/mnras/sty3249)] accordingly.


## Licence

Copyright 2021–22. All rights reserved.
