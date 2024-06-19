"""Generate the window convolution coefficients for a convolved multipole
from a list of window function and unwindowed multipoles for three-point
clustering statistics.

Examples
--------

.. code-block:: console

    $ python gen_wc_coeffs.py 000 -u 000 110 220 -w 000 110 220

"""
from argparse import ArgumentParser
from itertools import product

import numpy as np

from triumvirate.winconv import (
    Multipole,
    WinConvTerm,
    calc_threept_winconv_coeff,
)


def configure():
    """Configure parameters from command-line arguments.

    Returns
    -------
    args : argparse.Namespace
        The command line arguments.

    """
    dochead = __doc__.split('\n\n')[0]

    parser = ArgumentParser(description=dochead)

    parser.add_argument(
        'convolved_multipole',
        help="convolved multipole"
    )

    parser.add_argument(
        '-w', '--wmultipoles',
        action='extend',
        nargs='+',
        help="window function multipoles"
    )

    parser.add_argument(
        '-u', '--umultipoles',
        action='extend',
        nargs='+',
        help="unwindowed multipoles"
    )

    parser.add_argument(
        '-r', '--sort',
        choices=['c', 'u', 'w'],
        default='c',
        help="sort by 'c' (coefficient), 'u' (unwindowed multipole), "
             "or 'w' (window function multipole)"
    )
    parser.add_argument(
        '-s', '--save',
        help="path to file for saving the generated terms"
    )

    args = parser.parse_args()

    return args


def main():
    """Main function.

    """
    args = configure()

    # List multipoles.
    convolved_multipole = Multipole(args.convolved_multipole)
    window_multipoles = [Multipole(wmpole) for wmpole in args.wmultipoles]
    unwindowed_multipoles = [Multipole(umpole) for umpole in args.umultipoles]

    # Tabulate terms.
    terms = []
    for (wmpole, umpole) in product(window_multipoles, unwindowed_multipoles):
        wc_coeff = calc_threept_winconv_coeff(
            convolved_multipole, wmpole, umpole, symb=True
        )
        if wc_coeff != 0.:
            terms.append(WinConvTerm(wmpole, umpole, wc_coeff))

    # Sort terms.
    if args.sort == 'c':
        terms = sorted(
            terms, key=lambda t: (-t.coeff, t.index_Z, t.index_Q,)
        )
    elif args.sort == 'u':
        terms = sorted(
            terms, key=lambda t: (t.index_Z, t.index_Q, -t.coeff,)
        )
    elif args.sort == 'w':
        terms = sorted(
            terms, key=lambda t: (t.index_Q, t.index_Z, -t.coeff,)
        )
    else:
        raise ValueError(f"Invalid sort option: '{args.sort}'.")

    # Save terms.
    if args.save:
        np.save(args.save, terms, allow_pickle=True)
        if args.save.endswith('.npy'):
            print(f"Saved terms to file: {args.save}")
        else:
            print(f"Saved terms to file: {args.save}.npy")

    # Print terms.
    maxlen_coeff = max(len(str(term.coeff)) for term in terms)
    fmt_str = f"{{:}}: {{: >{maxlen_coeff}}} Q_{{{{{{:}}}}}} Z_{{{{{{:}}}}}}"

    for term in terms:
        print(
            fmt_str.format(
                convolved_multipole.abstr,
                str(term.coeff), term.multipole_Q.abstr, term.multipole_Z.abstr
            )
        )


if __name__ == '__main__':
    main()
