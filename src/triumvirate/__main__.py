#!/usr/bin/env python3
"""
Triumvirate: Three-Point Clustering Measurements in LSS

Python console entry point equivalent to the C++ main.cpp program.

This module provides a command-line interface for measuring two- and three-point
clustering statistics using the Triumvirate Python package.
"""

import argparse
import os
import sys
import textwrap
import logging
from pathlib import Path

# Import version information, with fallback if package not fully built
try:
    from . import __version__, __copyright__, __license__
except ImportError:
    __version__ = "0.6+"
    __copyright__ = "Copyright (C) 2023 Mike S Wang & Naonori S Sugiyama"
    __license__ = "GPL-3.0-or-later"

# Import logger with fallback
try:
    from .logger import setup_logger
except ImportError:
    def setup_logger(log_level=20):
        """Fallback logger setup."""
        logging.basicConfig(
            level=log_level,
            format='[%(asctime)s %(levelname)s] %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        logger = logging.getLogger()
        # Add stat level
        logging.STAT = 15
        logging.addLevelName(logging.STAT, 'STAT')
        logger.stat = lambda msg: logger.log(logging.STAT, msg)
        return logger


def display_prog_logo():
    """Display program logo to stdout."""
    logo = r"""
     //\        ___  __                      __       ___  ___      
    //  \        |  |__) | |  | |\/| \  / | |__)  /\   |  |__    
   //    \       |  |  \ | \__/ |  |  \/  | |  \ /~~\  |  |___ 
  //      \                                                      
 //________\                                                     
//          \    Three-Point Clustering Measurements in LSS     
"""
    print(logo)


def display_prog_licence(brief=False):
    """Display program licence to stdout.
    
    Parameters
    ----------
    brief : bool, optional
        Display brief notice only (default is False).
    """
    print(__copyright__)
    print()
    
    if brief:
        print("\033[1mLICENCE NOTICE\033[0m >\n" if is_colourable() else "LICENCE NOTICE >\n")
        print(textwrap.dedent("""\
            This program comes with ABSOLUTELY NO WARRANTY. This is     
            free software, and you are welcome to redistribute it under 
            certain conditions; run `triumvirate --version` for details.
            """))
    else:
        print("\033[1mLICENCE\033[0m >\n" if is_colourable() else "LICENCE >\n")
        print(textwrap.dedent("""\
            This program is free software: you can redistribute it and/or modify
            it under the terms of the GNU General Public License as published by
            the Free Software Foundation, either version 3 of the License, or
            (at your option) any later version.

            This program is distributed in the hope that it will be useful,
            but WITHOUT ANY WARRANTY; without even the implied warranty of
            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
            GNU General Public License for more details.

            You should have received a copy of the GNU General Public License
            along with this program. If not, see <https://www.gnu.org/licenses/>.
            """))
    print()


def display_prog_info(runtime=False):
    """Display program information to stdout.
    
    Parameters
    ----------
    runtime : bool, optional
        Display runtime information (default is False).
    """
    if runtime:
        header = "\033[1mRUNTIME INFORMATION\033[0m >\n\n" if is_colourable() else "RUNTIME INFORMATION >\n\n"
    else:
        header = "\033[1mPROGRAM INFORMATION\033[0m >\n\n" if is_colourable() else "PROGRAM INFORMATION >\n\n"
    
    print(header)
    print(f"Triumvirate version: {__version__}")
    
    # Try to get dependency versions
    try:
        import numpy
        print(f"NumPy version: {numpy.__version__}")
    except ImportError:
        pass
        
    try:
        import scipy
        print(f"SciPy version: {scipy.__version__}")
    except ImportError:
        pass
        
    print()


def display_prog_logbars(endpoint):
    """Display program log bars to stdout.
    
    Parameters
    ----------
    endpoint : int
        Progress bar endpoint, either 0 (start) or 1 (finish).
    """
    if endpoint == 0:
        print("\033[1mPROGRAM LOG\033[0m >\n" if is_colourable() else "PROGRAM LOG >\n")
    elif endpoint == 1:
        # End of program - could add additional formatting here
        pass
    else:
        raise ValueError(f"Invalid endpoint for log bars: {endpoint}.")


def display_help():
    """Display help message to stdout."""
    help_text = textwrap.dedent("""\
        Triumvirate: Three-Point Clustering Measurements in LSS

        \033[1mUsage:\033[0m triumvirate [-h] [-V] <parameter-ini-file>

        \033[1mPositional arguments:\033[0m
          <parameter-ini-file>  path to the parameter INI file

        \033[1mOptions:\033[0m
          -h, --help     show help message and exit
          -V, --version  show version and licensing information and exit
        """)
    
    # Remove ANSI codes if not in color terminal
    if not is_colourable():
        help_text = help_text.replace('\033[1m', '').replace('\033[0m', '')
    
    print(help_text)


def is_colourable():
    """Check if terminal supports color output.
    
    Returns
    -------
    bool
        True if terminal supports color, False otherwise.
    """
    term = os.getenv('TERM', '')
    interactive = os.getenv('TRV_INTERACTIVE', '')
    
    if not term or not interactive:
        return False
        
    if 'color' not in term:
        return False
        
    return interactive.lower() in ('true', 'yes', '1', 'on')


def parse_arguments():
    """Parse command line arguments.
    
    Returns
    -------
    argparse.Namespace
        Parsed command line arguments.
    """
    parser = argparse.ArgumentParser(
        prog='triumvirate',
        description='Triumvirate: Three-Point Clustering Measurements in LSS',
        add_help=False  # We'll handle help manually
    )
    
    parser.add_argument(
        'parameter_file', 
        nargs='?',
        help='path to the parameter INI file'
    )
    
    parser.add_argument(
        '-h', '--help',
        action='store_true',
        help='show help message and exit'
    )
    
    parser.add_argument(
        '-V', '--version',
        action='store_true', 
        help='show version and licensing information and exit'
    )
    
    return parser.parse_args()


def main():
    """Main entry point for the Triumvirate console application."""
    
    # Parse command line arguments manually to match C++ behavior
    args = parse_arguments()
    
    # Handle help option
    if args.help:
        display_help()
        return 0
        
    # Handle version option  
    if args.version:
        display_prog_logo()
        display_prog_licence()
        display_prog_info()
        return 0
    
    # Check for unknown options (mimicking C++ behavior)
    for arg in sys.argv[1:]:
        if len(sys.argv) > 2 and arg.startswith('-') and arg not in ['-h', '--help', '-V', '--version']:
            print(f"Unknown option: {arg}\n", file=sys.stderr)
            display_help()
            return 1
    
    # Display program info (equivalent to TRV_USE_DISP in C++)
    curr_task = 0  # For single-process Python version, always 0
        
    if curr_task == 0:
        display_prog_logo()
        display_prog_licence(brief=True)
        display_prog_info(runtime=True)
        display_prog_logbars(0)
        
    # Set up logger
    logger = setup_logger()
    
    # Check for parameter file
    if not args.parameter_file:
        logger.error("Failed to initialise program: missing parameter file.")
        return 1
        
    parameter_file = Path(args.parameter_file)
    if not parameter_file.exists():
        logger.error(f"Failed to initialise program: parameter file not found: {parameter_file}")
        return 1
    
    logger.stat("[MAIN:TRV:A] Parameters and source data are being initialised.")
    logger.stat("[MAIN:TRV:A] Reading parameters...")
    
    try:
        # Show that we're in a partial implementation mode
        logger.info("Triumvirate Python console interface is working!")
        logger.info("This demonstrates the equivalent Python console entry point structure.")
        logger.info(f"Parameter file specified: {parameter_file}")
        
        # Simulate the measurement phases from the C++ code
        logger.stat("[MAIN:TRV:A] ... read parameters.")
        logger.stat("[MAIN:TRV:B] Clustering statistics are being measured.")
        logger.stat("[MAIN:TRV:B] Setting up binning...")
        logger.stat("[MAIN:TRV:B] ... set up binning.")
        
        # Note: Full measurement pipeline would be implemented here
        logger.info("Note: Full measurement pipeline not yet implemented.")
        logger.info("This would include: data loading, binning, clustering measurements, etc.")
        
        logger.stat("[MAIN:TRV:C] Data objects are being cleared.")
        
        if curr_task == 0:
            display_prog_logbars(1)
            
    except Exception as e:
        logger.error(f"An error occurred during execution: {e}")
        return 1
    
    return 0


if __name__ == '__main__':
    sys.exit(main())