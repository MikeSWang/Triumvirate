"""
Three-Point Clustering Statistics (:mod:`~triumvirate.threept`)
==========================================================================

Measure three-point clustering statistics from catalogues.

Local plane-parallel estimators
-------------------------------

.. autosummary::
    compute_bispec
    compute_3pcf

Global plane-parallel estimators
--------------------------------

.. autosummary::
    compute_bispec_in_gpp_box
    compute_3pcf_in_gpp_box

Window function estimator
-------------------------

.. autosummary::
    compute_3pcf_window

"""
import warnings
from pathlib import Path

import numpy as np

from triumvirate._threept import (
    _calc_bispec_normalisation_from_mesh,
    _calc_bispec_normalisation_from_particles,
    _compute_bispec,
    # _compute_bispec_for_los_choice,
    _compute_bispec_in_gpp_box,
    _compute_3pcf,
    _compute_3pcf_in_gpp_box,
    _compute_3pcf_window,
)
from triumvirate.dataobjs import Binning
from triumvirate.parameters import (
    _modify_measurement_parameters,
    _modify_sampling_parameters,
    fetch_paramset_template,
    ParameterSet,
)


def _amalgamate_parameters(paramset=None, params_sampling=None,
                           degrees=None, wa_orders=None, binning=None,
                           form=None, idx_bin=None, **type_kwargs):
    """Amalgamate a parameter set with overriding sampling parameters
    and measurement parameters.

    Parameters
    ----------
    paramset : :class:`~triumvirate.parameters.ParameterSet`, optional
        Full parameter set.  If `None` (default), some other parameters
        should be provided (see below).
    params_sampling : dict, optional
        Dictionary containing a subset of the following entries
        for sampling parameters---

        - 'alignment': {'centre', 'pad'};
        - 'boxsize': sequence of [float, float, float];
        - 'ngrid': sequence of [int, int, int];
        - 'assignment': {'ngp', 'cic', 'tsc', 'pcs'};
        - 'interlace': bool;

        and one and only one of the following when 'alignment' is 'pad'---

        - 'boxpad': float;
        - 'gridpad': float.

        This will override corresponding entries in `paramset`.
    degrees : tuple of (int, int, int) or str of length 3, optional
        Multipole degrees either as a tuple ('ell1', 'ell2', 'ELL') or
        as a string of length 3.  If not `None` (default), this will
        override ``paramset['degrees']`` entries.  If a string, multipole
        degrees are assumed to be single-digit integers.
    wa_orders : tuple of (int, int) or str of length 2, optional
        Wide-angle correction orders either as a tuple ('i_wa', 'j_wa') or
        as a string of length 2.  If not `None` (default), this will
        override ``paramset['wa_orders']`` entries.  If a string,
        multipole degrees are assumed to be single-digit integers.
    binning : :class:`~.triumvirate.dataobjs.Binning`, optional
        Binning (default is `None`).
    form : {'diag', 'full'}, optional
        Binning form of the measurements.  If not `None` (default),
        this will override ``paramset['form']``.
    idx_bin : int, optional
        Fixed bin index for the first coordinate dimension when binning
        `form` is 'full'.  If not `None` (default), this will
        override ``paramset['idx_bin']``.
    **type_kwargs
        `catalogue_type` and `statistic_type` parameters to be filled in.

    Returns
    -------
    paramset : :class:`~triumvirate.parameters.ParameterSet`, optional
        Amalgamated full parameter set.

    Raises
    ------
    ValueError
        When `paramset` is `None` while `params_sampling`, `degrees`,
        `binning` or `form` is also `None`, or while `form` is 'full' and
        `idx_bin` is `None`.

    """
    if paramset is None and params_sampling is None:
        raise ValueError(
            "Either `paramset` or `sampling_params` must be provided."
        )
    if paramset is None and None in [degrees, binning, form]:
        raise ValueError(
            "`degrees`, `binning` and `form` must be provided "
            "when `paramset` is None."
        )
    if paramset is None and form == 'full' and idx_bin is None:
        raise ValueError(
            "`idx_bin` must be provided when `paramset` is None "
            "and `form` is 'full'."
        )

    if paramset is None:
        paramset, defaults = fetch_paramset_template('dict', ret_defaults=True)
    else:
        # paramset = paramset
        defaults = None

    if params_sampling is not None:
        paramset, defaults = _modify_sampling_parameters(
            paramset, params_sampling,
            params_default=defaults, ret_defaults=True
        )

    params_measure = {}
    if degrees is not None:
        ell1, ell2, ELL = degrees
        params_measure.update({
            'ell1': int(ell1), 'ell2': int(ell2), 'ELL': int(ELL)
        })
    if wa_orders is not None:
        i_wa, j_wa = wa_orders
        params_measure.update({
            'i_wa': int(i_wa), 'j_wa': int(j_wa)
        })
    if binning is not None:
        params_measure.update({
            'bin_min': binning.bin_min,
            'bin_max': binning.bin_max,
            'num_bins': binning.num_bins,
        })
    if form is not None:
        params_measure.update({'form': form})
    if idx_bin is not None:
        params_measure.update({'idx_bin': idx_bin})

    paramset, defaults = _modify_measurement_parameters(
        paramset, params_measure,
        params_default=defaults, ret_defaults=True
    )

    # Implicitly validate parameters.
    paramset.update(type_kwargs)

    if defaults:
        warnings.warn(
            "The following parameter default values "
            f"are unchanged: {defaults}. "
            "Not all parameters are necessarily used."
        )

    return paramset


def _get_measurement_filename(paramset):
    """Get output measurement filename.

    Parameters
    ----------
    paramset : :class:`~triumvirate.parameters.ParameterSet`
        Parameter set.

    Returns
    -------
    str
        Output measurement filename.

    Raises
    ------
    ValueError
        When `paramset` indicates the measurements are not
        three-point statistics.

    """
    multipole = "{ell1}{ell2}{ELL}".format(**paramset['degrees'])

    if paramset['form'] == 'diag':
        binform = "_diag"
    if paramset['form'] == 'full':
        binform = "_bin{:d}".format(paramset['idx_bin'])

    output_tag = paramset['tags']['output'] or ""

    if paramset['statistic_type'] == 'bispec':
        return "bk{}{}{}".format(multipole, binform, output_tag)
    if paramset['statistic_type'] == '3pcf':
        return "zeta{}{}{}".format(multipole, binform, output_tag)
    if paramset['statistic_type'] == '3pcf-win':
        return "zetaw{}{}{}".format(multipole, binform, output_tag)
    if paramset['statistic_type'] == '3pcf-win-wa':
        wa_tag = "_wa{i}{j}".format(**paramset['wa_orders'])
        return "zetaw{}{}{}{}".format(multipole, binform, wa_tag, output_tag)

    raise ValueError(
        "`paramset` 'statistic_type' does not correspond to a "
        "recognised three-point statistic."
    )


def _print_measurement_header(paramset, norm_factor, norm_factor_alt):
    """Print three-point statistic measurement header including
    sampling parameters, normalisation factors and data table columns.

    Parameters
    ----------
    paramset : :class:`~triumvirate.parameters.ParameterSet`
        Parameter set.
    norm_factor, norm_factor_alt : float
        Normalisation factor used and the alternative.

    Returns
    -------
    text_header : str
        Measurement information as a header string.

    Raises
    ------
    ValueError
        When `paramset` indicates the measurements are not
        three-point statistics.

    """
    if paramset['norm_convention'] == 'particle':
        norm_convention = 'particle'
        norm_convention_alt = 'mesh'
    if paramset['norm_convention'] == 'mesh':
        norm_convention = 'mesh'
        norm_convention_alt = 'particle'

    if paramset['npoint'] != '3pt':
        raise ValueError(
            "`paramset` 'statistic_type' does not correspond to a "
            "recognised three-point statistic."
        )

    multipole = "{ell1}{ell2}{ELL}".format(**paramset['degrees'])
    if paramset['space'] == 'fourier':
        datatab_colnames = [
            "k1_cen", "k1_eff", "k2_cen", "k2_eff", "nmodes",
            "Re{{bk{}_raw}}".format(multipole),
            "Im{{bk{}_raw}}".format(multipole),
            "Re{{bk{}_shot}}".format(multipole),
            "Im{{bk{}_shot}}".format(multipole)
        ]
    if paramset['space'] == 'config':
        datatab_colnames = [
            "r1_cen", "r1_eff", "r2_cen", "r2_eff", "npairs",
            "Re{{zeta{}_raw}}".format(multipole),
            "Im{{zeta{}_raw}}".format(multipole),
            "Re{{zeta{}_shot}}".format(multipole),
            "Im{{zeta{}_shot}}".format(multipole)
        ]

    text_lines = [
        "Box size: ({:.3f}, {:.3f}, {:.3f})".format(
            *[paramset['boxsize'][ax] for ax in ['x', 'y', 'z']]
        ),
        "Box alignment: {}".format(paramset['alignment']),
        "Mesh number: ({:d}, {:d}, {:d})".format(
            *[paramset['ngrid'][ax] for ax in ['x', 'y', 'z']]
        ),
        "Mesh assignment and interlacing: {}, {}".format(
            paramset['assignment'], paramset['interlace']
        ),
        "Normalisation factor: "
        "{:.9e} ({}-based, used), {:.9e} ({}-based, alternative)".format(
            norm_factor, norm_convention,
            norm_factor_alt, norm_convention_alt
        ),
        ", ".join([
            "[{:d}] {}".format(colidx, colname)
            for colidx, colname in enumerate(datatab_colnames)
        ])
    ]

    text_header = "\n".join(text_lines)

    return text_header


def _assemble_measurement_datatab(measurements, paramset):
    """Assemble measurement data table.

    Parameters
    ----------
    measurements : dict
        Measurement results.
    paramset : :class:`~triumvirate.parameters.ParameterSet`
        Parameter set.

    Returns
    -------
    datatab : :class:`numpy.ndarray`
        Column-major measurement data table.

    Raises
    ------
    ValueError
        When `paramset` indicates the measurements are not
        three-point statistics.

    """
    if paramset['npoint'] != '3pt':
        raise ValueError(
            "`paramset` 'statistic_type' does not correspond to a "
            "recognised three-point statistic."
        )
    if paramset['space'] == 'fourier':
        datatab = np.transpose([
            measurements['k1_bin'], measurements['k1_eff'],
            measurements['k2_bin'], measurements['k2_eff'],
            measurements['nmodes'],
            measurements['bk_raw'].real, measurements['bk_raw'].imag,
            measurements['bk_shot'].real, measurements['bk_shot'].imag,
        ])
    if paramset['space'] == 'config':
        datatab = np.transpose([
            measurements['r1_bin'], measurements['r1_eff'],
            measurements['r2_bin'], measurements['r2_eff'],
            measurements['npairs'],
            measurements['zeta_raw'].real, measurements['zeta_raw'].imag,
            measurements['zeta_shot'].real, measurements['zeta_shot'].imag,
        ])

    return datatab


# ========================================================================
# Survey statistics
# ========================================================================

def _compute_3pt_stats_survey_like(threept_algofunc,
                                   catalogue_data, catalogue_rand,
                                   los_data=None, los_rand=None,
                                   paramset=None, params_sampling=None,
                                   degrees=None, binning=None,
                                   form=None, idx_bin=None, types=None,
                                   save=False, logger=None):
    """Compute three-point statistics from survey-like data and random
    catalogues in the local plane-parallel approximation.

    Parameters
    ----------
    threept_algofunc : callable
        Three-point statistic algorithmic function.
    catalogue_data : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Data-source catalogue.
    catalogue_rand : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Random-source catalogue.
    los_data : (N, 3) array of float, optional
        Specified lines of sight for the data-source catalogue.
        If `None` (default), this is automatically computed using
        :meth:`~triumvirate.catalogue.ParticleCatalogue.compute_los`.
    los_rand : (N, 3) array of float, optional
        Specified lines of sight for the random-source catalogue.
        If `None` (default), this is automatically computed using
        :meth:`~triumvirate.catalogue.ParticleCatalogue.compute_los`.
    paramset : :class:`~triumvirate.parameters.ParameterSet`, optional
        Full parameter set.  If `None` (default), `degrees`, `binning`,
        `form` and `params_sampling` should be provided; `idx_bin`
        should be provided when `form` is 'full'.
    params_sampling : dict, optional
        Dictionary containing a subset of the following entries
        for sampling parameters---

        - 'alignment': {'centre', 'pad'};
        - 'boxsize': sequence of [float, float, float];
        - 'ngrid': sequence of [int, int, int];
        - 'assignment': {'ngp', 'cic', 'tsc', 'pcs'};
        - 'interlace': bool;

        and one and only one of the following when 'alignment' is 'pad'---

        - 'boxpad': float;
        - 'gridpad': float.

        This will override corresponding entries in `paramset`.
    degrees : tuple of (int, int, int) or str of length 3, optional
        Multipole degrees either as a tuple ('ell1', 'ell2', 'ELL') or
        as a string of length 3.  If not `None` (default), this will
        override ``paramset['degrees']`` entries.  If a string, multipole
        degrees are assumed to be single-digit integers.
    binning : :class:`~triumvirate.dataobjs.Binning`, optional
        Binning for the measurements.  If `None` (default),
        this is constructed from `paramset`.
    form : {'diag', 'full'}, optional
        Binning form of the measurements.  If not `None` (default),
        this will override ``paramset['form']``.
    idx_bin : int, optional
        Fixed bin index for the first coordinate dimension
        when binning `form` is 'full'.  If not `None` (default),
        this will override ``paramset['idx_bin']``.
    types : dict, optional
        `catalogue_type` and `statistic_type` (default is `None`).
        This should be set by the caller of this function.
    save : {'.txt', '.npz', False}, optional
        If not `False` (default), save the measurements as a '.txt' file
        or in '.npz' format.  The save path is determined from `paramset`
        (if unset, a default file in the current working directory is
        used).
    logger : :class:`logging.Logger`, optional
        Logger (default is `None`).

    Returns
    -------
    results : dict of {str: :class:`numpy.ndarray`}
        Measurement results.

    """
    # --------------------------------------------------------------------
    # Initialisation
    # --------------------------------------------------------------------

    # -- Parameters ------------------------------------------------------

    paramset = _amalgamate_parameters(
        paramset=paramset, params_sampling=params_sampling,
        degrees=degrees, binning=binning, form=form, idx_bin=idx_bin, **types
    )

    if isinstance(paramset, dict):  # likely redundant but safe
        paramset = ParameterSet(param_dict=paramset, logger=logger)

    if logger:
        logger.info("Parameter set have been initialised.")

    # -- Data ------------------------------------------------------------

    # Set up binning.
    if binning is None:
        binning = Binning.from_parameter_set(paramset)

    if logger:
        logger.info("Binning has been initialised.")

    # Set up lines of sight.
    if los_data is None:
        los_data = catalogue_data.compute_los()
    if los_rand is None:
        los_rand = catalogue_rand.compute_los()

    los_data = np.ascontiguousarray(los_data)
    los_rand = np.ascontiguousarray(los_rand)

    if logger:
        logger.info("Lines of sight have been initialised.")

    # Set up box alignment.
    if paramset['alignment'] == 'centre':
        catalogue_data.centre(
            [paramset['boxsize'][ax] for ax in ['x', 'y', 'z']],
            catalogue_ref=catalogue_rand
        )
    if paramset['alignment'] == 'pad':
        if paramset['padscale'] == 'box':
            kwargs = {'boxsize_pad': paramset['padfactor']}
        if paramset['padscale'] == 'grid':
            kwargs = {
                'ngrid': [paramset['ngrid'][ax] for ax in ['x', 'y', 'z']],
                'ngrid_pad': paramset['padfactor']
            }
        catalogue_data.pad(
            [paramset['boxsize'][ax] for ax in ['x', 'y', 'z']],
            catalogue_ref=catalogue_rand, **kwargs
        )

    if logger:
        logger.info("Catalogues have been aligned.")

    # --------------------------------------------------------------------
    # Measurements
    # --------------------------------------------------------------------

    # Prepare catalogues.
    if logger:
        logger.info(
            "Preparing catalogue for clustering algorithm...",
            cpp_state='start'
        )

    particles_data = \
        catalogue_data._convert_to_cpp_catalogue(verbose=paramset['verbose'])
    particles_rand = \
        catalogue_rand._convert_to_cpp_catalogue(verbose=paramset['verbose'])

    if logger:
        logger.info(
            "... prepared catalogue for clustering algorithm.",
            cpp_state='end'
        )

    # Set up constants.
    alpha = catalogue_data.wtotal / catalogue_rand.wtotal

    if logger:
        logger.info("Alpha contrast: %.6e.", alpha)

    if paramset['norm_convention'] == 'particle':
        norm_factor = _calc_bispec_normalisation_from_particles(
            particles_rand, alpha
        )
        norm_factor_alt = _calc_bispec_normalisation_from_mesh(
            particles_rand, paramset, alpha
        )
    if paramset['norm_convention'] == 'mesh':
        norm_factor = _calc_bispec_normalisation_from_mesh(
            particles_rand, paramset, alpha
        )
        norm_factor_alt = _calc_bispec_normalisation_from_particles(
            particles_rand, alpha
        )

    if logger:
        logger.info(
            "Normalisation factors: %.6e (used), %.6e (alternative).",
            norm_factor, norm_factor_alt
        )

    # Perform measurement.
    if logger:
        logger.info("Measuring clustering statistics...", cpp_state='start')

    results = threept_algofunc(
        particles_data, particles_rand, los_data, los_rand,
        paramset, binning, norm_factor
    )

    if logger:
        logger.info("... measured clustering statistics.", cpp_state='end')

    if save:
        odirpath = paramset['directories']['measurements'] or ""
        header = "\n".join([
            catalogue_data.write_attrs_as_header(catalogue_ref=catalogue_rand),
            _print_measurement_header(paramset, norm_factor, norm_factor_alt),
        ])
        if save.lower() == '.txt':
            datatab = _assemble_measurement_datatab(results, paramset)
            datafmt = '\t'.join(
                ['%.9e'] * 4 + ['%10d'] + ['% .9e'] * (datatab.shape[-1] - 5)
            )
            ofilename = _get_measurement_filename(paramset)
            ofilepath = Path(odirpath, ofilename).with_suffix('.txt')
            np.savetxt(
                ofilepath, datatab, fmt=datafmt, header=header, delimiter='\t'
            )
        elif save.lower().endswith('.npz'):
            results.update({'header': header})
            ofilename = _get_measurement_filename(paramset)
            ofilepath = Path(odirpath, ofilename).with_suffix('.npz')
            np.savez(ofilepath, **results)
        else:
            raise ValueError(
                f"Unrecognised save format for measurements: {save}."
            )

        if logger:
            logger.info("Measurements saved to %s.", ofilepath)

    return results


def compute_bispec(catalogue_data, catalogue_rand,
                   los_data=None, los_rand=None,
                   degrees=None, binning=None, form=None, idx_bin=None,
                   sampling_params=None,
                   paramset=None,
                   save=False, logger=None):
    """Compute bispectrum from survey-like data and random catalogues
    in the local plane-parallel approximation.

    Parameters
    ----------
    catalogue_data : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Data-source catalogue.
    catalogue_rand : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Random-source catalogue.
    los_data : (N, 3) array of float, optional
        Specified lines of sight for the data-source catalogue.
        If `None` (default), this is automatically computed using
        :meth:`~triumvirate.catalogue.ParticleCatalogue.compute_los`.
    los_rand : (N, 3) array of float, optional
        Specified lines of sight for the random-source catalogue.
        If `None` (default), this is automatically computed using
        :meth:`~triumvirate.catalogue.ParticleCatalogue.compute_los`.
    degrees : tuple of (int, int, int) or str of length 3, optional
        Multipole degrees either as a tuple ('ell1', 'ell2', 'ELL') or
        as a string of length 3.  If not `None` (default), this will
        override ``paramset['degrees']`` entries.  If a string, multipole
        degrees are assumed to be single-digit integers.
    binning : :class:`~triumvirate.dataobjs.Binning`, optional
        Binning for the measurements.  If `None` (default),
        this is constructed from `paramset`.
    form : {'diag', 'full'}, optional
        Binning form of the measurements.  If not `None` (default),
        this will override ``paramset['form']``.
    idx_bin : int, optional
        Fixed bin index for the first coordinate dimension
        when binning `form` is 'full'.  If not `None` (default),
        this will override ``paramset['idx_bin']``.
    sampling_params : dict, optional
        Dictionary containing a subset of the following entries
        for sampling parameters---

        - 'alignment': {'centre', 'pad'};
        - 'boxsize': sequence of [float, float, float];
        - 'ngrid': sequence of [int, int, int];
        - 'assignment': {'ngp', 'cic', 'tsc', 'pcs'};
        - 'interlace': bool;

        and one and only one of the following when 'alignment' is 'pad'---

        - 'boxpad': float;
        - 'gridpad': float.

        This will override corresponding entries in `paramset`.
    paramset : :class:`~triumvirate.parameters.ParameterSet`, optional
        Full parameter set (default is `None`).  This is used in lieu of
        `degrees`, `binning`, `form`, `idx_bin` or `sampling_params`.
    save : {'.txt', '.npz', False}, optional
        If not `False` (default), save the measurements as a '.txt' file
        or in '.npz' format.  The save path is determined from `paramset`
        (if unset, a default file in the current working directory is
        used).
    logger : :class:`logging.Logger`, optional
        Logger (default is `None`).

    Returns
    -------
    results : dict of {str: :class:`numpy.ndarray`}
        Measurement results.

    Examples
    --------
    Specify multipole `degrees` (0, 0, 2) and bispectrum 'full'
    non-diagonal `form` with the first wavenumber fixed in bin 5,
    which override the parameters supplied by `paramset`.

    >>> results = compute_bispec(
    ...     catalogue_data, catalogue_rand,
    ...     degrees=(0, 0, 2),
    ...     form='full',
    ...     idx_bin=5,
    ...     paramset=paramset
    ... )

    See more analogous examples in
    :func:`~triumvirate.twopt.compute_powspec`.

    """
    # if logger:
    #     logger.info(
    #         "Measuring bispectrum from paired survey-type catalogues...",
    #         cpp_state='start'
    #     )

    results = _compute_3pt_stats_survey_like(
        _compute_bispec,
        catalogue_data, catalogue_rand,
        los_data=los_data, los_rand=los_rand,
        paramset=paramset, params_sampling=sampling_params,
        degrees=degrees, binning=binning, form=form, idx_bin=idx_bin,
        types={'catalogue_type': 'survey', 'statistic_type': 'bispec'},
        save=save, logger=logger
    )

    # if logger:
    #     logger.info(
    #         "... measured bispectrum from paired survey-type catalogues.",
    #         cpp_state='end'
    #     )

    return results


def compute_3pcf(catalogue_data, catalogue_rand,
                 los_data=None, los_rand=None,
                 degrees=None, binning=None, form=None, idx_bin=None,
                 sampling_params=None,
                 paramset=None,
                 save=False, logger=None):
    """Compute three-point correlation function from survey-like
    data and random catalogues in the local plane-parallel approximation.

    Parameters
    ----------
    catalogue_data : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Data-source catalogue.
    catalogue_rand : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Random-source catalogue.
    los_data : (N, 3) array of float, optional
        Specified lines of sight for the data-source catalogue.
        If `None` (default), this is automatically computed using
        :meth:`~triumvirate.catalogue.ParticleCatalogue.compute_los`.
    los_rand : (N, 3) array of float, optional
        Specified lines of sight for the random-source catalogue.
        If `None` (default), this is automatically computed using
        :meth:`~triumvirate.catalogue.ParticleCatalogue.compute_los`.
    degrees : tuple of (int, int, int) or str of length 3, optional
        Multipole degrees either as a tuple ('ell1', 'ell2', 'ELL') or
        as a string of length 3.  If not `None` (default), this will
        override ``paramset['degrees']`` entries.  If a string, multipole
        degrees are assumed to be single-digit integers.
    binning : :class:`~triumvirate.dataobjs.Binning`, optional
        Binning for the measurements.  If `None` (default),
        this is constructed from `paramset`.
    form : {'diag', 'full'}, optional
        Binning form of the measurements.  If not `None` (default),
        this will override ``paramset['form']``.
    idx_bin : int, optional
        Fixed bin index for the first coordinate dimension
        when binning `form` is 'full'.  If not `None` (default),
        this will override ``paramset['idx_bin']``.
    sampling_params : dict, optional
        Dictionary containing a subset of the following entries
        for sampling parameters---

        - 'alignment': {'centre', 'pad'};
        - 'boxsize': sequence of [float, float, float];
        - 'ngrid': sequence of [int, int, int];
        - 'assignment': {'ngp', 'cic', 'tsc', 'pcs'};
        - 'interlace': bool;

        and one and only one of the following when 'alignment' is 'pad'---

        - 'boxpad': float;
        - 'gridpad': float.

        This will override corresponding entries in `paramset`.
    paramset : :class:`~triumvirate.parameters.ParameterSet`, optional
        Full parameter set (default is `None`).  This is used in lieu of
        `degrees`, `binning`, `form`, `idx_bin` or `sampling_params`.
    save : {'.txt', '.npz', False}, optional
        If not `False` (default), save the measurements as a '.txt' file
        or in '.npz' format.  The save path is determined from `paramset`
        (if unset, a default file in the current working directory is
        used).
    logger : :class:`logging.Logger`, optional
        Logger (default is `None`).

    Returns
    -------
    results : dict of {str: :class:`numpy.ndarray`}
        Measurement results.

    Examples
    --------
    See analogous examples in :func:`~triumvirate.threept.compute_bispec`.

    """
    # if logger:
    #     logger.info(
    #         "Measuring three-point correlation function "
    #         "from paired survey-type catalogues...",
    #         cpp_state='start'
    #     )

    results = _compute_3pt_stats_survey_like(
        _compute_3pcf,
        catalogue_data, catalogue_rand,
        los_data=los_data, los_rand=los_rand,
        paramset=paramset, params_sampling=sampling_params,
        degrees=degrees, binning=binning, form=form, idx_bin=idx_bin,
        types={'catalogue_type': 'survey', 'statistic_type': '3pcf'},
        save=save, logger=logger
    )

    # if logger:
    #     logger.info(
    #         "... measured three-point correlation function "
    #         "from paired survey-type catalogues.",
    #         cpp_state='end'
    #     )

    return results

# def compute_bispec_for_los_choice(catalogue_data, catalogue_rand, los_choice,
#                                   los_data=None, los_rand=None,
#                                   degrees=None, binning=None,
#                                   form=None, idx_bin=None,
#                                   sampling_params=None,
#                                   paramset=None,
#                                   save=False, logger=None):
#     """Compute bispectrum from data and random catalogues.
#
#     Parameters
#     ----------
#     catalogue_data : :class:`~triumvirate.catalogue.ParticleCatalogue`
#         Data-source catalogue.
#     catalogue_rand : :class:`~triumvirate.catalogue.ParticleCatalogue`
#         Random-source catalogue.
#     los_choice : {0, 1, 2}, int
#         Line-of-sight choice.
#     los_data : (N, 3) array of float, optional
#         Specified lines of sight for the data-source catalogue.
#         If `None` (default), this is automatically computed using
#         :meth:`~triumvirate.catalogue.ParticleCatalogue.compute_los`.
#     los_rand : (N, 3) array of float, optional
#         Specified lines of sight for the random-source catalogue.
#         If `None` (default), this is automatically computed using
#         :meth:`~triumvirate.catalogue.ParticleCatalogue.compute_los`.
#     degrees : tuple of (int, int, int) or str of length 3, optional
#         Multipole degrees either as a tuple ('ell1', 'ell2', 'ELL') or
#         as a string of length 3.  If not `None` (default), this will
#         override ``paramset['degrees']`` entries.  If a string, multipole
#         degrees are assumed to be single-digit integers.
#     binning : :class:`~triumvirate.dataobjs.Binning`, optional
#         Binning for the measurements.  If `None` (default),
#         this is constructed from `paramset`.
#     form : {'diag', 'full'}, optional
#         Binning form of the measurements.  If not `None` (default),
#         this will override ``paramset['form']``.
#     idx_bin : int, optional
#         Fixed bin index for the first coordinate dimension
#         when binning `form` is 'full'.  If not `None` (default),
#         this will override ``paramset['idx_bin']``.
#     sampling_params : dict, optional
#         Dictionary containing a subset of the following entries
#         for sampling parameters---
#
#         - 'alignment': {'centre', 'pad'};
#         - 'boxsize': sequence of [float, float, float];
#         - 'ngrid': sequence of [int, int, int];
#         - 'assignment': {'ngp', 'cic', 'tsc', 'pcs'};
#         - 'interlace': bool;
#
#         and one and only one of the following when 'alignment' is 'pad'---
#
#         - 'boxpad': float;
#         - 'gridpad': float.
#
#         This will override corresponding entries in `paramset`.
#     paramset : :class:`~triumvirate.parameters.ParameterSet`, optional
#         Full parameter set (default is `None`).  This is used in lieu of
#         `degrees`, `binning`, `form`, `idx_bin` or `sampling_params`.
#     save : {'.txt', '.npz', False}, optional
#         If not `False` (default), save the measurements as a '.txt' file
#         or in '.npz' format.  The save path is determined from `paramset`
#         (if unset, a default file in the current working directory is
#         used).
#     logger : :class:`logging.Logger`, optional
#         Logger (default is `None`).
#
#     Returns
#     -------
#     results : dict of {str: :class:`numpy.ndarray`}
#         Measurement results.
#
#     """
#     # --------------------------------------------------------------------
#     # Initialisation
#     # --------------------------------------------------------------------

#     # -- Parameters ------------------------------------------------------

#     paramset = _amalgamate_parameters(
#         paramset=paramset, params_sampling=sampling_params,
#         degrees=degrees, binning=binning, form=form, idx_bin=idx_bin,
#         catalogue_type='survey', statistic_type='bispec'
#     )

#     if isinstance(paramset, dict):  # likely redundant but safe
#         paramset = ParameterSet(param_dict=paramset, logger=logger)

#     if logger:
#         logger.info("Parameter set have been initialised.")

#     # -- Data ------------------------------------------------------------

#     # Set up binning.
#     if binning is None:
#         binning = Binning.from_parameter_set(paramset)

#     if logger:
#         logger.info("Binning has been initialised.")

#     # Set up lines of sight.
#     if los_data is None:
#         los_data = catalogue_data.compute_los()
#     if los_rand is None:
#         los_rand = catalogue_rand.compute_los()

#     los_data = np.ascontiguousarray(los_data)
#     los_rand = np.ascontiguousarray(los_rand)

#     if logger:
#         logger.info("Lines of sight have been initialised.")

#     # Set up box alignment.
#     if paramset['alignment'] == 'centre':
#         catalogue_data.centre(
#             [paramset['boxsize'][ax] for ax in ['x', 'y', 'z']],
#             catalogue_ref=catalogue_rand
#         )
#     if paramset['alignment'] == 'pad':
#         if paramset['padscale'] == 'box':
#             kwargs = {'boxsize_pad': paramset['padfactor']}
#         if paramset['padscale'] == 'grid':
#             kwargs = {
#                 'ngrid': [paramset['ngrid'][ax] for ax in ['x', 'y', 'z']],
#                 'ngrid_pad': paramset['padfactor']
#             }
#         catalogue_data.pad(
#             [paramset['boxsize'][ax] for ax in ['x', 'y', 'z']],
#             catalogue_ref=catalogue_rand, **kwargs
#         )

#     if logger:
#         logger.info("Catalogues have been aligned.")

#     # --------------------------------------------------------------------
#     # Measurements
#     # --------------------------------------------------------------------

#     # Prepare catalogues.
#     if logger:
#         logger.info(
#             "Preparing catalogue for clustering algorithm...",
#             cpp_state='start'
#         )

#     particles_data = \
#         catalogue_data._convert_to_cpp_catalogue(verbose=paramset['verbose'])
#     particles_rand = \
#         catalogue_rand._convert_to_cpp_catalogue(verbose=paramset['verbose'])

#     if logger:
#         logger.info(
#             "... prepared catalogue for clustering algorithm.",
#             cpp_state='end'
#         )

#     # Set up constants.
#     alpha = catalogue_data.wtotal / catalogue_rand.wtotal

#     if logger:
#         logger.info("Alpha contrast: %.6e.", alpha)

#     if paramset['norm_convention'] == 'particle':
#         norm_factor = _calc_bispec_normalisation_from_particles(
#             particles_rand, alpha
#         )
#         norm_factor_alt = _calc_bispec_normalisation_from_mesh(
#             particles_rand, paramset, alpha
#         )
#     if paramset['norm_convention'] == 'mesh':
#         norm_factor = _calc_bispec_normalisation_from_mesh(
#             particles_rand, paramset, alpha
#         )
#         norm_factor_alt = _calc_bispec_normalisation_from_particles(
#             particles_rand, alpha
#         )

#     if logger:
#         logger.info(
#             "Normalisation factors: %.6e (used), %.6e (alternative).",
#             norm_factor, norm_factor_alt
#         )

#     # Perform measurement.
#     if logger:
#         logger.info("Measuring clustering statistics...", cpp_state='start')

#     results = _compute_bispec_for_los_choice(
#         particles_data, particles_rand, los_data, los_rand, los_choice,
#         paramset, binning, norm_factor
#     )

#     if logger:
#         logger.info("... measured clustering statistics.", cpp_state='end')

#     if save:
#         odirpath = paramset['directories']['measurements'] or ""
#         header = "\n".join([
#             catalogue_data.write_attrs_as_header(catalogue_ref=catalogue_rand),  # noqa: E501
#             _print_measurement_header(paramset, norm_factor, norm_factor_alt),  # noqa: E501
#         ])
#         if save.lower() == '.txt':
#             datatab = _assemble_measurement_datatab(results, paramset)
#             datafmt = '\t'.join(
#                 ['%.9e'] * 4 + ['%10d'] + ['% .9e'] * (datatab.shape[-1] - 5)
#             )
#             ofilename = _get_measurement_filename(paramset)
#             ofilepath = Path(odirpath, ofilename).with_suffix('.txt')
#             np.savetxt(
#                 ofilepath, datatab,
#                 fmt=datafmt, header=header, delimiter='\t'
#             )
#         elif save.lower().endswith('.npz'):
#             results.update({'header': header})
#             ofilename = _get_measurement_filename(paramset)
#             ofilepath = Path(odirpath, ofilename).with_suffix('.npz')
#             np.savez(ofilepath, **results)
#         else:
#             raise ValueError(
#                 f"Unrecognised save format for measurements: {save}."
#             )

#         if logger:
#             logger.info("Measurements saved to %s.", ofilepath)

#     return results


# ========================================================================
# Simulation statistics
# ========================================================================

def _compute_3pt_stats_sim_like(threept_algofunc, catalogue_data,
                                paramset=None, params_sampling=None,
                                degrees=None, binning=None, types=None,
                                form=None, idx_bin=None,
                                save=False, logger=None):
    """Compute three-point statistics from a simulation-box catalogue
    in the global plane-parallel approximation.

    Parameters
    ----------
    threept_algofunc : callable
        Three-point statistic algorithmic function.
    catalogue_data : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Data-source catalogue.
    paramset : :class:`~triumvirate.parameters.ParameterSet`, optional
        Full parameter set.  If `None` (default), `degrees`, `binning`,
        `form` and `params_sampling` should be provided; `idx_bin`
        should be provided when `form` is 'full'.
    params_sampling : dict, optional
        Dictionary containing a subset of the following entries
        for sampling parameters---

        - 'alignment': {'centre', 'pad'};
        - 'boxsize': sequence of [float, float, float];
        - 'ngrid': sequence of [int, int, int];
        - 'assignment': {'ngp', 'cic', 'tsc', 'pcs'};
        - 'interlace': bool;

        and one and only one of the following when 'alignment' is 'pad'---

        - 'boxpad': float;
        - 'gridpad': float.

        This will override corresponding entries in `paramset`.
    degrees : tuple of (int, int, int) or str of length 3, optional
        Multipole degrees either as a tuple ('ell1', 'ell2', 'ELL') or
        as a string of length 3.  If not `None` (default), this will
        override ``paramset['degrees']`` entries.  If a string, multipole
        degrees are assumed to be single-digit integers.
    binning : :class:`~triumvirate.dataobjs.Binning`, optional
        Binning for the measurements.  If `None` (default),
        this is constructed from `paramset`.
    types : dict, optional
        `catalogue_type` and `statistic_type` (default is `None`).
        This should be set by the caller of this function.
    form : {'diag', 'full'}, optional
        Binning form of the measurements.  If not `None` (default),
        this will override ``paramset['form']``.
    idx_bin : int, optional
        Fixed bin index for the first coordinate dimension
        when binning `form` is 'full'.  If not `None` (default),
        this will override ``paramset['idx_bin']``.
    save : {'.txt', '.npz', False}, optional
        If not `False` (default), save the measurements as a '.txt' file
        or in '.npz' format.  The save path is determined from `paramset`
        (if unset, a default file in the current working directory is
        used).
    logger : :class:`logging.Logger`, optional
        Logger (default is `None`).

    Returns
    -------
    results : dict of {str: :class:`numpy.ndarray`}
        Measurement results.

    """
    # --------------------------------------------------------------------
    # Initialisation
    # --------------------------------------------------------------------

    # -- Parameters ------------------------------------------------------

    paramset = _amalgamate_parameters(
        paramset=paramset, params_sampling=params_sampling,
        degrees=degrees, binning=binning, form=form, idx_bin=idx_bin, **types
    )

    if isinstance(paramset, dict):  # likely redundant but safe
        paramset = ParameterSet(param_dict=paramset, logger=logger)

    if logger:
        logger.info("Parameter set have been initialised.")

    # -- Data ------------------------------------------------------------

    # Set up binning.
    if binning is None:
        binning = Binning.from_parameter_set(paramset)

    if logger:
        logger.info("Binning has been initialised.")

    # Set up box alignment.
    catalogue_data.periodise(
        [paramset['boxsize'][ax] for ax in ['x', 'y', 'z']]
    )

    if logger:
        logger.info("Catalogue box has been periodised.")

    # --------------------------------------------------------------------
    # Measurements
    # --------------------------------------------------------------------

    # Prepare catalogues.
    if not catalogue_data['nz'].any():
        catalogue_data.compute_mean_density(
            boxsize=list(paramset['boxsize'].values())
        )
        if logger:
            logger.info(
                "Inserted missing 'nz' field "
                "based on particle count and boxsize."
            )

    if logger:
        logger.info(
            "Preparing catalogue for clustering algorithm...",
            cpp_state='start'
        )

    particles_data = \
        catalogue_data._convert_to_cpp_catalogue(verbose=paramset['verbose'])

    if logger:
        logger.info(
            "... prepared catalogue for clustering algorithm.",
            cpp_state='end'
        )

    # Set up constants.
    if paramset['norm_convention'] == 'particle':
        norm_factor = _calc_bispec_normalisation_from_particles(
            particles_data, alpha=1.
        )
        norm_factor_alt = _calc_bispec_normalisation_from_mesh(
            particles_data, paramset, alpha=1.
        )
    if paramset['norm_convention'] == 'mesh':
        norm_factor = _calc_bispec_normalisation_from_mesh(
            particles_data, paramset, alpha=1.
        )
        norm_factor_alt = _calc_bispec_normalisation_from_particles(
            particles_data, alpha=1.
        )

    if logger:
        logger.info(
            "Normalisation factors: %.6e (used), %.6e (alternative).",
            norm_factor, norm_factor_alt
        )

    # Perform measurement.
    if logger:
        logger.info("Measuring clustering statistics...", cpp_state='start')

    results = threept_algofunc(particles_data, paramset, binning, norm_factor)

    if logger:
        logger.info("... measured clustering statistics.", cpp_state='end')

    if save:
        odirpath = paramset['directories']['measurements'] or ""
        header = "\n".join([
            catalogue_data.write_attrs_as_header(),
            _print_measurement_header(paramset, norm_factor, norm_factor_alt),
        ])
        if save.lower() == '.txt':
            datatab = _assemble_measurement_datatab(results, paramset)
            datafmt = '\t'.join(
                ['%.9e'] * 4 + ['%10d'] + ['% .9e'] * (datatab.shape[-1] - 5)
            )
            ofilename = _get_measurement_filename(paramset)
            ofilepath = Path(odirpath, ofilename).with_suffix('.txt')
            np.savetxt(
                ofilepath, datatab, fmt=datafmt, header=header, delimiter='\t'
            )
        elif save.lower().endswith('.npz'):
            results.update({'header': header})
            ofilename = _get_measurement_filename(paramset)
            ofilepath = Path(odirpath, ofilename).with_suffix('.npz')
            np.savez(ofilepath, **results)
        else:
            raise ValueError(
                f"Unrecognised save format for measurements: {save}."
            )

        if logger:
            logger.info("Measurements saved to %s.", ofilepath)

    return results


def compute_bispec_in_gpp_box(catalogue_data,
                              degrees=None, binning=None,
                              form=None, idx_bin=None,
                              sampling_params=None,
                              paramset=None,
                              save=False, logger=None):
    """Compute bispectrum from a simulation-box catalogue
    in the global plane-parallel approximation.

    Parameters
    ----------
    catalogue_data : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Data-source catalogue.
    degrees : tuple of (int, int, int) or str of length 3, optional
        Multipole degrees either as a tuple ('ell1', 'ell2', 'ELL') or
        as a string of length 3.  If not `None` (default), this will
        override ``paramset['degrees']`` entries.  If a string, multipole
        degrees are assumed to be single-digit integers.
    binning : :class:`~triumvirate.dataobjs.Binning`, optional
        Binning for the measurements.  If `None` (default),
        this is constructed from `paramset`.
    form : {'diag', 'full'}, optional
        Binning form of the measurements.  If not `None` (default),
        this will override ``paramset['form']``.
    idx_bin : int, optional
        Fixed bin index for the first coordinate dimension
        when binning `form` is 'full'.  If not `None` (default),
        this will override ``paramset['idx_bin']``.
    sampling_params : dict, optional
        Dictionary containing a subset of the following entries
        for sampling parameters---

        - 'alignment': {'centre', 'pad'};
        - 'boxsize': sequence of [float, float, float];
        - 'ngrid': sequence of [int, int, int];
        - 'assignment': {'ngp', 'cic', 'tsc', 'pcs'};
        - 'interlace': bool;

        and one and only one of the following when 'alignment' is 'pad'---

        - 'boxpad': float;
        - 'gridpad': float.

        This will override corresponding entries in `paramset`.
    paramset : :class:`~triumvirate.parameters.ParameterSet`, optional
        Full parameter set (default is `None`).  This is used in lieu of
        `degrees`, `binning`, `form`, `idx_bin` or `sampling_params`.
    save : {'.txt', '.npz', False}, optional
        If not `False` (default), save the measurements as a '.txt' file
        or in '.npz' format.  The save path is determined from `paramset`
        (if unset, a default file in the current working directory is
        used).
    logger : :class:`logging.Logger`, optional
        Logger (default is `None`).

    Returns
    -------
    results : dict of {str: :class:`numpy.ndarray`}
        Measurement results.

    Examples
    --------
    See analogous examples in :func:`~triumvirate.threept.compute_bispec`
    (though without line-of-sight arguments).

    """
    # if logger:
    #     logger.info(
    #         "Measuring bispectrum from a simulation-box catalogue "
    #         "in the global plane-parallel approximation...",
    #         cpp_state='start'
    #     )

    results = _compute_3pt_stats_sim_like(
        _compute_bispec_in_gpp_box,
        catalogue_data,
        paramset=paramset, params_sampling=sampling_params,
        degrees=degrees, binning=binning, form=form, idx_bin=idx_bin,
        types={'catalogue_type': 'sim', 'statistic_type': 'bispec'},
        save=save, logger=logger
    )

    # if logger:
    #     logger.info(
    #         "... measured bispectrum from a simulation-box catalogue "
    #         "in the global plane-parallel approximation.",
    #         cpp_state='end'
    #     )

    return results


def compute_3pcf_in_gpp_box(catalogue_data,
                            degrees=None, binning=None,
                            form=None, idx_bin=None,
                            sampling_params=None,
                            paramset=None,
                            save=False, logger=None):
    """Compute three-point correlation function from a simulation-box
    catalogue in the global plane-parallel approximation.

    Parameters
    ----------
    catalogue_data : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Data-source catalogue.
    degrees : tuple of (int, int, int) or str of length 3, optional
        Multipole degrees either as a tuple ('ell1', 'ell2', 'ELL') or
        as a string of length 3.  If not `None` (default), this will
        override ``paramset['degrees']`` entries.  If a string, multipole
        degrees are assumed to be single-digit integers.
    binning : :class:`~triumvirate.dataobjs.Binning`, optional
        Binning for the measurements.  If `None` (default),
        this is constructed from `paramset`.
    form : {'diag', 'full'}, optional
        Binning form of the measurements.  If not `None` (default),
        this will override ``paramset['form']``.
    idx_bin : int, optional
        Fixed bin index for the first coordinate dimension
        when binning `form` is 'full'.  If not `None` (default),
        this will override ``paramset['idx_bin']``.
    sampling_params : dict, optional
        Dictionary containing a subset of the following entries
        for sampling parameters---

        - 'alignment': {'centre', 'pad'};
        - 'boxsize': sequence of [float, float, float];
        - 'ngrid': sequence of [int, int, int];
        - 'assignment': {'ngp', 'cic', 'tsc', 'pcs'};
        - 'interlace': bool;

        and one and only one of the following when 'alignment' is 'pad'---

        - 'boxpad': float;
        - 'gridpad': float.

        This will override corresponding entries in `paramset`.
    paramset : :class:`~triumvirate.parameters.ParameterSet`, optional
        Full parameter set (default is `None`).  This is used in lieu of
        `degree`, `binning`, `form`, `idx_bin` or `sampling_params`.
    save : {'.txt', '.npz', False}, optional
        If not `False` (default), save the measurements as a '.txt' file
        or in '.npz' format.  The save path is determined from `paramset`
        (if unset, a default file in the current working directory is
        used).
    logger : :class:`logging.Logger`, optional
        Logger (default is `None`).

    Returns
    -------
    results : dict of {str: :class:`numpy.ndarray`}
        Measurement results.

    Examples
    --------
    See analogous examples in :func:`~triumvirate.threept.compute_bispec`
    (though without line-of-sight arguments).

    """
    # if logger:
    #     logger.info(
    #         "Measuring three-point correlation function "
    #         "from a simulation-box catalogue "
    #         "in the global plane-parallel approximation...",
    #         cpp_state='start'
    #     )

    results = _compute_3pt_stats_sim_like(
        _compute_3pcf_in_gpp_box,
        catalogue_data,
        paramset=paramset, params_sampling=sampling_params,
        degrees=degrees, binning=binning, form=form, idx_bin=idx_bin,
        types={'catalogue_type': 'sim', 'statistic_type': '3pcf'},
        save=save, logger=logger
    )

    # if logger:
    #     logger.info(
    #         "... measured three-point correlation function "
    #         "from a simulation-box catalogue "
    #         "in the global plane-parallel approximation.",
    #         cpp_state='end'
    #     )

    return results


# ========================================================================
# Window statistics
# ========================================================================

def compute_3pcf_window(catalogue_rand, los_rand=None,
                        degrees=None, wa_orders=None,
                        binning=None, form=None, idx_bin=None,
                        sampling_params=None,
                        paramset=None,
                        save=False, logger=None):
    """Compute three-point correlation function window
    from a random catalogue.

    Parameters
    ----------
    catalogue_rand : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Random-source catalogue.
    los_rand : (N, 3) array of float, optional
        Specified lines of sight for the random-source catalogue.
        If `None` (default), this is automatically computed using
        :meth:`~triumvirate.catalogue.ParticleCatalogue.compute_los`.
    degrees : tuple of (int, int, int) or str of length 3, optional
        Multipole degrees either as a tuple ('ell1', 'ell2', 'ELL') or
        as a string of length 3.  If not `None` (default), this will
        override ``paramset['degrees']`` entries.  If a string, multipole
        degrees are assumed to be single-digit integers.
    wa_orders : tuple of (int, int) or str of length 2, optional
        Wide-angle correction orders either as a tuple ('i_wa', 'j_wa') or
        as a string of length 2.  If not `None` (default), this will
        override ``paramset['wa_orders']`` entries.  If a string,
        multipole degrees are assumed to be single-digit integers.
    binning : :class:`~triumvirate.dataobjs.Binning`, optional
        Binning for the measurements.  If `None` (default),
        this is constructed from `paramset`.
    form : {'diag', 'full'}, optional
        Binning form of the measurements.  If not `None` (default),
        this will override ``paramset['form']``.
    idx_bin : int, optional
        Fixed bin index for the first coordinate dimension
        when binning `form` is 'full'.  If not `None` (default),
        this will override ``paramset['idx_bin']``.
    sampling_params : dict, optional
        Dictionary containing a subset of the following entries
        for sampling parameters---

        - 'alignment': {'centre', 'pad'};
        - 'boxsize': sequence of [float, float, float];
        - 'ngrid': sequence of [int, int, int];
        - 'assignment': {'ngp', 'cic', 'tsc', 'pcs'};
        - 'interlace': bool;

        and one and only one of the following when 'alignment' is 'pad'---

        - 'boxpad': float;
        - 'gridpad': float.

        This will override corresponding entries in `paramset`.
    paramset : :class:`~triumvirate.parameters.ParameterSet`, optional
        Full parameter set (default is `None`).  This is used in lieu of
        `degrees`, `binning`, `form`, `idx_bin` or `sampling_params`.
    save : {'.txt', '.npz', False}, optional
        If not `False` (default), save the measurements as a '.txt' file
        or in '.npz' format.  The save path is determined from `paramset`
        (if unset, a default file in the current working directory is
        used).
    logger : :class:`logging.Logger`, optional
        Logger (default is `None`).

    Returns
    -------
    results : dict of {str: :class:`numpy.ndarray`}
        Measurement results.

    Examples
    --------
    See analogous examples in :func:`~triumvirate.threept.compute_bispec`.

    """
    # --------------------------------------------------------------------
    # Initialisation
    # --------------------------------------------------------------------

    # -- Parameters ------------------------------------------------------

    if wa_orders is not None:
        wide_angle = True
    else:
        wide_angle = False

    paramset = _amalgamate_parameters(
        paramset=paramset, params_sampling=sampling_params,
        degrees=degrees, wa_orders=wa_orders,
        binning=binning, form=form, idx_bin=idx_bin,
        catalogue_type='random', statistic_type='3pcf-win'
    )

    if isinstance(paramset, dict):  # likely redundant but safe
        paramset = ParameterSet(param_dict=paramset, logger=logger)

    if logger:
        logger.info("Parameter set have been initialised.")

    # -- Data ------------------------------------------------------------

    # Set up binning.
    if binning is None:
        binning = Binning.from_parameter_set(paramset)

    if logger:
        logger.info("Binning has been initialised.")

    # Set up lines of sight.
    if los_rand is None:
        los_rand = catalogue_rand.compute_los()
    los_rand = np.ascontiguousarray(los_rand)

    if logger:
        logger.info("Lines of sight have been initialised.")

    # Set up box alignment.
    catalogue_rand.centre(
        [paramset['boxsize'][ax] for ax in ['x', 'y', 'z']]
    )

    if logger:
        logger.info("Catalogues have been aligned.")

    # --------------------------------------------------------------------
    # Measurements
    # --------------------------------------------------------------------

    # Prepare catalogues.
    particles_rand = \
        catalogue_rand._convert_to_cpp_catalogue(verbose=paramset['verbose'])

    # Set up constants.
    if paramset['norm_convention'] == 'particle':
        norm_factor = _calc_bispec_normalisation_from_particles(
            particles_rand, alpha=1.
        )
        norm_factor_alt = _calc_bispec_normalisation_from_mesh(
            particles_rand, paramset, alpha=1.
        )
    if paramset['norm_convention'] == 'mesh':
        norm_factor = _calc_bispec_normalisation_from_mesh(
            particles_rand, paramset, alpha=1.
        )
        norm_factor_alt = _calc_bispec_normalisation_from_particles(
            particles_rand, alpha=1.
        )

    if logger:
        logger.info(
            "Normalisation factors: %.6e (used), %.6e (alternative).",
            norm_factor, norm_factor_alt
        )

    # Perform measurement.
    if logger:
        logger.info(
            "Measuring window function statistics...", cpp_state='start'
        )

    results = _compute_3pcf_window(
        particles_rand, los_rand,
        paramset, binning, alpha=1., norm_factor=norm_factor,
        wide_angle=wide_angle
    )

    if logger:
        logger.info(
            "... measured window function statistics.", cpp_state='end'
        )

    if save:
        odirpath = paramset['directories']['measurements'] or ""
        header = "\n".join([
            catalogue_rand.write_attrs_as_header(),
            _print_measurement_header(paramset, norm_factor, norm_factor_alt),
        ])
        if save.lower() == '.txt':
            datatab = _assemble_measurement_datatab(results, paramset)
            datafmt = '\t'.join(
                ['%.9e'] * 4 + ['%10d'] + ['% .9e'] * (datatab.shape[-1] - 5)
            )
            ofilename = _get_measurement_filename(paramset)
            ofilepath = Path(odirpath, ofilename).with_suffix('.txt')
            np.savetxt(
                ofilepath, datatab, fmt=datafmt, header=header, delimiter='\t'
            )
        elif save.lower().endswith('.npz'):
            results.update({'header': header})
            ofilename = _get_measurement_filename(paramset)
            ofilepath = Path(odirpath, ofilename).with_suffix('.npz')
            np.savez(ofilepath, **results)
        else:
            raise ValueError(
                f"Unrecognised save format for measurements: {save}."
            )

        if logger:
            logger.info("Measurements saved to %s.", ofilepath)

    return results
