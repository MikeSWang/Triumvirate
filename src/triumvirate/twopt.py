"""
Two-Point Clustering Statistics (:mod:`~triumvirate.twopt`)
==========================================================================

Measure two-point clustering statistics from catalogues.

Local plane-parallel estimators
-------------------------------

.. autosummary::
    compute_powspec
    compute_corrfunc

Global plane-parallel estimators
--------------------------------

.. autosummary::
    compute_powspec_in_gpp_box
    compute_corrfunc_in_gpp_box

Window function estimator
-------------------------

.. autosummary::
    compute_corrfunc_window

"""
import warnings
from pathlib import Path

import numpy as np

from triumvirate._twopt import (
    _calc_powspec_normalisation_from_mesh,
    _calc_powspec_normalisation_from_particles,
    _compute_corrfunc,
    _compute_corrfunc_in_gpp_box,
    _compute_corrfunc_window,
    _compute_powspec,
    _compute_powspec_in_gpp_box,
)
from triumvirate.dataobjs import Binning
from triumvirate.parameters import (
    _modify_measurement_parameters,
    _modify_sampling_parameters,
    fetch_paramset_template,
    ParameterSet,
)


def _amalgamate_parameters(paramset=None, params_sampling=None,
                           degree=None, binning=None, **type_kwargs):
    """Amalgamate a parameter set with overriding sampling parameters
    and the measured multipole degree and coordinate binning.

    Parameters
    ----------
    paramset : :class:`~triumvirate.parameters.ParameterSet`, optional
        Full parameter set.  If `None` (default), the other parameters
        must be provided.
    params_sampling : dict, optional
        Dictionary containing a subset of the following entries
        for sampling parameters---

        - 'boxalign': {'centre', 'pad'};
        - 'boxsize': sequence of [float, float, float];
        - 'ngrid': sequence of [int, int, int];
        - 'assignment': {'ngp', 'cic', 'tsc', 'pcs'};
        - 'interlace': bool;

        and one and only one of the following when 'boxalign' is 'pad'---

        - 'boxpad': float;
        - 'gridpad': float.

        This will override corresponding entries in `paramset`.
    degree : int, optional
        Multipole degree.  If not `None` (default), this will override
        ``paramset['degrees']['ELL']``.
    binning : :class:`~.triumvirate.dataobjs.Binning`, optional
        Binning (default is `None`).
    **type_kwargs
        `catalogue_type` and `statistic_type` parameters to be filled in.

    Returns
    -------
    paramset : :class:`~triumvirate.parameters.ParameterSet`, optional
        Amalgamated full parameter set.

    Raises
    ------
    ValueError
        When `paramset` is `None` while `params_sampling`, `degree` or
        `binning` is also `None`.

    """
    if paramset is None and params_sampling is None:
        raise ValueError(
            "Either `paramset` or `sampling_params` must be provided."
        )
    if paramset is None and degree is None:
        raise ValueError("Either `paramset` or `degree` must be provided.")
    if paramset is None and binning is None:
        raise ValueError("Either `paramset` or `binning` must be provided.")

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
    if degree is not None:
        params_measure.update({'ELL': degree})
    if binning is not None:
        params_measure.update({
            'bin_min': binning.bin_min,
            'bin_max': binning.bin_max,
            'num_bins': binning.num_bins,
        })

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
        two-point statistics.

    """
    multipole = paramset['degrees']['ELL']

    output_tag = paramset['tags']['output'] or ""

    if paramset['statistic_type'] == 'powspec':
        return "pk{:d}{}".format(multipole, output_tag)
    if paramset['statistic_type'] == '2pcf':
        return "xi{:d}{}".format(multipole, output_tag)
    if paramset['statistic_type'] == '2pcf-win':
        return "xiw{:d}{}".format(multipole, output_tag)

    raise ValueError(
        "`paramset` 'statistic_type' does not correspond to a "
        "recognised two-point statistic."
    )


def _print_measurement_header(paramset, norm_factor, norm_factor_alt):
    """Print two-point statistic measurement header including
    sampling parameters, normalisation factors and data table columns.

    Parameters
    ----------
    paramset : :class:`~triumvirate.parameters.ParameterSet`
        Parameter set.
    norm_factor, norm_factor_alt : double
        Normalisation factor used and the alternative.

    Returns
    -------
    text_header : str
        Measurement information as a header string.

    Raises
    ------
    ValueError
        When `paramset` indicates the measurements are not
        two-point statistics.

    """
    if paramset['norm_convention'] == 'particle':
        norm_convention = 'particle'
        norm_convention_alt = 'mesh'
    if paramset['norm_convention'] == 'mesh':
        norm_convention = 'mesh'
        norm_convention_alt = 'particle'

    if paramset['npoint'] != '2pt':
        raise ValueError(
            "`paramset` 'statistic_type' does not correspond to a "
            "recognised two-point statistic."
        )
    if paramset['space'] == 'fourier':
        datatab_colnames = [
            "k_cen", "k_eff", "nmodes",
            "Re{{pk{:d}_raw}}".format(paramset['degrees']['ELL']),
            "Im{{pk{:d}_raw}}".format(paramset['degrees']['ELL']),
            "Re{{pk{:d}_shot}}".format(paramset['degrees']['ELL']),
            "Im{{pk{:d}_shot}}".format(paramset['degrees']['ELL'])
        ]
    if paramset['space'] == 'config':
        datatab_colnames = [
            "r_cen", "r_eff", "npairs",
            "Re{{xi{:d}}}".format(paramset['degrees']['ELL']),
            "Im{{xi{:d}}}".format(paramset['degrees']['ELL'])
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
        two-point statistics.

    """
    if paramset['npoint'] != '2pt':
        raise ValueError(
            "Measurement header being printed is for two-point statistics."
        )
    if paramset['space'] == 'fourier':
        datatab = np.transpose([
            measurements['kbin'], measurements['keff'], measurements['nmodes'],
            measurements['pk_raw'].real, measurements['pk_raw'].imag,
            measurements['pk_shot'].real, measurements['pk_shot'].imag,
        ])
    if paramset['space'] == 'config':
        datatab = np.transpose([
            measurements['rbin'], measurements['reff'], measurements['npairs'],
            measurements['xi'].real, measurements['xi'].imag,
        ])

    return datatab


# ========================================================================
# Survey statistics
# ========================================================================

def _compute_2pt_stats_survey_like(twopt_algofunc,
                                   catalogue_data, catalogue_rand,
                                   los_data=None, los_rand=None,
                                   paramset=None, params_sampling=None,
                                   degree=None, binning=None, types=None,
                                   save=False, logger=None):
    """Compute two-point statistics from survey-like data and random
    catalogues in the local plane-parallel approximation.

    Parameters
    ----------
    twopt_algofunc : callable
        Two-point statistic algorithmic function.
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
        Full parameter set.  If `None` (default), `degree`, `binning`
        and `params_sampling` must be provided.
    params_sampling : dict, optional
        Dictionary containing a subset of the following entries
        for sampling parameters---

        - 'boxalign': {'centre', 'pad'};
        - 'boxsize': sequence of [float, float, float];
        - 'ngrid': sequence of [int, int, int];
        - 'assignment': {'ngp', 'cic', 'tsc', 'pcs'};
        - 'interlace': bool;

        and one and only one of the following when 'boxalign' is 'pad'---

        - 'boxpad': float;
        - 'gridpad': float.

        This will override corresponding entries in `paramset`.
    degree : int, optional
        Multipole degree.  If not `None` (default), this will override
        ``paramset['degrees']['ELL']``.
    binning : :class:`~triumvirate.dataobjs.Binning`, optional
        Binning for the measurements.  If `None` (default), this is
        constructed from `paramset`.
    types : dict, optional
        `catalogue_type` and `statistic_type` (default is `None`).
        This should be set by the caller of this function.
    save : {'.txt', '.npz', False}, optional
        If not `False` (default), save the measurements as a '.txt' file
        or in '.npz' format. The save path is determined from `paramset`
        (if unset, a default file in the current working directory is
        used).
    logger : :class:`logging.Logger`, optional
        Logger (default is `None`).

    Returns
    -------
    results : dict of {str: :class:`numpy.ndarray`}
        Measurement results.

    Raises
    ------
    ValueError
        When `paramset` is `None` while `params_sampling`, `degree` or
        `binning` is also `None`.

    """
    # --------------------------------------------------------------------
    # Initialisation
    # --------------------------------------------------------------------

    # -- Parameters ------------------------------------------------------

    paramset = _amalgamate_parameters(
        paramset=paramset, params_sampling=params_sampling,
        degree=degree, binning=binning, **types
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
        norm_factor = _calc_powspec_normalisation_from_particles(
            particles_rand, alpha
        )
        norm_factor_alt = _calc_powspec_normalisation_from_mesh(
            particles_rand, paramset, alpha
        )
    if paramset['norm_convention'] == 'mesh':
        norm_factor = _calc_powspec_normalisation_from_mesh(
            particles_rand, paramset, alpha
        )
        norm_factor_alt = _calc_powspec_normalisation_from_particles(
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

    results = twopt_algofunc(
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
                ['%.9e'] * 2 + ['%10d'] + ['% .9e'] * (datatab.shape[-1] - 3)
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


def compute_powspec(catalogue_data, catalogue_rand,
                    los_data=None, los_rand=None,
                    degree=None, binning=None, sampling_params=None,
                    paramset=None,
                    save=False, logger=None):
    """Compute power spectrum from survey-like data and random catalogues
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
    degree : int, optional
        Multipole degree.  If not `None` (default), this will override
        ``paramset['degrees']['ELL']``.
    binning : :class:`~triumvirate.dataobjs.Binning`, optional
        Binning for the measurements.  If `None` (default), this is
        constructed from `paramset`.
    sampling_params : dict, optional
        Dictionary containing a subset of the following entries
        for sampling parameters---

        - 'boxalign': {'centre', 'pad'};
        - 'boxsize': sequence of [float, float, float];
        - 'ngrid': sequence of [int, int, int];
        - 'assignment': {'ngp', 'cic', 'tsc', 'pcs'};
        - 'interlace': bool;

        and one and only one of the following when 'boxalign' is 'pad'---

        - 'boxpad': float;
        - 'gridpad': float.

        This will override corresponding entries in `paramset`.
    paramset : :class:`~triumvirate.parameters.ParameterSet`, optional
        Full parameter set (default is `None`).  This is used in lieu of
        `degree`, `binning` or `sampling_params`.
    save : {'.txt', '.npz', False}, optional
        If not `False` (default), save the measurements as a '.txt' file
        or in '.npz' format. The save path is determined from `paramset`
        (if unset, a default file in the current working directory is
        used).
    logger : :class:`logging.Logger`, optional
        Logger (default is `None`).

    Returns
    -------
    results : dict of {str: :class:`numpy.ndarray`}
        Measurement results.

    Raises
    ------
    ValueError
        When `paramset` is `None` but `degree`, `binning` or
        `sampling_params` is also `None`.

    Examples
    --------
    Specify line-of-sight vectors explicitly.

    >>> results = compute_powspec(
    ...     catalogue_data, catalogue_rand,
    ...     los_data=np.ones((1e3, 3)), los_rand=np.ones((1e4, 3)),
    ...     paramset=None
    ... )

    Specify multipole `degree` 2, customised
    :class:`~triumvirate.dataobjs.Binning` object ``binning``.
    Whether `paramset` provided or not, relevant parameters are overriden
    by the supplied keyword arguments.

    >>> results = compute_powspec(
    ...     catalogue_data, catalogue_rand,
    ...     degree=2,
    ...     binning=binning,
    ...     sampling_params={
    ...         'boxsize': [1000., 1500., 1000.],
    ...         'ngrid': [256, 256, 256],
    ...         # 'boxalign' at default initial value in `ParameterSet`
    ...         # 'assignment' at default initial value in `ParameterSet`
    ...     }
    ... )

    """
    # if logger:
    #     logger.info(
    #         "Measuring power spectrum from paired survey-type catalogues...",
    #         cpp_state='start'
    #     )

    results = _compute_2pt_stats_survey_like(
        _compute_powspec,
        catalogue_data, catalogue_rand,
        los_data=los_data, los_rand=los_rand,
        paramset=paramset, params_sampling=sampling_params,
        degree=degree, binning=binning,
        types={'catalogue_type': 'survey', 'statistic_type': 'powspec'},
        save=save, logger=logger
    )

    # if logger:
    #     logger.info(
    #         "... measured power spectrum from paired survey-type catalogues.",  # noqa: E501
    #         cpp_state='end'
    #     )

    return results


def compute_corrfunc(catalogue_data, catalogue_rand,
                     los_data=None, los_rand=None,
                     degree=None, binning=None, sampling_params=None,
                     paramset=None,
                     save=False, logger=None):
    """Compute correlation function from survey-like data and random
    catalogues in the local plane-parallel approximation.

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
    degree : int, optional
        Multipole degree.  If not `None` (default), this will override
        ``paramset['degrees']['ELL']``.
    binning : :class:`~triumvirate.dataobjs.Binning`, optional
        Binning for the measurements.  If `None` (default), this is
        constructed from `paramset`.
    sampling_params : dict, optional
        Dictionary containing a subset of the following entries
        for sampling parameters---

        - 'boxalign': {'centre', 'pad'};
        - 'boxsize': sequence of [float, float, float];
        - 'ngrid': sequence of [int, int, int];
        - 'assignment': {'ngp', 'cic', 'tsc', 'pcs'};
        - 'interlace': bool;

        and one and only one of the following when 'boxalign' is 'pad'---

        - 'boxpad': float;
        - 'gridpad': float.

        This will override corresponding entries in `paramset`.
    paramset : :class:`~triumvirate.parameters.ParameterSet`, optional
        Full parameter set (default is `None`).  This is used
        in lieu of `degree`, `binning` or `sampling_params`.
    save : {'.txt', '.npz', False}, optional
        If not `False` (default), save the measurements as a '.txt' file
        or in '.npz' format. The save path is determined from `paramset`
        (if unset, a default file in the current working directory is
        used).
    logger : :class:`logging.Logger`, optional
        Logger (default is `None`).

    Returns
    -------
    results : dict of {str: :class:`numpy.ndarray`}
        Measurement results.

    Raises
    ------
    ValueError
        When `paramset` is `None` but `degree`, `binning` or
        `sampling_params` is also `None`.

    Examples
    --------
    See analogous examples in :func:`~triumvirate.twopt.compute_powspec`.

    """
    # if logger:
    #     logger.info(
    #         "Measuring two-point correlation function "
    #         "from paired survey-type catalogues...",
    #         cpp_state='start'
    #     )

    results = _compute_2pt_stats_survey_like(
        _compute_corrfunc,
        catalogue_data, catalogue_rand,
        los_data=los_data, los_rand=los_rand,
        paramset=paramset, params_sampling=sampling_params,
        degree=degree, binning=binning,
        types={'catalogue_type': 'survey', 'statistic_type': '2pcf'},
        save=save, logger=logger
    )

    # if logger:
    #     logger.info(
    #         "... measured two-point correlation function "
    #         "from paired survey-type catalogues.",
    #         cpp_state='end'
    #     )

    return results


# ========================================================================
# Simulation statistics
# ========================================================================

def _compute_2pt_stats_sim_like(twopt_algofunc, catalogue_data,
                                paramset=None, params_sampling=None,
                                degree=None, binning=None, types=None,
                                save=False, logger=None):
    """Compute two-point statistics from a simulation-box catalogue
    in the global plane-parallel approximation.

    Parameters
    ----------
    twopt_algofunc : callable
        Two-point statistic algorithmic function.
    catalogue_data : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Data-source catalogue.
    paramset : :class:`~triumvirate.parameters.ParameterSet`, optional
        Full parameter set.  If `None` (default), `degree`, `binning`
        and `sampling_params` must be provided.
    params_sampling : dict, optional
        Dictionary containing a subset of the following entries
        for sampling parameters---

        - 'boxalign': {'centre', 'pad'};
        - 'boxsize': sequence of [float, float, float];
        - 'ngrid': sequence of [int, int, int];
        - 'assignment': {'ngp', 'cic', 'tsc', 'pcs'};
        - 'interlace': bool;

        and one and only one of the following when 'boxalign' is 'pad'---

        - 'boxpad': float;
        - 'gridpad': float.

        This will override corresponding entries in `paramset`.
    degree : int, optional
        Multipole degree.  If not `None` (default), this will override
        ``paramset['degrees']['ELL']``.
    binning : :class:`~triumvirate.dataobjs.Binning`, optional
        Binning for the measurements.  If `None` (default), this is
        constructed from `paramset`.
    types : dict, optional
        `catalogue_type` and `statistic_type` (default is `None`).
        This should be set by the caller of this function.
    save : {'.txt', '.npz', False}, optional
        If not `False` (default), save the measurements as a '.txt' file
        or in '.npz' format. The save path is determined from `paramset`
        (if unset, a default file in the current working directory is
        used).
    logger : :class:`logging.Logger`, optional
        Logger (default is `None`).

    Returns
    -------
    results : dict of {str: :class:`numpy.ndarray`}
        Measurement results.

    Raises
    ------
    ValueError
        When `paramset` is `None` while `params_sampling`, `degree` or
        `binning` is also `None`.

    """
    # --------------------------------------------------------------------
    # Initialisation
    # --------------------------------------------------------------------

    # -- Parameters ------------------------------------------------------

    paramset = _amalgamate_parameters(
        paramset=paramset, params_sampling=params_sampling,
        degree=degree, binning=binning, **types
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

    particles_data =  \
        catalogue_data._convert_to_cpp_catalogue(verbose=paramset['verbose'])

    if logger:
        logger.info(
            "... prepared catalogue for clustering algorithm.",
            cpp_state='end'
        )

    # Set up constants.
    if paramset['norm_convention'] == 'particle':
        norm_factor = _calc_powspec_normalisation_from_particles(
            particles_data, alpha=1.
        )
        norm_factor_alt = _calc_powspec_normalisation_from_mesh(
            particles_data, paramset, alpha=1.
        )
    if paramset['norm_convention'] == 'mesh':
        norm_factor = _calc_powspec_normalisation_from_mesh(
            particles_data, paramset, alpha=1.
        )
        norm_factor_alt = _calc_powspec_normalisation_from_particles(
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

    results = twopt_algofunc(particles_data, paramset, binning, norm_factor)

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
                ['%.9e'] * 2 + ['%10d'] + ['% .9e'] * (datatab.shape[-1] - 3)
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


def compute_powspec_in_gpp_box(catalogue_data,
                               degree=None, binning=None, sampling_params=None,
                               paramset=None,
                               save=False, logger=None):
    """Compute power spectrum from a simulation-box catalogue
    in the global plane-parallel approximation.

    Parameters
    ----------
    catalogue_data : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Data-source catalogue.
    degree : int, optional
        Multipole degree.  If not `None` (default), this will override
        ``paramset['degrees']['ELL']``.
    binning : :class:`~triumvirate.dataobjs.Binning`, optional
        Binning for the measurements.  If `None` (default), this is
        constructed from `paramset`.
    sampling_params : dict, optional
        Dictionary containing a subset of the following entries
        for sampling parameters---

        - 'boxalign': {'centre', 'pad'};
        - 'boxsize': sequence of [float, float, float];
        - 'ngrid': sequence of [int, int, int];
        - 'assignment': {'ngp', 'cic', 'tsc', 'pcs'};
        - 'interlace': bool;

        and one and only one of the following when 'boxalign' is 'pad'---

        - 'boxpad': float;
        - 'gridpad': float.

        This will override corresponding entries in `paramset`.
    paramset : :class:`~triumvirate.parameters.ParameterSet`, optional
        Full parameter set (default is `None`).  This is used
        in lieu of `degree`, `binning` or `sampling_params`.
    save : {'.txt', '.npz', False}, optional
        If not `False` (default), save the measurements as a '.txt' file
        or in '.npz' format. The save path is determined from `paramset`
        (if unset, a default file in the current working directory is
        used).
    logger : :class:`logging.Logger`, optional
        Logger (default is `None`).

    Returns
    -------
    results : dict of {str: :class:`numpy.ndarray`}
        Measurement results.

    Raises
    ------
    ValueError
        When `paramset` is `None` but `degree`, `binning` or
        `sampling_params` is also `None`.

    Examples
    --------
    See analogous examples in :func:`~triumvirate.twopt.compute_powspec`
    (though without line-of-sight arguments).

    """
    # if logger:
    #     logger.info(
    #         "Measuring power spectrum from a simulation-box catalogue "
    #         "in the global plane-parallel approximation...",
    #         cpp_state='start'
    #     )

    results = _compute_2pt_stats_sim_like(
        _compute_powspec_in_gpp_box, catalogue_data,
        paramset=paramset, params_sampling=sampling_params,
        degree=degree, binning=binning,
        types={'catalogue_type': 'sim', 'statistic_type': 'powspec'},
        save=save, logger=logger
    )

    # if logger:
    #     logger.info(
    #         "... measured power spectrum from a simulation-box catalogue "
    #         "in the global plane-parallel approximation.",
    #         cpp_state='end'
    #     )

    return results


def compute_corrfunc_in_gpp_box(catalogue_data,
                                degree=None, binning=None,
                                sampling_params=None,
                                paramset=None,
                                save=False, logger=None):
    """Compute correlation function from a simulation-box catalogue
    in the global plane-parallel approximation.

    Parameters
    ----------
    catalogue_data : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Data-source catalogue.
    degree : int, optional
        Multipole degree.  If not `None` (default), this will override
        ``paramset['degrees']['ELL']``.
    binning : :class:`~triumvirate.dataobjs.Binning`, optional
        Binning for the measurements.  If `None` (default), this is
        constructed from `paramset`.
    sampling_params : dict, optional
        Dictionary containing a subset of the following entries
        for sampling parameters---

        - 'boxalign': {'centre', 'pad'};
        - 'boxsize': sequence of [float, float, float];
        - 'ngrid': sequence of [int, int, int];
        - 'assignment': {'ngp', 'cic', 'tsc', 'pcs'};
        - 'interlace': bool;

        and one and only one of the following when 'boxalign' is 'pad'---

        - 'boxpad': float;
        - 'gridpad': float.

        This will override corresponding entries in `paramset`.
    paramset : :class:`~triumvirate.parameters.ParameterSet`, optional
        Full parameter set (default is `None`).  This is used
        in lieu of `degree`, `binning` or `sampling_params`.
    save : {'.txt', '.npz', False}, optional
        If not `False` (default), save the measurements as a '.txt' file
        or in '.npz' format. The save path is determined from `paramset`
        (if unset, a default file in the current working directory is
        used).
    logger : :class:`logging.Logger`, optional
        Logger (default is `None`).

    Returns
    -------
    results : dict of {str: :class:`numpy.ndarray`}
        Measurement results.

    Raises
    ------
    ValueError
        When `paramset` is `None` but `degree`, `binning` or
        `sampling_params` is also `None`.

    Examples
    --------
    See analogous examples in :func:`~triumvirate.twopt.compute_powspec`
    (though without line-of-sight arguments).

    """
    # if logger:
    #     logger.info(
    #         "Measuring two-point correlation function "
    #         "from a simulation-box catalogue "
    #         "in the global plane-parallel approximation...",
    #         cpp_state='start'
    #     )

    results = _compute_2pt_stats_sim_like(
        _compute_corrfunc_in_gpp_box, catalogue_data,
        paramset=paramset, params_sampling=sampling_params,
        degree=degree, binning=binning,
        types={'catalogue_type': 'sim', 'statistic_type': '2pcf'},
        save=save, logger=logger
    )

    # if logger:
    #     logger.info(
    #         "... measured two-point correlation function "
    #         "from a simulation-box catalogue "
    #         "in the global plane-parallel approximation.",
    #         cpp_state='end'
    #     )

    return results


# ========================================================================
# Window statistics
# ========================================================================

def compute_corrfunc_window(catalogue_rand, los_rand=None,
                            degree=None, binning=None, sampling_params=None,
                            paramset=None,
                            save=False, logger=None):
    """Compute correlation function window from a random catalogue.

    Parameters
    ----------
    catalogue_rand : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Random-source catalogue.
    los_rand : (N, 3) array of float, optional
        Specified lines of sight for the random-source catalogue.
        If `None` (default), this is automatically computed using
        :meth:`~triumvirate.catalogue.ParticleCatalogue.compute_los`.
    degree : int, optional
        Multipole degree.  If not `None` (default), this will override
        ``paramset['degrees']['ELL']``.
    binning : :class:`~triumvirate.dataobjs.Binning`, optional
        Binning for the measurements.  If `None` (default), this is
        constructed from `paramset`.
    sampling_params : dict, optional
        Dictionary containing a subset of the following entries
        for sampling parameters---

        - 'boxalign': {'centre', 'pad'};
        - 'boxsize': sequence of [float, float, float];
        - 'ngrid': sequence of [int, int, int];
        - 'assignment': {'ngp', 'cic', 'tsc', 'pcs'};
        - 'interlace': bool;

        and one and only one of the following when 'boxalign' is 'pad'---

        - 'boxpad': float;
        - 'gridpad': float.

        This will override corresponding entries in `paramset`.
    paramset : :class:`~triumvirate.parameters.ParameterSet`, optional
        Full parameter set (default is `None`).  This is used
        in lieu of `degree`, `binning` or `sampling_params`.
    save : {'.txt', '.npz', False}, optional
        If not `False` (default), save the measurements as a '.txt' file
        or in '.npz' format. The save path is determined from `paramset`
        (if unset, a default file in the current working directory is
        used).
    logger : :class:`logging.Logger`, optional
        Logger (default is `None`).

    Returns
    -------
    results : dict of {str: :class:`numpy.ndarray`}
        Measurement results.

    Raises
    ------
    ValueError
        When `paramset` is `None` but `degree`, `binning` or
        `sampling_params` is also `None`.

    Examples
    --------
    See analogous examples in :func:`~triumvirate.twopt.compute_powspec`.

    """
    # --------------------------------------------------------------------
    # Initialisation
    # --------------------------------------------------------------------

    # -- Parameters ------------------------------------------------------

    paramset = _amalgamate_parameters(
        paramset=paramset, params_sampling=sampling_params,
        degree=degree, binning=binning,
        catalogue_type='random', statistic_type='2pcf-win'
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
    particles_rand =  \
        catalogue_rand._convert_to_cpp_catalogue(verbose=paramset['verbose'])

    # Set up constants.
    if paramset['norm_convention'] == 'particle':
        norm_factor = _calc_powspec_normalisation_from_particles(
            particles_rand, alpha=1.
        )
        norm_factor_alt = _calc_powspec_normalisation_from_mesh(
            particles_rand, paramset, alpha=1.
        )
    if paramset['norm_convention'] == 'mesh':
        norm_factor = _calc_powspec_normalisation_from_mesh(
            particles_rand, paramset, alpha=1.
        )
        norm_factor_alt = _calc_powspec_normalisation_from_particles(
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

    results = _compute_corrfunc_window(
        particles_rand, los_rand, paramset, binning,
        alpha=1., norm_factor=norm_factor
    )

    if logger:
        logger.info(
            "... measured window function statistics.", cpp_state='end'
        )

    if save:
        odirpath = paramset['directories']['measurements']
        header = "\n".join([
            catalogue_rand.write_attrs_as_header(),
            _print_measurement_header(paramset, norm_factor, norm_factor_alt),
        ])
        if save.lower() == '.txt':
            datatab = _assemble_measurement_datatab(results, paramset)
            datafmt = '\t'.join(
                ['%.9e'] * 2 + ['%10d'] + ['% .9e'] * (datatab.shape[-1] - 3)
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
