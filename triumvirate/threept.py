"""
Three-Point Correlator Measurements (:mod:`~triumvirate.threept`)
==========================================================================

Measuring three-point correlator statistics from catalogues.

"""
import numpy as np

from catalogue import _prepare_catalogue
from parameters import InvalidParameter
from _threept import (
    _calc_bispec_normalisation_from_mesh,
    _calc_bispec_normalisation_from_particles,
    _compute_bispec,
    # _compute_bispec_for_los_choice,
    _compute_bispec_in_box,
    _compute_3pcf,
    _compute_3pcf_in_box,
    _compute_3pcf_window,
)


def compute_bispec(catalogue_data, catalogue_rand, params,
                   los_data=None, los_rand=None,
                   box_align='centre', ngrid_pad=None,
                   save=False, logger=None):
    """Compute bispectrum from data and random catalogues.

    Parameters
    ----------
    catalogue_data : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Data-source catalogue.
    catalogue_rand : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Random-source catalogue.
    params : :class:`~triumvirate.parameters.ParameterSet`
        Measurement parameters.
    los_data : (N, 3) array of float, optional
        Specified lines of sight for the data-source catalogue.
        If `None` (default), this is automatically computed using
        :meth:`~triumvirate.catalogue.ParticleCatalogue.compute_los`.
    los_rand : (N, 3) array of float, optional
        Specified lines of sight for the random-source catalogue.
        If `None` (default), this is automatically computed using
        :meth:`~triumvirate.catalogue.ParticleCatalogue.compute_los`.
    box_align : {'centre', 'pad'}, optional
        Alignment of the catalogue(s) inside the box.  The random-source
        catalogue is used as the reference to 'centre' (default) or
        'pad' (away from the box origin corner) both catalogues.
    ngrid_pad : (array of) float, optional
        Number of grids as padding for the catalogues away from the
        origin corner of the box.  If `None` (default), the default
        padding factor is assumed (see
        :meth:`triumvirate.catalogue.ParticleCatalogue.pad`).
    save : bool, optional
        If `True` (default is `False`), measurement results are
        automatically saved to an output file specified from `params`.
    logger : :class:`logging.Logger`, optional
        Logger (default is `None`).

    Returns
    -------
    results : dict of {str: :class:`numpy.ndarray`}
        Measurement results.

    """
    # Prepare catalogues.
    if los_data is None:
        los_data = catalogue_data.compute_los()
    if los_rand is None:
        los_rand = catalogue_rand.compute_los()

    los_data = np.ascontiguousarray(los_data)
    los_rand = np.ascontiguousarray(los_rand)

    if box_align.lower() == 'centre':
        catalogue_data.centre(
            [params['boxsize'][axis] for axis in ['x', 'y', 'z']],
            catalogue_ref=catalogue_rand
        )
    elif box_align.lower() == 'pad':
        kwargs = {} if ngrid_pad is None else {'ngrid_pad': ngrid_pad}
        catalogue_data.pad(
            [params['boxsize'][axis] for axis in ['x', 'y', 'z']],
            [params['ngrid'][axis] for axis in ['x', 'y', 'z']],
            catalogue_ref=catalogue_rand,
            **kwargs
        )
    else:
        raise ValueError("`box_alignment` must be 'centre' or 'pad'.")

    particles_data = _prepare_catalogue(catalogue_data)
    particles_rand = _prepare_catalogue(catalogue_rand)

    # Compute auxiliary quantities.
    kbin = np.ascontiguousarray(
        np.linspace(*params['range'], num=params['dim'])
    )

    alpha = catalogue_data.wtotal / catalogue_rand.wtotal

    try:
        logger.info("Calculating normalisation...", cpp_state='start')
    except (AttributeError, TypeError):
        pass

    if params['norm_convention'] == 'mesh':
        norm = _calc_bispec_normalisation_from_mesh(
            particles_rand, params, alpha
        )
    elif params['norm_convention'] == 'particle':
        norm = _calc_bispec_normalisation_from_particles(particles_rand, alpha)
    else:
        raise InvalidParameter("Invalid `norm_convention` parameter.")

    try:
        logger.info("... calculated normalisation.", cpp_state='end')
    except (AttributeError, TypeError):
        pass

    if logger:
        logger.info("Alpha contrast: %.6e.", alpha)
        logger.info("Normalisation constant: %.6e.", norm)

    # Perform measurement.
    try:
        logger.info("Making measurements...", cpp_state='start')
    except (AttributeError, TypeError):
        pass

    results = _compute_bispec(
        particles_data, particles_rand, los_data, los_rand,
        params, kbin, alpha, norm, save=save
    )

    try:
        logger.info("... made measurements.", cpp_state='end')
    except (AttributeError, TypeError):
        pass

    return results

def compute_3pcf(catalogue_data, catalogue_rand, params,
                 los_data=None, los_rand=None,
                 box_align='centre', ngrid_pad=None,
                 save=False, logger=None):
    """Compute three-point correlation function from data and
    random catalogues.

    Parameters
    ----------
    catalogue_data : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Data-source catalogue.
    catalogue_rand : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Random-source catalogue.
    params : :class:`~triumvirate.parameters.ParameterSet`
        Measurement parameters.
    los_data : (N, 3) array of float, optional
        Specified lines of sight for the data-source catalogue.
        If `None` (default), this is automatically computed using
        :meth:`~triumvirate.catalogue.ParticleCatalogue.compute_los`.
    los_rand : (N, 3) array of float, optional
        Specified lines of sight for the random-source catalogue.
        If `None` (default), this is automatically computed using
        :meth:`~triumvirate.catalogue.ParticleCatalogue.compute_los`.
    box_align : {'centre', 'pad'}, optional
        Alignment of the catalogue(s) inside the box.  The random-source
        catalogue is used as the reference to 'centre' (default) or
        'pad' (away from the box origin corner) both catalogues.
    ngrid_pad : (array of) float, optional
        Number of grids as padding for the catalogues away from the
        origin corner of the box.  If `None` (default), the default
        padding factor is assumed (see
        :meth:`triumvirate.catalogue.ParticleCatalogue.pad`).
    save : bool, optional
        If `True` (default is `False`), measurement results are
        automatically saved to an output file specified from `params`.
    logger : :class:`logging.Logger`, optional
        Logger (default is `None`).

    Returns
    -------
    results : dict of {str: :class:`numpy.ndarray`}
        Measurement results.

    """
    # Prepare catalogues.
    if los_data is None:
        los_data = catalogue_data.compute_los()
    if los_rand is None:
        los_rand = catalogue_rand.compute_los()

    los_data = np.ascontiguousarray(los_data)
    los_rand = np.ascontiguousarray(los_rand)

    if box_align.lower() == 'centre':
        catalogue_data.centre(
            [params['boxsize'][axis] for axis in ['x', 'y', 'z']],
            catalogue_ref=catalogue_rand
        )
    elif box_align.lower() == 'pad':
        kwargs = {} if ngrid_pad is None else {'ngrid_pad': ngrid_pad}
        catalogue_data.pad(
            [params['boxsize'][axis] for axis in ['x', 'y', 'z']],
            [params['ngrid'][axis] for axis in ['x', 'y', 'z']],
            catalogue_ref=catalogue_rand,
            **kwargs
        )
    else:
        raise ValueError("`box_alignment` must be 'centre' or 'pad'.")

    particles_data = _prepare_catalogue(catalogue_data)
    particles_rand = _prepare_catalogue(catalogue_rand)

    # Compute auxiliary quantities.
    rbin = np.ascontiguousarray(
        np.linspace(*params['range'], num=params['dim'])
    )

    alpha = catalogue_data.wtotal / catalogue_rand.wtotal

    try:
        logger.info("Calculating normalisation...", cpp_state='start')
    except (AttributeError, TypeError):
        pass

    if params['norm_convention'] == 'mesh':
        norm = _calc_bispec_normalisation_from_mesh(
            particles_rand, params, alpha
        )
    elif params['norm_convention'] == 'particle':
        norm = _calc_bispec_normalisation_from_particles(particles_rand, alpha)
    else:
        raise InvalidParameter("Invalid `norm_convention` parameter.")

    try:
        logger.info("... calculated normalisation.", cpp_state='end')
    except (AttributeError, TypeError):
        pass

    if logger:
        logger.info("Alpha contrast: %.6e.", alpha)
        logger.info("Normalisation constant: %.6e.", norm)

    # Perform measurement.
    try:
        logger.info("Making measurements...", cpp_state='start')
    except (AttributeError, TypeError):
        pass

    results = _compute_3pcf(
        particles_data, particles_rand, los_data, los_rand,
        params, rbin, alpha, norm, save=save
    )

    try:
        logger.info("... made measurements.", cpp_state='end')
    except (AttributeError, TypeError):
        pass

    return results

def compute_bispec_in_box(catalogue_data, params, save=False, logger=None):
    """Compute power spectrum in a box.

    Parameters
    ----------
    catalogue_data : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Data-source catalogue.
    params : :class:`~triumvirate.parameters.ParameterSet`
        Measurement parameters.
    save : bool, optional
        If `True` (default is `False`), measurement results are
        automatically saved to an output file specified from `params`.
    logger : :class:`logging.Logger`, optional
        Logger (default is `None`).

    Returns
    -------
    results : dict of {str: :class:`numpy.ndarray`}
        Measurement results.

    """
    # Prepare catalogues.
    catalogue_data.periodise(
        [params['boxsize'][axis] for axis in ['x', 'y', 'z']]
    )

    particles_data = _prepare_catalogue(catalogue_data)

    # Compute auxiliary quantities.
    kbin = np.ascontiguousarray(
        np.linspace(*params['range'], num=params['dim'])
    )

    try:
        logger.info("Calculating normalisation...", cpp_state='start')
    except (AttributeError, TypeError):
        pass

    if params['norm_convention'] == 'mesh':
        norm = _calc_bispec_normalisation_from_mesh(
            particles_data, params, alpha=1.)
    elif params['norm_convention'] == 'particle':
        norm = _calc_bispec_normalisation_from_particles(
            particles_data, alpha=1.
        )
    else:
        raise InvalidParameter("Invalid `norm_convention` parameter.")

    try:
        logger.info("... calculated normalisation.", cpp_state='end')
    except (AttributeError, TypeError):
        pass

    if logger:
        logger.info("Normalisation constant: %.6e.", norm)

    # Perform measurement.
    try:
        logger.info("Making measurements...", cpp_state='start')
    except (AttributeError, TypeError):
        pass

    results = _compute_bispec_in_box(
        particles_data, params, kbin, norm, save=save
    )

    try:
        logger.info("... made measurements.", cpp_state='end')
    except (AttributeError, TypeError):
        pass

    return results

def compute_3pcf_in_box(catalogue_data, params, save=False, logger=None):
    """Compute correlation function in a box.

    Parameters
    ----------
    catalogue_data : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Data-source catalogue.
    params : :class:`~triumvirate.parameters.ParameterSet`
        Measurement parameters.
    save : bool, optional
        If `True` (default is `False`), measurement results are
        automatically saved to an output file specified from `params`.
    logger : :class:`logging.Logger`, optional
        Logger (default is `None`).

    Returns
    -------
    results : dict of {str: :class:`numpy.ndarray`}
        Measurement results.

    """
    # Prepare catalogues.
    catalogue_data.periodise(
        [params['boxsize'][axis] for axis in ['x', 'y', 'z']]
    )

    particles_data = _prepare_catalogue(catalogue_data)

    # Compute auxiliary quantities.
    rbin = np.ascontiguousarray(
        np.linspace(*params['range'], num=params['dim'])
    )

    try:
        logger.info("Calculating normalisation...", cpp_state='start')
    except (AttributeError, TypeError):
        pass

    if params['norm_convention'] == 'mesh':
        norm = _calc_bispec_normalisation_from_mesh(
            particles_data, params, alpha=1.)
    elif params['norm_convention'] == 'particle':
        norm = _calc_bispec_normalisation_from_particles(
            particles_data, alpha=1.
        )
    else:
        raise InvalidParameter("Invalid `norm_convention` parameter.")

    try:
        logger.info("... calculated normalisation.", cpp_state='end')
    except (AttributeError, TypeError):
        pass

    if logger:
        logger.info("Normalisation constant: %.6e.", norm)

    # Perform measurement.
    try:
        logger.info("Making measurements...", cpp_state='start')
    except (AttributeError, TypeError):
        pass

    results = _compute_3pcf_in_box(
        particles_data, params, rbin, norm, save=save
    )

    try:
        logger.info("... made measurements.", cpp_state='end')
    except (AttributeError, TypeError):
        pass

    return results

def compute_3pcf_window(catalogue_rand, params, los_rand=None, wide_angle=False,
                        save=False, logger=None):
    """Compute correlation function window from a random catalogue.

    Parameters
    ----------
    catalogue_rand : :class:`~triumvirate.catalogue.ParticleCatalogue`
        Random-source catalogue.
    params : :class:`~triumvirate.parameters.ParameterSet`
        Measurement parameters.
    los_rand : (N, 3) array of float, optional
        Specified lines of sight for the random-source catalogue.
        If `None` (default), this is automatically computed using
        :meth:`~triumvirate.catalogue.ParticleCatalogue.compute_los`.
    wide_angle : bool, optional
        If `True` (default is `False`), wide-angle correction terms
        are computed.
    save : bool, optional
        If `True` (default is `False`), measurement results are
        automatically saved to an output file specified from `params`.
    logger : :class:`logging.Logger`, optional
        Logger (default is `None`).

    Returns
    -------
    results : dict of {str: :class:`numpy.ndarray`}
        Measurement results.

    """
    # Prepare catalogues.
    if los_rand is None:
        los_rand = catalogue_rand.compute_los()
    los_rand = np.ascontiguousarray(los_rand)

    catalogue_rand.centre(
        [params['boxsize'][axis] for axis in ['x', 'y', 'z']]
    )

    particles_rand = _prepare_catalogue(catalogue_rand)

    # Compute auxiliary quantities.
    rbin = np.ascontiguousarray(
        np.linspace(*params['range'], num=params['dim'])
    )

    try:
        logger.info("Calculating normalisation...", cpp_state='start')
    except (AttributeError, TypeError):
        pass

    if params['norm_convention'] == 'mesh':
        norm = _calc_bispec_normalisation_from_mesh(
            particles_rand, params, alpha=1.)
    elif params['norm_convention'] == 'particle':
        norm = _calc_bispec_normalisation_from_particles(
            particles_rand, alpha=1.
        )
    else:
        raise InvalidParameter("Invalid `norm_convention` parameter.")

    try:
        logger.info("... calculated normalisation.", cpp_state='end')
    except (AttributeError, TypeError):
        pass

    if logger:
        logger.info("Normalisation constant: %.6e.", norm)

    # Perform measurement.
    try:
        logger.info("Making measurements...", cpp_state='start')
    except (AttributeError, TypeError):
        pass

    results = _compute_3pcf_window(
        particles_rand, los_rand,
        params, rbin, alpha=1., norm=norm, wide_angle=wide_angle, save=save
    )

    try:
        logger.info("... made measurements.", cpp_state='end')
    except (AttributeError, TypeError):
        pass

    return results

# def compute_bispec_for_los_choice(catalogue_data, catalogue_rand, params,
#                    los_choice, los_data=None, los_rand=None,
#                    box_align='centre', ngrid_pad=None,
#                    save=False, logger=None):
#     """Compute bispectrum from data and random catalogues.
#
#     Parameters
#     ----------
#     catalogue_data : :class:`~triumvirate.catalogue.ParticleCatalogue`
#         Data-source catalogue.
#     catalogue_rand : :class:`~triumvirate.catalogue.ParticleCatalogue`
#         Random-source catalogue.
#     params : :class:`~triumvirate.parameters.ParameterSet`
#         Measurement parameters.
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
#     box_align : {'centre', 'pad'}, optional
#         Alignment of the catalogue(s) inside the box.  The random-source
#         catalogue is used as the reference to 'centre' (default) or
#         'pad' (away from the box origin corner) both catalogues.
#     ngrid_pad : (array of) float, optional
#         Number of grids as padding for the catalogues away from the
#         origin corner of the box.  If `None` (default), the default
#         padding factor is assumed (see
#         :meth:`triumvirate.catalogue.ParticleCatalogue.pad`).
#     save : bool, optional
#         If `True` (default is `False`), measurement results are
#         automatically saved to an output file specified from `params`.
#     logger : :class:`logging.Logger`, optional
#         Logger (default is `None`).
#
#     Returns
#     -------
#     results : dict of {str: :class:`numpy.ndarray`}
#         Measurement results.
#
#     """
#     # Prepare catalogues.
#     if los_data is None:
#         los_data = catalogue_data.compute_los()
#     if los_rand is None:
#         los_rand = catalogue_rand.compute_los()
#
#     los_data = np.ascontiguousarray(los_data)
#     los_rand = np.ascontiguousarray(los_rand)
#
#     if box_align.lower() == 'centre':
#         catalogue_data.centre(
#             [params['boxsize'][axis] for axis in ['x', 'y', 'z']],
#             catalogue_ref=catalogue_rand
#         )
#     elif box_align.lower() == 'pad':
#         kwargs = {} if ngrid_pad is None else {'ngrid_pad': ngrid_pad}
#         catalogue_data.pad(
#             [params['boxsize'][axis] for axis in ['x', 'y', 'z']],
#             [params['ngrid'][axis] for axis in ['x', 'y', 'z']],
#             catalogue_ref=catalogue_rand,
#             **kwargs
#         )
#     else:
#         raise ValueError("`box_alignment` must be 'centre' or 'pad'.")
#
#     particles_data = _prepare_catalogue(catalogue_data)
#     particles_rand = _prepare_catalogue(catalogue_rand)
#
#     # Compute auxiliary quantities.
#     kbin = np.ascontiguousarray(
#         np.linspace(*params['range'], num=params['dim'])
#     )
#
#     alpha = catalogue_data.wtotal / catalogue_rand.wtotal
#
#     try:
#         logger.info("Calculating normalisation...", cpp_state='start')
#     except (AttributeError, TypeError):
#         pass
#
#     if params['norm_convention'] == 'mesh':
#         norm = _calc_bispec_normalisation_from_mesh(
#             particles_rand, params, alpha
#         )
#     elif params['norm_convention'] == 'particle':
#         norm = _calc_bispec_normalisation_from_particles(
#             particles_rand, alpha
#         )
#     else:
#         raise InvalidParameter("Invalid `norm_convention` parameter.")
#
#     try:
#         logger.info("... calculated normalisation.", cpp_state='end')
#     except (AttributeError, TypeError):
#         pass
#
#     if logger:
#         logger.info("Alpha contrast: %.6e.", alpha)
#         logger.info("Normalisation constant: %.6e.", norm)
#
#     # Perform measurement.
#     try:
#         logger.info("Making measurements...", cpp_state='start')
#     except (AttributeError, TypeError):
#         pass
#
#     results = _compute_bispec_for_los_choice(
#         particles_data, particles_rand, los_data, los_rand,
#         params, kbin, alpha, norm, los_choice, save=save
#     )
#
#     try:
#         logger.info("... made measurements.", cpp_state='end')
#     except (AttributeError, TypeError):
#         pass
#
#     return results
