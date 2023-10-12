"""
Hankel-Like Transforms (:mod:`~triumvirate.transforms`)
==========================================================================

Configure the program logger.

.. autosummary::
    SphericalBesselTransform

"""
import numpy as np

from triumvirate._fftlog import HankelTransform


class SphericalBesselTransform:
    """Spherical Bessel Transform.

    Parameters
    ----------
    degree : int
        Degree of the spherical Bessel transform.
    bias : int
        Power-law bias index.
    sample_pts : array of float
        Pre-transform sample points.
    pivot : float, optional
        Pivot value (default is 1.).  When `lowring` is `True`, this is
        adjusted if it is non-zero, or otherwise directly calculated.
    lowring : bool, optional
        Low-ringing condition (default is `True`).

    Attributes
    ----------
    degree : int
        Degree of the spherical Bessel transform.
    bias : int
        Power-law bias index.
    size : int
        Sample size.
    pivot : float
        Pivot value used.

    """

    def __init__(self, degree, bias, sample_pts, pivot=1., lowring=True):

        self.degree = degree
        self.bias = bias

        self._fbht = HankelTransform(
            degree + 1./2, bias, sample_pts,
            kr_c=pivot, lowring=lowring
        )

        self._logres = self._fbht._logres
        self._pre_sampts = np.asarray(self._fbht._pre_sampts)
        self._post_sampts = np.asarray(self._fbht._post_sampts)

    @property
    def size(self):
        """Sample size.

        """
        return self._fbht._nsamp

    @property
    def pivot(self):
        """Pivot value.

        """
        return self._fbht._pivot

    def transform(self, pre_samples):
        """Transform samples at initialised sample points.

        Parameters
        ----------
        pre_samples : array_like
            Pre-transform samples.

        Returns
        -------
        post_sampts, post_samples : array_like
            Post-transform samples and sample values.

        """
        pre_samples *= self._pre_sampts ** (3./2)

        post_sampts, post_samples = self._fbht.transform(pre_samples)

        post_samples *= (2*np.pi / self._post_sampts) ** (3./2)

        return post_sampts, post_samples

    def transform_cosmo_multipoles(self, direction, pre_samples):
        """Transform cosmological multipoles between configuration and
        Fourier spaces.

        Parameters
        ----------
        direction : {'forward', 'backward', 1, -1}
            Transform direction: 'forward' or 1 from configuration to
            Fourier space, 'backward' or -1 from Fourier to configuration
            space.
        pre_samples : array of float
            Pre-transform samples.

        Returns
        -------
        post_sampts, post_samples : array of float
            Post-transform samples and sample values.

        """
        if direction in {'forward', 1}:
            parity = (-1j) ** self.degree
            prefactor = 1.
        elif direction in {'backward', -1}:
            parity = (1j) ** self.degree
            prefactor = (2*np.pi) ** (-3)
        else:
            raise ValueError(
                "Transform direction must be 'forward'/1 or 'backward'/-1."
            )

        post_sampts, post_samples = self.transform(pre_samples)

        post_samples *= prefactor * parity

        return post_sampts, post_samples
