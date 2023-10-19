"""
Hankel-Like Transforms (:mod:`~triumvirate.transforms`)
==========================================================================

Configure the program logger.

.. autosummary::
    SphericalBesselTransform

"""
import numpy as np

from ._fftlog import HankelTransform


class SphericalBesselTransform:
    """Spherical Bessel Transform.

    Parameters
    ----------
    degree : int
        Degree of the spherical Bessel transform.
    bias : int
        Power-law bias index.
    sample_pts : array of float
        Pre-transform sample points in 1-d.  Must be log-linearly spaced.
        Must be even in length if extrapolation is used.
    pivot : float, optional
        Pivot value (default is 1.).  When `lowring` is `True`, this is
        adjusted if it is non-zero, or otherwise directly calculated.
    lowring : bool, optional
        Low-ringing condition (default is `True`).
    extrap : int, optional
        Extrapolation method (default is 0) with the following
        options:

        - 0: none;
        - 1: extrapolate by constant padding;
        - 2: extrapolate linearly;
        - 3: extrapolate log-linearly.

        Any extrapolation doubles the sample size in effect, and
        assumes the pre-transform samples to be real.

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

    def __init__(self, degree, bias, sample_pts, pivot=1.,
                 lowring=True, extrap=0):

        self.degree = degree
        self.bias = bias

        self._fbht = HankelTransform(
            degree + 1./2, bias, sample_pts,
            kr_c=pivot, lowring=lowring, extrap=extrap
        )
        self._lowring = lowring
        self._extrap = extrap

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
            Pre-transform samples in 1-d.  Assumed to be real
            if extrapolation is used.

        Returns
        -------
        post_sampts, post_samples : array_like
            Post-transform sample points and samples.

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
            Pre-transform samples in 1-d.  Assumed to be real
            if extrapolation is used.

        Returns
        -------
        post_sampts, post_samples : array of float
            Post-transform sample points and samples.

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


class DoubleSphericalBesselTransform:
    """Double spherical Bessel Transform.

    Parameters
    ----------
    degrees : (int, int)
        Degrees of the double spherical Bessel transform.
    biases : (int, int)
        Power-law bias indices.
    sample_pts : array of float
        Pre-transform sample points in 1-d for both dimensions.
        Must be log-linearly spaced.  Must be even in length if
        extrapolation is used.
    pivot : float, optional
        Pivot value (default is 1.).  When `lowring` is `True`, this is
        adjusted if it is non-zero, or otherwise directly calculated.
    lowring : bool, optional
        Low-ringing condition (default is `True`).
    extrap : int, optional
        Extrapolation method (default is 0) with the following
        options:

        - 0: none;
        - 1: extrapolate by constant padding;
        - 2: extrapolate linearly;
        - 3: extrapolate log-linearly.

        Any extrapolation doubles the sample size in effect, and
        assumes the pre-transform samples to be real.

    Attributes
    ----------
    degrees : (int, int)
        Degrees of the spherical Bessel transform.
    biases : (int, int)
        Power-law bias indices.
    size : int
        Sample size for either dimension.
    pivot : float
        Pivot value used for both dimensions.

    """
    # Representative transform, either vertical (0, across rows)
    # or horizontal (1, across columns).
    _rep = 0

    def __init__(self, degrees, biases, sample_pts, pivot=1.,
                 lowring=True, extrap=0):

        self.degrees = degrees
        self.biases = biases

        self._fbsjt = tuple(
            SphericalBesselTransform(
                degree_, bias_, sample_pts,
                pivot=pivot, lowring=lowring, extrap=extrap
            )
            for (degree_, bias_) in zip(degrees, biases)
        )
        self._lowring = lowring
        self._extrap = extrap

        self._logres = self._fbsjt[self._rep]._logres
        self._pre_sampts = np.asarray(self._fbsjt[self._rep]._pre_sampts)
        self._post_sampts = np.asarray(self._fbsjt[self._rep]._post_sampts)

    @property
    def size(self):
        """Sample size.

        """
        return self._fbsjt[self._rep]._nsamp

    @property
    def pivot(self):
        """Pivot value.

        """
        return self._fbsjt[self._rep]._pivot

    def transform(self, pre_samples):
        """Transform samples at initialised sample points.

        Parameters
        ----------
        pre_samples : array_like
            Pre-transform samples in 2-d.  Assumed to be real
            if extrapolation is used.

        Returns
        -------
        post_sampts, post_samples : array_like
            Post-transform sample points and samples in 2-d, with
            `post_sampts` in :func:`numpy.meshgrid` 'ij'-indexing format.

        """
        inter_samples = []
        for pre_samples_row in pre_samples:
            _, inter_samples_row = self._fbsjt[1].transform(pre_samples_row)
            inter_samples.append(inter_samples_row)

        post_samples_T = []
        for inter_samples_col in np.asarray(inter_samples).T:
            _, post_samples_col = self._fbsjt[0].transform(inter_samples_col)
            post_samples_T.append(post_samples_col)

        post_samples = np.asarray(post_samples_T).T

        post_sampts = np.meshgrid(
            self._post_sampts, self._post_sampts, indexing='ij'
        )

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
            Pre-transform samples in 2-d.  Assumed to be real
            if extrapolation is used.

        Returns
        -------
        post_sampts, post_samples : array of float
            Post-transform sample points and samples in 2-d, with
            `post_sampts` in :func:`numpy.meshgrid` 'ij'-indexing format.

        """
        if direction in {'forward', 1}:
            parity = (-1j) ** sum(self.degrees)
            prefactor = 1.
        elif direction in {'backward', -1}:
            parity = (1j) ** sum(self.degrees)
            prefactor = (2*np.pi) ** (-3 * len(self.degrees))
        else:
            raise ValueError(
                "Transform direction must be 'forward'/1 or 'backward'/-1."
            )

        post_sampts, post_samples = self.transform(pre_samples)

        post_samples *= prefactor * parity

        return post_sampts, post_samples
