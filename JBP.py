from astropy.cosmology import FLRW
from astropy import units as u
from astropy.utils.misc import isiterable
import numpy as np

class JBP(FLRW):
    """FLRW cosmology with a JBP dark energy equation of state and curvature.

    The equation for the dark energy equation of state uses the
    JBP form as described in Jassal et al. (2004):
    :math:`w(z) = w_0 + w_a (a-a^2) = w_0 + w_a z / (1+z)^2`.

    Parameters
    ----------
    H0 : float or `~astropy.units.Quantity`
        Hubble constant at z = 0. If a float, must be in [km/sec/Mpc]

    Om0 : float
        Omega matter: density of non-relativistic matter in units of the
        critical density at z=0.

    Ode0 : float
        Omega dark energy: density of dark energy in units of the critical
        density at z=0.

    w0 : float, optional
        Dark energy equation of state at z=0 (a=1). This is pressure/density
        for dark energy in units where c=1.

    wa : float, optional
        Negative derivative of the dark energy equation of state with respect
        to the scale factor. A cosmological constant has w0=-1.0 and wa=0.0.

    Tcmb0 : float or scalar `~astropy.units.Quantity`, optional
        Temperature of the CMB z=0. If a float, must be in [K].
        Default: 2.725 [K]. Setting this to zero will turn off both photons
        and neutrinos (even massive ones).

    Neff : float, optional
        Effective number of Neutrino species. Default 3.04.

    m_nu : `~astropy.units.Quantity`, optional
        Mass of each neutrino species. If this is a scalar Quantity, then all
        neutrino species are assumed to have that mass. Otherwise, the mass of
        each species. The actual number of neutrino species (and hence the
        number of elements of m_nu if it is not scalar) must be the floor of
        Neff. Typically this means you should provide three neutrino masses
        unless you are considering something like a sterile neutrino.

    Ob0 : float or None, optional
        Omega baryons: density of baryonic matter in units of the critical
        density at z=0.  If this is set to None (the default), any
        computation that requires its value will raise an exception.

    name : str, optional
        Name for this cosmological object.

    Examples
    --------
    >>> from astropy.cosmology import w0waCDM
    >>> cosmo = JBP_CDMCDM(H0=70, Om0=0.3, Ode0=0.7, w0=-0.9, wa=0.2)

    The comoving distance in Mpc at redshift z:

    >>> z = 0.5
    >>> dc = cosmo.comoving_distance(z)
    """

    def __init__(self, H0, Om0, Ode0, w0=-1., wa=0., Tcmb0=2.725,
                 Neff=3.04, m_nu=u.Quantity(0.0, u.eV), Ob0=None, name=None):

        FLRW.__init__(self, H0, Om0, Ode0, Tcmb0, Neff, m_nu, name=name,
                      Ob0=Ob0)
        self._w0 = float(w0)
        self._wa = float(wa)

        # Please see "Notes about speeding up integrals" for discussion
        # about what is being done here.
        '''
        if self._Tcmb0.value == 0:
            self._inv_efunc_scalar = scalar_inv_efuncs.w0wacdm_inv_efunc_norel
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0, self._Ok0,
                                           self._w0, self._wa)
        elif not self._massivenu:
            self._inv_efunc_scalar = scalar_inv_efuncs.w0wacdm_inv_efunc_nomnu
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0, self._Ok0,
                                           self._Ogamma0 + self._Onu0,
                                           self._w0, self._wa)
        else:
            self._inv_efunc_scalar = scalar_inv_efuncs.w0wacdm_inv_efunc
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0, self._Ok0,
                                           self._Ogamma0, self._neff_per_nu,
                                           self._nmasslessnu,
                                           self._nu_y_list, self._w0,
                                           self._wa)
        '''
        

    @property
    def w0(self):
        """ Dark energy equation of state at z=0"""
        return self._w0

    @property
    def wa(self):
        """ Negative derivative of dark energy equation of state w.r.t. a"""
        return self._wa

    def w(self, z):
        """Returns dark energy equation of state at redshift ``z``.

        Parameters
        ----------
        z : array-like
          Input redshifts.

        Returns
        -------
        w : ndarray, or float if input scalar
          The dark energy equation of state

        Notes
        ------
        The dark energy equation of state is defined as
        :math:`w(z) = P(z)/\\rho(z)`, where :math:`P(z)` is the
        pressure at redshift z and :math:`\\rho(z)` is the density
        at redshift z, both in units where c=1.  Here this is
        :math:`w(z) = w_0 + w_a (a - a^2) = w_0 + w_a \\frac{z}{(1+z)^2}`.
        """

        if isiterable(z):
            z = np.asarray(z)

        return self._w0 + self._wa * z / ( (1.0 + z)**2)

    def de_density_scale(self, z):
        """ Evaluates the redshift dependence of the dark energy density.

        Parameters
        ----------
        z : array-like
          Input redshifts.

        Returns
        -------
        I : ndarray, or float if input scalar
          The scaling of the energy density of dark energy with redshift.

        Notes
        -----
        The scaling factor, I, is defined by :math:`\\rho(z) = \\rho_0 I`,
        and in this case is given by

        .. math::

          I = \\left(1 + z\\right)^{3 \\left(1 + w_0 + w_a\\right)}
          \exp \\left(-3 w_a \\frac{z}{1+z}\\right)

        """
        if isiterable(z):
            z = np.asarray(z)
        zp1 = 1.0 + z
        return zp1 ** (3 * (1 + self._w0)) * \
            np.exp(1.5 * self._wa * z**2 / zp1**2)

    def __repr__(self):
        retstr = "{0}H0={1:.3g}, Om0={2:.3g}, "\
                 "Ode0={3:.3g}, w0={4:.3g}, wa={5:.3g}, Tcmb0={6:.4g}, "\
                 "Neff={7:.3g}, m_nu={8}, Ob0={9:s})"
        return retstr.format(self._namelead(), self._H0, self._Om0,
                             self._Ode0, self._w0, self._wa,
                             self._Tcmb0, self._Neff, self.m_nu,
                             _float_or_none(self._Ob0))