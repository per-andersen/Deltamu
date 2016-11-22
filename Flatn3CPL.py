from astropy import units as u
from n3CPL import n3CPL
import scalar_inv_efuncs

class Flatn3CPL(n3CPL):
    """FLRW cosmology with a n=3 nCDM dark energy equation of state and no curvature.

    The equation for the dark energy equation of state uses the
    nCPL form as described in Pantazis et al. (2016):
    :math:`w(z) = w_0 + w_a (1-a)^n = w_0 + w_a z^n / (1+z)^n`.

    Parameters
    ----------

    H0 : float or `~astropy.units.Quantity`
        Hubble constant at z = 0. If a float, must be in [km/sec/Mpc]

    Om0 : float
        Omega matter: density of non-relativistic matter in units of the
        critical density at z=0.

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
    >>> from astropy.cosmology import Flatw0waCDM
    >>> cosmo = Flatn3CPL(H0=70, Om0=0.3, w0=-0.9, wa=0.2)

    The comoving distance in Mpc at redshift z:

    >>> z = 0.5
    >>> dc = cosmo.comoving_distance(z)
    """

    def __init__(self, H0, Om0, w0=-1., wa=0., Tcmb0=2.725,
                 Neff=3.04, m_nu=u.Quantity(0.0, u.eV), Ob0=None, name=None):

        n3CPL.__init__(self, H0, Om0, 0.0, w0=w0, wa=wa, Tcmb0=Tcmb0,
                         Neff=Neff, m_nu=m_nu, name=name, Ob0=Ob0)
        # Do some twiddling after the fact to get flatness
        self._Ode0 = 1.0 - self._Om0 - self._Ogamma0 - self._Onu0
        self._Ok0 = 0.0

        # Please see "Notes about speeding up integrals" for discussion
        # about what is being done here.
        '''
        if self._Tcmb0.value == 0:
            self._inv_efunc_scalar = scalar_inv_efuncs.fw0wacdm_inv_efunc_norel
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0,
                                           self._w0, self._wa)
        elif not self._massivenu:
            self._inv_efunc_scalar = scalar_inv_efuncs.fw0wacdm_inv_efunc_nomnu
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0,
                                           self._Ogamma0 + self._Onu0,
                                           self._w0, self._wa)
        else:
            self._inv_efunc_scalar = scalar_inv_efuncs.fw0wacdm_inv_efunc
            self._inv_efunc_scalar_args = (self._Om0, self._Ode0,
                                           self._Ogamma0, self._neff_per_nu,
                                           self._nmasslessnu,
                                           self._nu_y_list, self._w0,
                                           self._wa)
        '''
        
        self._inv_efunc_scalar = scalar_inv_efuncs.fn3cdm_inv_efunc
        self._inv_efunc_scalar_args = (self._Om0, self._Ode0,
                                           self._Ogamma0, self._w0,
                                           self._wa)
        
    def __repr__(self):
        retstr = "{0}H0={1:.3g}, Om0={2:.3g}, "\
                 "w0={3:.3g}, Tcmb0={4:.4g}, Neff={5:.3g}, m_nu={6}, "\
                 "Ob0={7:s})"
        return retstr.format(self._namelead(), self._H0, self._Om0, self._w0,
                             self._Tcmb0, self._Neff, self.m_nu,
                             _float_or_none(self._Ob0))
