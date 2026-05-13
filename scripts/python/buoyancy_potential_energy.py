import gsw
import numpy as np
from scipy import interpolate

def buoyancy_potential_energy(p, SA, CT, Peq, interp=True):
    """
    Buoyancy Potential Energy from accurate energetic formulations

    Parameters
    ----------
    p: Sea pressure (dbar)
    SA: Absolute Salinity (g kg-1)
    CT: Conservative Temperature (°C)
    Peq: Pressure of the equilibrium position of the parcel (dbar)
    interp: If True, Peq is inserted into the pressure vector and SA, CT are
            linearly interpolated to ensure exact thermodynamic consistency (default True)

    Returns
    -------
    BPE_mass: BPE per unit mass (J kg-1)
    BPE_volume: BPE per unit volume (J m-3)
    p: Pressure
    
    Notes
    -----
    TEOS-10 implementation for BPE follows https://doi.org/10.5194/egusphere-2025-3359-RC1
    """
    p = np.asarray(p)
    SA = np.asarray(SA)
    CT = np.asarray(CT)

	# If there is no data at Peq, it is interpolated
    if interp and not np.any(np.isclose(p, Peq)):
        p_ext = np.sort(np.append(p, Peq))
        SA = np.interp(p_ext, p, SA)
        CT = np.interp(p_ext, p, CT)
        p = p_ext

    idx = np.where(np.isclose(p, Peq))[0]
    if idx.size == 0:
        izref = np.argmin(np.abs(p - Peq))
    else:
        izref = idx[0]

    SA_eq = SA[izref]
    CT_eq = CT[izref]
    Peq = p[izref]

    # Computation of BPE
    e1 = gsw.enthalpy(SA_eq, CT_eq, p)
    e2 = gsw.enthalpy(SA_eq, CT_eq, Peq)
    e3 = gsw.enthalpy(35.16504, 0, p)
    e4 = gsw.enthalpy(35.16504, 0, Peq)

    gsdh = gsw.geo_strf_dyn_height(SA, CT, p, Peq)
    gsdh -= gsdh[izref]
    specvol = gsw.specvol(SA_eq, CT_eq, p)

    BPE_mass = e1 - e2 + gsdh - e3 + e4
    BPE_volume = BPE_mass / specvol

    return BPE_mass, BPE_volume, p

def buoyancy_potential_energy_approx(rho, z, zint):
	"""
    Buoyancy Potential Energy from approximate energetic formulations
    
    Parameters
    ----------
    Depth (z) and density (rho) data should go from the deepest to the shallowest and be column vectors
	rho: Sigma-0 Potential Density Anomaly profile (kg m-3)
	z: Heights (m) with negative units. The z-vector can be non-equidistant.
	zint: Reference height (m)

    Returns
    -------
	BPE: Buoyancy potential energy profile (J m-3)
	z: Heights (m)
	"""

	# If there is no data at zint, it is interpolated
	if np.nansum(z[z == zint]) == 0:
		f = interpolate.interp1d(z, rho)
		it = np.argwhere(z >= zint)[0]
		z = np.insert(z, it, zint)
		rho = np.insert(rho, it, f(zint))

	g = 9.81 # Acceleration of gravity
	nz  = len(z) # Amount of data in the vertical
	
	it = np.nonzero(z==zint)[0][0]
	dz = np.diff(z) 
	rho_int = rho[it] # Density at the reference height zint

    # Computation of BPE
	S = np.nan*np.zeros_like(z)
	S[it] = 0
	
	for i in range(0,it):
		sr=rho[i:it]+rho[i+1:it+1]
		sdz=dz[i:it]
		S[i] = -0.5*np.nansum(sr*sdz)
	for i in range(it+1,nz):
		sr=rho[it:i]+rho[it+1:i+1]
		sdz=dz[it:i]
		S[i] = 0.5*np.nansum(sr*sdz)

	BPE = (g*(z-zint)*rho_int)-(g*S)

	return BPE, z
