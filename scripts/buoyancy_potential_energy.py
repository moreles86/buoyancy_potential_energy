import numpy as np
from scipy import interpolate


"""
Function to calculate the BUOYANCY POTENTIAL ENERGY (BPE)
BPE quantifies the buoyancy potential energy required for a water parcel in hydrostatic equilibrium to be vertically displaced from a specific reference height.
"""
def buoyancy_potential_energy(rho, z, zint):
	"""
	Parameters:
    - Depth (z) and density (rho) data should go from the deepest to the shallowest and be column vectors
	- rho: Sigma-0 Potential Density Anomaly profile (kg*m^-3)
	- z: Heights (m) with negative units. The z-vector can be non-equidistant.
	- zint: Reference height (m)

	Return:
	- BPE: Buoyancy potential energy profile (J*m^-3)
	- z: Heights (m)
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
