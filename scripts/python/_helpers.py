import gsw
import numpy as np

def compute_ild_mld(p, SA, CT, lat, dT=0.2):
    """
    Compute isothermal layer depth (ILD) and mixed layer depth (MLD)

    Parameters
    ----------
    p: Sea pressure (dbar)
    SA: Absolute Salinity (g kg-1)
    CT: Conservative Temperature (°C)
    lat: Latitude (decimal degrees)
    dT: Temperature threshold relative to 10 m reference depth (°C, default 0.2)

    Returns
    -------
    ild_z: Isothermal layer depth (m)
    mld_z: Mixed layer depth (m)
    ild_p: Pressure at ILD (dbar)
    mld_p: Pressure at MLD (dbar)
    ild_at_bottom: True if ILD corresponds to the deepest available level
    mld_at_bottom: True if MLD corresponds to the deepest available level

    Notes
    -----
    ILD calculation follows https://doi.org/10.1029/2004JC002378
    MLD calculation follows https://doi.org/10.1029/2006JC003953
    """
    p = np.asarray(p)
    SA = np.asarray(SA)
    CT = np.asarray(CT)

    z = -gsw.z_from_p(p, lat)

    if z.min() > 10 or z.max() < 10:
        return np.nan, np.nan, np.nan, np.nan, False, False

    T10 = np.interp(10, z, CT)
    S10 = np.interp(10, z, SA)

    # ILD
    T_threshold = T10 - dT
    fT = CT - T_threshold

    ild_z = np.nan
    ild_p = np.nan
    ild_at_bottom = False

    for i in range(1, len(z)):
        if fT[i-1] > 0 and fT[i] <= 0:
            w = fT[i-1] / (fT[i-1] - fT[i])
            ild_z = z[i-1] + w * (z[i] - z[i-1])
            ild_p = p[i-1] + w * (p[i] - p[i-1])
            break

    if np.isnan(ild_z):
        ild_z = z.max()
        ild_p = p[-1]
        ild_at_bottom = True

    # MLD
    sigma = gsw.sigma0(SA, CT)
    sigma10 = gsw.sigma0(S10, T10)
    sigma_ref = gsw.sigma0(S10, T10 - dT)

    delta_sigma = sigma_ref - sigma10
    fRho = (sigma - sigma10) - delta_sigma

    mld_z = np.nan
    mld_p = np.nan
    mld_at_bottom = False

    for i in range(1, len(z)):
        if fRho[i-1] < 0 and fRho[i] >= 0:
            w = fRho[i-1] / (fRho[i-1] - fRho[i])
            mld_z = z[i-1] + w * (z[i] - z[i-1])
            mld_p = p[i-1] + w * (p[i] - p[i-1])
            break

    if np.isnan(mld_z):
        mld_z = z.max()
        mld_p = p[-1]
        mld_at_bottom = True

    return ild_z, mld_z, ild_p, mld_p, ild_at_bottom, mld_at_bottom