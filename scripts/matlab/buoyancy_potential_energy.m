function [BPE_mass, BPE_volume, p] = buoyancy_potential_energy(p, SA, CT, Peq)
% Buoyancy Potential Energy (BPE) from accurate energetic formulations
% BPE quantifies the buoyancy potential energy required for a water parcel in hydrostatic equilibrium to be vertically displaced from a specific reference height zint.

% Parameters
% ----------
% p: Sea pressure (dbar)
% SA: Absolute Salinity (g kg-1)
% CT: Conservative Temperature (°C)
% Peq: Pressure of the equilibrium position of the parcel (dbar)

% Returns
% ----------
% BPE_mass: BPE per unit mass (J kg-1)
% BPE_volume: BPE per unit volume (J m-3)
% p: Sea pressure (dbar)

% Notes
% ----------
% TEOS-10 implementation for BPE follows https://doi.org/10.5194/egusphere-2025-3359-RC1

% If there is no data at Peq, it is interpolated in SA and CT
if sum(find(p==Peq),'omitnan') == 0
    p_ext=[Peq; p]; p_ext=sort(p_ext,'ascend');

    % Interpolate SA and CT at Peq
    SA_ext= interp1(p,SA,p_ext,'linear','extrap');
    CT_ext= interp1(p,CT,p_ext,'linear','extrap');

    % Rename variables
    p=p_ext;
    SA=SA_ext;
    CT=CT_ext;
end

iPt=find(p==Peq);

% Computation of BPE
BPE_mass=NaN(length(p),1);
BPE_volume=NaN(length(p),1);
specvol=NaN(length(p),1); % Specific volume (m3 kg-1)

for k=1:length(p)
    iP2=k;
    e1 = gsw_enthalpy(SA(iPt),CT(iPt),p(iP2)); % Specific enthalpy [J/kg]
    e2 = gsw_enthalpy(SA(iPt),CT(iPt),p(iPt));
    gsdh = gsw_geo_strf_dyn_height(SA,CT,p,p(iPt));
    e3 = gsw_enthalpy(35.16504,0,p(iP2));
    e4 = gsw_enthalpy(35.16504,0,p(iPt));
    specvol(k) = gsw_specvol(SA(iPt),CT(iPt),p(iP2));

    BPE_mass(k)=(e1-e2+gsdh(iP2)-e3+e4);
    BPE_volume(k)=BPE_mass(k)/specvol(k);
end

end

