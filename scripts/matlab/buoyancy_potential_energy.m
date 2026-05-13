% Function to calculate the BUOYANCY POTENTIAL ENERGY (BPE)
% BPE quantifies the buoyancy potential energy required for a water parcel in hydrostatic equilibrium to be vertically displaced from a specific reference height.

function [BPE,z] = buoyancy_potential_energy(rho, z, zint)
% Parameters:
%   - Depth (z) and density (rho) data should go from the deepest to the shallowest and be column vectors
%   - rho: Sigma-0 Potential Density Anomaly profile (kg*m^-3)
%   - z: Heights (m) with negative units. The z-vector can be non-equidistant.
%   - zint: Reference height (m)
% Return:
%   - BPE: Buoyancy potential energy profile (J*m^-3)
%   - z: Heights (m)

% If there is no data at zint, it is interpolated
if sum(z(z == zint), 'omitnan') == 0
    f = @(x) interp1(z, rho, x, 'linear');
    it = find(z >= zint, 1);
    z = [z(1:it-1); zint; z(it:end)];
    rho = [rho(1:it-1); f(zint); rho(it:end)];
end

g = 9.81; % Acceleration of gravity
nz = length(rho); % Amount of data in the vertical

it = find(z == zint, 1); 
dz = diff(z); 
rho_int = rho(it); % Density at the reference height zint

% Computation of BPE
S = NaN(nz, 1);
S(it) = 0;

for i = 1:1:it-1
    sr = rho(i:it-1) + rho(i+1:it);
    sdz = dz(i:it-1);
    S(i) = -0.5*sum(sr.*sdz);
end
for i = it+1:1:nz
    sr = rho(it:i-1) + rho(it+1:i);
    sdz = dz(it:i-1);
    S(i) = 0.5*sum(sr.*sdz);
end

BPE = (g*(z-zint)*rho_int)-(g*S);

end

