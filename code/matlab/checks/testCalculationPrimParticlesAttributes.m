
% Chemical and physical constants
MOLAR_MASS_CARBON    = 12.011; % g mol-1
MOLAR_MASS_CACO3     = 100.1; % g mol-1
MOLAR_MASS_OPAL      = 67.3; % g mol-1
MOLAR_MASS_CLAY      = 389.34;
RHO_ORGMATTER        = 1.06;
RHO_TEP              = 0.800;
RHO_CALCITE          = 2.72; % g cm-3
RHO_OPAL             = 1.90; % g cm-3
RHO_CLAY             = 2.70;
C_FRAC_IN_OM         = 0.54; % after Anderson (1995), agrees with Lam et al. (2011)
C_FRAC_IN_TEP        = 0.95;
GRAVITY_CNT          = 9.81; % m s-2

% Mineral quotas PFTs
Si2C_diat = 0.13; % mol Si (mol C)-1
Calc2C_cocco_max = 0.50; % mol PIC (mol POC)-1
clay_quota = 1.0e-13; % mol clay part-1

% Carbon quotas PFTs (min to max)
C_quota_diat = [10e-12, 5000e-12]; % mol C cell
C_quota_flagel = [400e-12, 10000e-12]; % mol C cell-1
C_quota_cocco = [0.8e-12, 5e-12]; % mol C cell-1
C_quota_pico = [0.0020e-12, 0.10e-12]; % mol C cell-1
C_quota_tep = [10e-12, 1000e-12];

% Output arrays
particleTypes = {'diat', 'flagel', 'cocco', 'pico'};
nTypes = numel(particleTypes);
densityPft         = zeros(nTypes,1);
diameterPft        = zeros(nTypes,2);
solidVolumePft     = zeros(nTypes,2);
sinkingVelocityPft = zeros(nTypes,2);

for iPass = 1:2

    C_quota_part = [C_quota_diat(iPass),C_quota_flagel(iPass),C_quota_cocco(iPass),C_quota_pico(iPass)];
    Si_quota_part = [C_quota_diat(iPass)*Si2C_diat, 0, 0, 0]; % mol Si cell-1
    CaCO3_quota_part = [0, 0, C_quota_cocco(iPass)*Calc2C_cocco_max, 0]; % mol CaCO3 cell-1

    for iPart = 1:nTypes
        massOrgC = C_quota_part(iPart) * MOLAR_MASS_CARBON; % g organic C per particle
        massOrgMatter = massOrgC / C_FRAC_IN_OM; % g organic matter per particle
        massOpal = Si_quota_part(iPart) * MOLAR_MASS_OPAL; % g bSi per particle
        massCalcite = CaCO3_quota_part(iPart) * MOLAR_MASS_CACO3; % g CaCO3 per particle
        mass = massOrgMatter + massOpal + massCalcite;
        volumeCarbon = massOrgMatter/RHO_ORGMATTER;
        volumeOpal = massOpal/RHO_OPAL;
        volumeCalcite = massCalcite/RHO_CALCITE;
        solidVolumePft(iPart,iPass) = 1e12.*(volumeCarbon + volumeOpal + volumeCalcite); % cm3 --> um3
        radius = ((3.*solidVolumePft(iPart,iPass))./(4.*pi)).^(1/3); % um
        diameterPft(iPart,iPass) = 2.*radius;
        densityPft(iPart) = (mass./solidVolumePft(iPart,iPass)) .* 1d12; % g cm-3

        % Stokes
	    rhoP = 1d3 .* densityPft(iPart); % kg m-3
	    rhoW = 1d3 .* 1.0262; % kg m-3 (mean value global first 200 m)
	    deltaRho = abs(rhoP - rhoW);
	    diameterP = 2d0 .* radius.*1d-6; % m
	    mu = 0.0134 .* 0.1d0; % kg m-1 s-1         
	    w = deltaRho.*(diameterP.^2).*GRAVITY_CNT ./ (18.*mu); % m s-1
	    reynoldsNo = w.*rhoW.*diameterP ./ mu; % dimensionless 
    	sinkingVelocityPft(iPart,iPass) = w.*(3600*24);

    end
    
end

% Clay
massClay = clay_quota .* MOLAR_MASS_CLAY;
solidVolumeClay = 1e12.*(massClay/RHO_CLAY); % cm3 --> um3
radiusClay = ((3.*solidVolumeClay)./(4.*pi)).^(1/3); % um
diameterClay = 2.*radiusClay;
densityClay = (massClay./solidVolumeClay) .* 1d12; % g cm-3
rhoP = 1d3 .* densityClay; % kg m-3
rhoW = 1d3 .* 1.0262; % kg m-3 (mean value global first 200 m)
deltaRho = abs(rhoP - rhoW);
diameterP = 2d0 .* radiusClay.*1d-6; % m
mu = 0.0134 .* 0.1d0; % kg m-1 s-1         
w = deltaRho.*(diameterP.^2).*GRAVITY_CNT ./ (18.*mu); % m s-1
sinkingVelocityClay = w.*(3600*24);

% TEP
densityTep = zeros(2,1);
solidVolumeTep = zeros(2,1);
diameterTep = zeros(2,1);
sinkingVelocityTep = zeros(2,1);
for iPass = 1:2
    massTepC = C_quota_tep(iPass) * MOLAR_MASS_CARBON;
    massTep = massTepC / C_FRAC_IN_TEP; % g organic matter per TEP
    volumeTep = massTep/RHO_TEP;
    solidVolumeTep(iPass) = 1e12.*volumeTep; % cm3 --> um3
    radius = ((3.*solidVolumeTep(iPass))./(4.*pi)).^(1/3); % um
    diameterTep(iPass) = 2.*radius;
    densityTep(iPass) = (massTep./solidVolumeTep(iPass)) .* 1d12; % g cm-3
    rhoP = 1d3 .* densityTep(iPass); % kg m-3
    rhoW = 1d3 .* 1.0262; % kg m-3 (mean value global first 200 m)
    deltaRho = abs(rhoP - rhoW);
    diameterP = 2d0 .* radius.*1d-6; % m
    mu = 0.0134 .* 0.1d0; % kg m-1 s-1         
    w = deltaRho.*(diameterP.^2).*GRAVITY_CNT ./ (18.*mu); % m s-1
    reynoldsNo = w.*rhoW.*diameterP ./ mu; % dimensionless 
    sinkingVelocityTep(iPass) = w.*(3600*24);
end
