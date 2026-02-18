% Stokes test

particleradius = 250; % um
particledensity = 1.060; % g cm-3
waterRho = 1.026; % g cm-3
waterDynVisco = 0.0100; % g cm-1 s-1

[wwhite, wstokes, wcael] = settling_velocity_MATLAB(particleradius,...
    particledensity, waterRho, waterDynVisco);

% Formatted print
fprintf('wwhite  = %.1f m/d\nwstokes = %.1f m/d\nwcael   = %.1f m/d\n', ...
        wwhite, wstokes, wcael);

%% LOCAL FUNCTION

function [wwhite, wstokes, wcael] = settling_velocity_MATLAB(particleradius,...
    particledensity, waterRho, waterDynVisco)

GRAVITY_CNT = 9.81; % m s-2
SECONDS_PER_DAY = 24*3600;
rhoP      = 1d3 .* particledensity;        % kg m-3
rhoW      = 1d3 .* waterRho;               % kg m-3 
deltaRho  = abs(rhoP - rhoW);                         
diameterP = 2 .* particleradius .* 1d-6;  % m
mu        = waterDynVisco .* 0.1;        % g cm-1 s-1 -> g m-1 s-1

% Cael et al. (2021) approach
% McDonnell & Buesseler (2010) also argue against the unsuitability of Stokes
wcael = 129.*(2.*particleradius.*1d-3).^0.63d0; % m d-1, diameter in mm

% Stokes initial guess (magnitude)
wstokes = deltaRho .* (diameterP.^2) * GRAVITY_CNT ./ (18.*mu); % m s-1
reynoldsNo = wstokes.*rhoW.*diameterP ./ mu; % dimensionless 
	
% At higher Reynolds numbers
w = wstokes;      
if reynoldsNo >= 0.50

	for iVeloRecals = 1:50
		reynoldsNo = w.*rhoW.*diameterP ./ mu; % dimensionless

		if reynoldsNo <= 1000
			dragCoeff = (24 ./ reynoldsNo) + (6 ./ (1 + sqrt(reynoldsNo))) + 0.4;
		else
			dragCoeff = 0.4;
        end

		wNew = sqrt((4 .* deltaRho .* diameterP .* GRAVITY_CNT) / (3 .* rhoW .* dragCoeff)); % m s-1

		if all( abs(wstokes - wNew) <= max(1e-9, 1e-6.*wNew))
            break
        end
		w = wNew; % m s-1
    end
end

% Sign (negative: particle going upwards; positive: particle sinking)  
if (rhoP < rhoW)
	wwhite = -1.*w.*SECONDS_PER_DAY; % m d-1 (particle rising)
    wcael = -1.*wcael;
    wstokes = -1.*wstokes.*SECONDS_PER_DAY;
elseif (rhoP > rhoW)
	wwhite = w.*SECONDS_PER_DAY; 	 % m d-1 (particle sinking)
    wcael = wcael;
    wstokes = wstokes.*SECONDS_PER_DAY;
else 
	wwhite = 0; 					         % m d-1 (neutral buoyancy)
    wcael = 0;
    wstokes = 0;
end

end % function