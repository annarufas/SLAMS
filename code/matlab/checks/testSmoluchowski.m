Nzc = 16;
timeStep = 3600; % s
binVolume = 5; % m3
encounterKernelByZc = zeros(Nzc,1);
probEncPerPair = zeros(Nzc,1);
rateEncPerPair = zeros(Nzc,1);
zooDensityByZc = zeros(Nzc,1);
zooCountsByZc = zeros(Nzc,1);

encounterKernelByZc(:) = [8.03e-09, 8.09e-09, 8.16e-09, 8.25e-09, 8.38e-09, 8.56e-09,...
         8.77e-09, 9.06e-09, 9.47e-09, 9.98e-09, 1.07e-08, 1.17e-08,...
         1.29e-08, 1.47e-08, 1.73e-08, 3.14e-08]; % m3 s-1
zooDensityByZc(:) = [1e6, 6.3e5, 4e5, 2.5e5, 1.6e5, 1e5, 6.3e4, 4e4,...
          2.5e4, 1.6e4, 1e4, 6.3e3, 4e3, 2.5e3, 1.6e3, 1e3]; % num. m-3 (density)
zooCountsByZc(:) = zooDensityByZc(:) .* binVolume; % counts in the bin

rateEncPerPair(:) = encounterKernelByZc(:) .* zooDensityByZc; % s-1
probEncPerPair(:) = encounterKernelByZc(:) .* timeStep ./ binVolume; % unitless

% Route A (same as lambda in D)
probPairFail_est = exp( - sum( zooCountsByZc(:) .* probEncPerPair(:) ) ); % exact equivalence
probPairEncAtLeastOnce_est = 1 - probPairFail_est;

% Route B - Binomial form
probPairFail_prod = prod( (1 - probEncPerPair(:)) .^(zooCountsByZc(:)) );
probPairEncAtLeastOnce_prod = 1 - probPairFail_prod;

% Route C - equivalent to B, stable, written with logs
probPairFail_log = exp( sum( zooCountsByZc(:) .* log1p( -probEncPerPair(:) ) ) );
probPairEncAtLeastOnce_log = 1 - probPairFail_log;

% Route D - the correct form
numExpectedEncPerParticleTotal = zeros(Nzc,1);
numExpectedEncPerParticleTotal(:) = zooDensityByZc(:) .* encounterKernelByZc(:) .* timeStep; % lambda
probPairFail_poiss = exp(sum(-numExpectedEncPerParticleTotal(:)));
probPairEncAtLeastOnce_poiss = 1 - probPairFail_poiss; 

solidVolume = 456.89;
npp = 43;
primParticleMatVol = solidVolume / npp;
