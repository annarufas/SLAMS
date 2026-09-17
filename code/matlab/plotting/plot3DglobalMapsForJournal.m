
% My locations 
% lat = [31.6, 49, 79]; % BATS/OFP, PAP-SO, HAUSGARTEN
% lon = [-64.2, -16.5, 4.3];
% lat = [22.5, 0, 50]; % HOT/ALOHA, EqPac, OSP
% lon = [-158, -140, -145];
% lat = [22.5, 0, 50, 31.6, 49, 79];
% lon = [-158, -140, -145, -64.2, -16.5, 4.3];

lat = [22.5, 50, 31.6, 49];
lon = [-158, -145, -64.2, -16.5];

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.30 0.30])

% Orthographic projection centered over the Atlantic
% m_proj('ortho','lat',35','long',-35');
m_proj('ortho','lat',35','long',-90');

% Draw land in light grey
m_coast('patch',[0.8 0.8 0.8]);

% Hide grid labels
m_grid('linest','-','linewidth',0.20,'xticklabels',[],'yticklabels',[]);

% Hold on to add more elements
hold on

% Plot red dots
m_plot(lon, lat, 'o', ...
    'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerSize', 12, ...
    'LineWidth', 1);

% saveFigureInFolder('LOCALTS6','globalViewAtlantic')
% saveFigureInFolder('LOCALTS6','globalViewPacific')
saveFigureInFolder('.figures/LOCALTS6/','globalViewSome')