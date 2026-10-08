% function h2 = plot_meridional_ftle_map(modelstr, t0, plev, ht, lat, FTLE_Value, plane)
% plots a the FTLE map as a function of latitude and height, in a
% meridional plane given by plane [deg E].
% FTLE_Value should be dimension nlev x nlat.
% 
% S. Datta-Barua
% 7 Oct 2026 Preparing for publication to github.

function h2 = plot_meridional_ftle_map(modelstr, t0, plev, ht, lat, FTLE_Value, plane)

% First restrict to heights of interest to avoid edge effects and tropospheric LCSs which are dominant.
minht = 30e3; % km
maxht = 100e3; % km

% Eliminate FTLE values that are outside the height range of interest.
nanrows = find(ht(:,:,1) < minht | ht(:,:,1) > maxht);
FTLE_Value(nanrows) = nan;

fontsz = 8;
nlat = numel(lat);
theta = lat;


% Next plot vs altitude.
% This is the form used for 3D LCS plotting.
th_matrix = repmat(theta', size(ht,1),1);
ht_matrix = squeeze(ht(:,:,1))/1e3;
DT = delaunay(th_matrix, ht_matrix);

h2 = trisurf(DT, th_matrix, ht_matrix, FTLE_Value, 'FaceColor', 'interp', 'EdgeColor', 'none');

set(gca, 'CLim', [1e-4 1.5e-4])
xlabel('Latitude [deg]')
grid on
ylabel('Geopotential Height [km]')
view([0, 90])
h = title([modelstr ' FTLE Map in the meridional plane of ' num2str(plane) ' deg']);
set(h,'FontSize', fontsz);
hc = colorbar;
ylabel(hc, 'FTLE Value [1/s]');
axis([0 90 minht/1e3 maxht/1e3])

