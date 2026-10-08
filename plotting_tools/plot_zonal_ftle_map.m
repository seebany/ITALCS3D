% function h2 = plot_zonal_ftle_map(modelstr, t0, plev, ht, lon, FTLE_Value, plane, tracer_pos)
% plots an FTLE map in the longitude-height plane at a given latitude of
% plane [deg N].
%
% S. Datta-Barua
% 7 Oct 2026 Preparing for publication to github.

function h2 = plot_zonal_ftle_map(modelstr, t0, plev, ht, lon, FTLE_Value, plane, tracer_pos)

fontsz = 8;
% First restrict to heights of interest to avoid edge effects and
% tropospheric LCSs which would dominate the colorbar.
minht = 30e3;
maxht = 100e3;

% Use this form for 3D LCSs.
nanrows = find(ht(:,:,1) < minht | ht(:,:,1) > maxht);
FTLE_Value(nanrows) = nan;

% Need theta to be 1 x nlon.
theta = lon';

% Next plot vs altitude.
th_matrix = repmat(theta, size(ht,1),1);
ht_matrix = squeeze(ht(:,:,1))/1e3;

% 2/16/26 SDB If lons have been wrapped, need to augment ht_matrix.
if max(theta) - min(theta) > 357.5 & size(th_matrix) ~= size(ht_matrix)
	ht_matrix = cat(2, ht_matrix, ht_matrix(:,1));
end

DT = delaunay(th_matrix, ht_matrix);
h2 = trisurf(DT, th_matrix, ht_matrix, FTLE_Value, 'FaceColor', 'interp', 'EdgeColor', 'none');
set(gca, 'CLim', [1e-4 1.5e-4])
h = xlabel('Longitude [deg]');
set(h,'FontSize', fontsz);
grid on
h = ylabel('Geopotential Height [km]');
set(h,'FontSize', fontsz);

view([0, 90])
h = title([modelstr ' FTLE Map in the zonal plane of ' num2str(plane) ' deg']);
set(h,'FontSize', fontsz);

hc = colorbar;
ylabel(hc, 'FTLE Value [1/s]');
axis([0 360 minht/1e3 maxht/1e3]);
