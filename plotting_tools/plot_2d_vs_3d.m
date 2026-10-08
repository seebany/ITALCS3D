% This script loads and plots subfigures to compare the horizontal flows
% using ITALCS (2D) and ITALCS3Datmo with the hourly assimilated
% WACCMX+DART.
%
% S. Datta-Barua
% 7 Oct 2026 Preparing for publication to github.

clear

FTLE_flag = 1;
tracer_flag = 1;

% Specify the path where the ITALCS3D run result netcdf files were placed.
run_folder = '/path/to/italcs3d/output/';%'/Users/seebany/Dropbox/projects/atmoLCS/atmoLCS_2/from_pleiades/data_for_publication/2026_ITALCS3D/'; %

% FTLE run file for Fig 3a.
fname{1} = [ run_folder 'WACCMX+DART_ITALCS_200901_24hr_90km_hrlyassimFTLE.nc'];
% File for Fig 3b.
fname{2} = [ run_folder 'WACCMX+DART_ITALCS3Datmo_20090125_24hr_90km_hrlyassim_uvw110FTLE.nc'];
% File for Fig 3c.
fname{3} = [ run_folder 'WACCMX+DART_ITALCS3Datmo_20090125_24hr_90km_hrlyassim_uvw110_JcartFTLE.nc'];
% File for Fig 3d.
fname{4} = [ run_folder 'WACCMX+DART_ITALCS3Datmo_20090125_24hr_90km_hrlyassim_uvw111_JcartFTLE.nc'];

%---------------------------
% Specify whether to use J (set flag to 0) or Jtilde (set flag to 1) for
% Figs 3 [a b c d].
J_in_Cartesian_flag = [0 0 1 1]; % True if you want to transform to Jtilde, otherwise false.

prefixstr = 'abcd';

tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');
for subplot_idx = 1:4
    nexttile;
    lat = ncread(fname{subplot_idx}, 'lat');
    lon = ncread(fname{subplot_idx}, 'lon');
    clear FTLE_Value tracer_pos
    FTLE_Value = ncread(fname{subplot_idx}, 'FTLE');

    try
    	tracer_pos = ncread(fname{subplot_idx}, 'tracer_pos');
    catch
    	tracerlat = ncread(fname{subplot_idx}, 'tracerlat');
    	tracerlon = ncread(fname{subplot_idx}, 'tracerlon');
    	tracer_pos(:,2,:) = tracerlat;
    	tracer_pos(:,1,:) = tracerlon;
    end

    if ~exist('modelstr', 'var')
    	modelstr = 'WACCMX+DART';
    end
    if ~exist('t0', 'var')
        year = 2009;
        mon = 1;
        day = 25;
        hr = 0;
        minute = 0;
        sec = 0;
        t0 = datenum([year, mon, day, hr, minute, sec]);
    end
    if ~exist('trueht', 'var')
    	trueht = 90; %km
    	tarht = 90e3; %m
    end
    if ~exist('regionstr', 'var')
    	regionstr = 'North Pole';
    	origin = [90, 0]; % plot north pole at the center, 0 UT local time at the bottom
    end
    try
    	ncread(fname{subplot_idx}, 'filterstr');
    catch
    	filterstr = '';
    end

    output_dir = [getenv('HOME') '/mfigs/lcs/' modelstr '/' ...
    	datestr(now,'yymmdd') '/'];
    if ~exist(output_dir,'dir')
    	mkdir(output_dir)
    end
    try
        integration_hours = ncread(fname{subplot_idx}, 'integration_hours');
    catch
        integration_hours = [24];%[24; 24; 24];
    end
    try
    	FTLE_flag = ncread(fname{subplot_idx}, 'FTLE_flag');
    catch
    	FTLE_flag = 1;
    end

    % Fig a: hourly WACCMX+DART data, run with ITALCS (2D) called by waccm_lcs_script.m.
    if FTLE_flag % Plot a polar map with the FTLE and continents, and colorbar.
    	if ndims(FTLE_Value) < 3 % FTLE_Value [96 x 145 x 126] nlat x nlon x nlev
    	    h = plot_polar_ftle_map(['(' prefixstr(subplot_idx) ') ' modelstr], t0(1), lat, lon, trueht, FTLE_Value(:,:,1), regionstr, integration_hours(1), J_in_Cartesian_flag(subplot_idx));
    	else
    	    try
        		ht = ncread(fname{subplot_idx}, 'ht');
        		htrow = find(abs(ht(1,1,:,1) - tarht) == min(abs(ht(1,1,:,1) - tarht)))
    	    end
    	    FTLE_Value_htslice = squeeze(FTLE_Value(:,:,htrow));
            h = plot_polar_ftle_map(['(' prefixstr(subplot_idx) ') ' modelstr], t0(1), lat, lon, tarht, FTLE_Value_htslice, regionstr, integration_hours, J_in_Cartesian_flag(subplot_idx));
        end
    else % Plot a plain map with continents.
    	h = plot_polar_map(modelstr, t0(1), lat, lon, tarht, regionstr, integration_hours(i));
    end
    if tracer_flag
    	% Input arguments should be [ntimes x ntracers].
    	hA = plot_tracers_on_map(squeeze(tracer_pos(:,2,:)), squeeze(tracer_pos(:,1,:)));
    end
end

clear figfile
figfile = [output_dir datestr(now,'yymmdd') '_' ...
    modelstr '_' datestr(t0(1), 'yymmdd_hhMM') '_' ...
    num2str(integration_hours(1)) 'hr_' ...
    num2str(tarht/1e3) 'km_ftle' filterstr '2d-vs-3d.fig'];
saveas(gcf,[figfile(1:end-3) 'png'], 'png')
close
