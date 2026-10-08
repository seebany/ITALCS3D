% This script generates WACCMX+DART 3D FTLE maps and tracers using ITALCS3D
% and saves the results in netcdf files to be used for plotting.
%
% S. Datta-Barua
% 7 Oct 2026 Preparing for publication to github.

clear all
close all

%% Select settings here.  No other section should need editing.
% Specify the directory path containing the WACCMX+DART netcdf files.
input_dir = '/path/to/waccmxdart/data/';

% To save FTLE map and tracers as .nc file, specify output directory.
ncdir = '/path/to/ncfiles/outputs/';
outfiledir = [ncdir filesep datestr(today, 'yymmdd') '_runs/'];
if ~exist(outfiledir, 'dir')
    mkdir(outfiledir)
end

% Specify the UTC year, month, day of the start time of the event you want
% to analyze.
yymmddstr = '090125';%'090129';%
integration_hours = 24; %168; If generating data for Hovmoller Fig 8.

vel_nonzero_flag = [1 1 1]; %3x1 array of true/false for [u, v, w]. Use true for nonzero components.
J_in_Cartesian_flag = 1; % True if you want to transform to Jtilde, otherwise false.
FTLE_flag = 1; % 0 if you don't want to compute them.


% These options are fixed at these values in the 2026 ITALCS3D manuscript.
modelstr = 'WACCMX+DART';
highres_flag = 2; % 0 for standard resolution (6-hr assimilations used in the
% first ITALCS project), 1 for high (0.25 deg) res, 2 if you want the hourly assimilations at standard spatial resolution.
tracer_flag = 1;
plev_flag = 0; % True if you want grid to be plev, false if you want height.
filter_flag = 0;
timedirection_flag = 1;%-1;
quiver_flag = 1;

% The rest of the script loads the WACCMX+DART nc files, runs ITALCS3D, and
% saves the outputs in nc files.
%% --------------------------
% No need to change anything below here.
if J_in_Cartesian_flag
    Jstr = '_Jcart';
else
    Jstr = '';
end

if filter_flag
    filterstr = '_filtered';
else
    filterstr = '';%'_backward';
end
if timedirection_flag == 1
    timedirstr = '';
elseif timedirection_flag == -1
    timedirstr = '_backward';
end

% Select WACCMX+DART type of run.
if highres_flag == 2
    filterstr = 'hrlyassim';
elseif highres_flag == 1
    filterstr = 'highres';
else
    filterstr = '6hrassim';
end

% Initialize tracers along a latitude.
dlon = 10;
if tracer_flag
    tracer_start(:,1) = [0:dlon:355]; % deg lon
    tracer_start(:,2) = 65; % deg lat
    % SDB 3/11/26 Using this to check apparent vertical oscillation.
    % tracer_lons(:,1) = [0:dlon:359]; % deg lon
    %tracer_lons(:,1) = [200:dlon:300]; % deg lon
    %[tracer_grid_lon,tracer_grid_lat] = ndgrid(tracer_lons(:,1),[65]);% 75 85]);
    %tracer_start(:,1) = reshape(tracer_grid_lon, numel(tracer_grid_lon),1); % deg lon
    %tracer_start(:,2) = reshape(tracer_grid_lat, numel(tracer_grid_lat),1); % deg lon
else
    tracer_start = [];
end

year = 2009;
mon = 1;
day = str2num(yymmddstr(5:6));
hr = 0;
minute = 0;
sec = 0;
if tracer_flag
    lat_plane = tracer_start(1,2); % latitude, deg N.
else
    lat_plane = 65;
end
lon_plane = mod(-45, 360); % longitude, positive deg E.
ht_flag = '';%'gph';%'gmh';
tarht = 90e3; % m
if tracer_flag
    tracer_start(:,3) = tarht; % m
end

tau = integration_hours*3600; % s
t0 = datenum([year, mon, day, hr, minute, sec]);

%--------------------------
%% First select nc files to be used.

date = datenum([year, mon, day]);
t_load_datenum = date + [0:integration_hours]/24;
if highres_flag == 2
    for j = 1:numel(t_load_datenum)
        % Hourly assimilation of data.
        ncfilename{j} = ['waccmx_dart_h2o_' num2str(year) '.cam.h1.ensmean.' datestr(t_load_datenum(j), 'yyyy-mm-dd-') sprintf('%05d', round(rem(t_load_datenum(j),1)*24*3600,0)) '.nc'];
        filelist{j} = [input_dir ncfilename{j}];
    end % for j = 1:numel(t_load_datenum)
end % if hourly data set

%% Get the flow data from the desired files(s).
[t_datenum, lon, lat, plev, flow, ht, u, v, w] = waccm_nc2flowmat3d(filelist, filter_flag, plev_flag, highres_flag, vel_nonzero_flag);

%% Compute the FTLE Values.
num_dt = 10; % number of subdivisions dt of the data cadence Delta t.
[FTLE_Value, tracer_pos] = compute_3d_ftle_map(t_datenum, lon, lat, ht, flow, t0, tau, num_dt, timedirection_flag, FTLE_flag, tracer_flag, tracer_start, vel_nonzero_flag, J_in_Cartesian_flag);


%% save netcdf file of outputs. 
run_conditions= [ modelstr '_ITALCS3Datmo_' datestr(t0(1), 'yyyymmdd') '_' ...
    num2str(integration_hours(1)) 'hr_' ... %regionstr '_' ...
    ht_flag num2str(tarht/1e3) 'km_' filterstr '_uvw' num2str(vel_nonzero_flag(1)), ...
    num2str(vel_nonzero_flag(2)), num2str(vel_nonzero_flag(3)) Jstr];
fname=[outfiledir run_conditions 'FTLE.nc'];

% Clear any temp file.
if exist('tempnc.nc', 'file')
    !rm tempnc.nc
end

% write flags and run configurations to file.
nccreate('tempnc.nc', 'modelstr', 'Datatype', 'char', 'Dimensions', {'nchar_modelstr' numel(modelstr)});
ncwrite('tempnc.nc', 'modelstr', modelstr);
nccreate('tempnc.nc', 'integration_hours');
ncwrite('tempnc.nc', 'integration_hours', integration_hours);
nccreate('tempnc.nc', 'filterstr', 'Datatype', 'char', 'Dimensions', {'nchar_filterstr' numel(filterstr)});
ncwrite('tempnc.nc', 'filterstr', filterstr);
nccreate('tempnc.nc', 'vel_nonzero_flag', 'Dimensions', {'nspatial_dim' 3});
ncwrite('tempnc.nc', 'vel_nonzero_flag', vel_nonzero_flag);
nccreate('tempnc.nc', 't0');
ncwrite('tempnc.nc', 't0', t0);
nccreate('tempnc.nc', 'tarht');
ncwrite('tempnc.nc', 'tarht', tarht);
nccreate('tempnc.nc', 'FTLE_flag')
ncwrite('tempnc.nc', 'FTLE_flag', FTLE_flag)
nccreate('tempnc.nc', 'tracer_flag')
ncwrite('tempnc.nc', 'tracer_flag', tracer_flag)
nccreate('tempnc.nc', 'highres_flag')
ncwrite('tempnc.nc', 'highres_flag', FTLE_flag)

% write dimensions and outputs.
sz_lat=size(lat);
sz_lon=size(lon);
sz_plev = size(plev);
sz_ht = size(ht);
sz_t0=size(t0);%YYYYMMDDHH);
sz_datenum = size(t_datenum);
nccreate('tempnc.nc','lat','Dimensions',{'nlat' sz_lat(1)});%,'Format','classic');
ncwrite('tempnc.nc','lat',lat);
nccreate('tempnc.nc','lon','Dimensions',{'nlon' sz_lon(1)});%,'Format','classic');
ncwrite('tempnc.nc','lon',lon);
nccreate('tempnc.nc', 'plev', 'Dimensions', {'nlev' sz_plev(1)});
ncwrite('tempnc.nc', 'plev', plev);
nccreate('tempnc.nc', 'ht', 'Dimensions', {'nlon' sz_lon(1) 'nlat' sz_lat(1) 'nlev' sz_plev(1) 'epoch' sz_datenum(1)});
ht_wrapped = cat(1, ht, ht(1,:,:,:));
ncwrite('tempnc.nc', 'ht', ht_wrapped);

% write FTLE values if they were computed.
if FTLE_flag
    if timedirection_flag
        sz_ftle=size(FTLE_Value);
        nccreate('tempnc.nc','FTLE','Dimensions',{'nlat' sz_ftle(1) 'nlon' sz_ftle(2) 'nht' sz_ftle(3)});%,'Format','classic');
        ncwrite('tempnc.nc','FTLE',FTLE_Value);
    else
        nccreate('tempnc.nc','backFTLE','Dimensions',{'nlat' sz_ftle(1) 'nlon' sz_ftle(2) 'nt0' sz_t0(1)});%,'Format','classic');
        ncwrite('tempnc.nc','backFTLE',backFTLE);
        nccreate('tempnc.nc','PRA','Dimensions',{'nlat' sz_ftle(1) 'nlon' sz_ftle(2) 'nt0' sz_t0(1)});%,'Format','classic');
        ncwrite('tempnc.nc','PRA',PRA);
    end
end
% write tracer positions if they were computed.
if tracer_flag
    sz_tracer = size(tracer_pos);
    nccreate('tempnc.nc','tracer_pos','Dimensions',{'epoch' sz_tracer(1) 'nspatial_dim' sz_tracer(2) 'ntracers' sz_tracer(3)});%,'Format','classic');
    ncwrite('tempnc.nc','tracer_pos',tracer_pos);
end

eval(['!mv tempnc.nc ' fname]);
ncdisp(fname);

disp(['data saved to: ', fname])


