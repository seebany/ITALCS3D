% function [FTLE_Value] = compute_3d_ftle_map(t_datenum, lat, lon, gph, 
%	flow, t0, IntegrationTime, timesteps, computedirection_flag)
% This script generates 3D WACCMX+DART LCSs 
% and is adapted from the compute_ftle_map.m 2d version.
% Inputs are:
%	t_datenum - 1 x ntimes array of datenums of times at which the flows are provided.
%	lon - nlon x 1 array of longitudes [deg].
%	lat - nlat x 1 array of latitudes [deg].
%	gph - nlon x nlat x nlev x ntime array of heights [m], not pressure levels [hPa?].
%	flow - ntimes x npts matrix of [time, u1, v1, w1,..., u_npts, v_npts, w_npts].
% 	t0 - 1 x 1 datenum of start time.
%	IntegrationTime - 1 x 1 number of seconds of integration from time t0 to tf.
%	num_dt - integer number of subdivisions of the time between each flow field.
%	timedirection_flag - 1 for forward time integration, 0 for backward.
%
% S. Datta-Barua
% 23 May 2025
% 8 Sept 2025 To fix the tracer calculation in the 3rd dimension make the z 
% 	dimension a regularly spaced grid in plev index units.
% 23 Jul 2026 Input flag about whether to transform J to Jtilde (the Cartesian form).

function [FTLE, tracer_pos] = compute_3d_ftle_map_pleiades(t_datenum, lon, lat, gph, flow, t0, IntegrationTime, num_dt, timedirection_flag, FTLE_flag, tracer_flag, tracer_start, vel_nonzero_flag, J_in_Cartesian_flag)

% Is the domain periodic?
period_flag = 0; % Ningchao used 1.
% Is the domain Euclidean (Cartesian)?
Euclidean_flag = 0; %Ningchao used -1.
% Along which dimensions is the domain closed? [x/lon, y/lat, z/lev]
CloseDomain_flag = [1 1 0];
%computedirection_flag = 1;
dim_active_flag = [1 1 1];
%FTLE_flag = 1;
%tracer_flag = 0;
%tracer_start = [];
FDFTLE_index = 0;
boundaries = [];
Delta_xyz = [];
fdmin_xyz = [];
fdmax_xyz = [];
% Ok to use +/-90 deg latitude in the boundaries here.
xmin = min(lon); % 0 deg for WACCMX
xmax = max(lon); % 357.5 deg for WACCMX, 360 after wrapping 2/12/26 SDB.
%nx = numel(lon); %145;
ymin = min(lat); % -90 deg for WACCMX
ymax = max(lat); % 90 deg for WACCM
%ny = numel(lat); %72;
% 11/13/25 SDB reducing number of points to troubleshoot running of qsub job.
ny = numel(lat);%/8; %9 points, a factor of 8 reduction
nx = numel(lon);%/8; % 18 points, a factor of 8 reduction
% Pressure levels are not regularly spaced, which makes dz incorrect when it is computed and used.
%zmin = min(lev);
%zmax = max(lev);
%zmin = min(log(lev));
%zmax = max(log(lev));
%nz = numel(lev);
%keyboard

secperday = 24*60*60;
%IntegrationTime = 15*3600; % Time interval in the time units of the data (usually s)
% that  is tf - t0.
tmax = t_datenum(end)*secperday + 86400*100; % I think Ningchao said set this to a large number so
% that it doesn't try to wrap around cyclically in time.
TimeStamp = mean(diff(t_datenum))*24*3600; % Data cadence $\Delta t$ in time units of the data (usually s).
%timesteps = 10; %30; % How many dt the data cadence should be subdivided into.

t_initial = (t0 - t_datenum(1))*secperday;%9*secperday:9*secperday %600*6*0:600:600*6*69;
tic
[FTLE,tracer_pos] = FTLE_FunctionCalls(flow, xmin, xmax, nx, ...
	ymin, ymax, ny, ... %zmin, zmax, nz, ...
	gph, ...%lev, ...
	tmax,TimeStamp, num_dt, t_initial, IntegrationTime, ...
	CloseDomain_flag, period_flag, timedirection_flag, ...
	Euclidean_flag, dim_active_flag, FTLE_flag, ...
	tracer_flag, tracer_start, ...
	FDFTLE_index, Delta_xyz, fdmin_xyz, fdmax_xyz, ...
	vel_nonzero_flag, J_in_Cartesian_flag);
%FTLE_Value = FTLEmap_LAB_WithDIRECTION(flow,xmin,xmax,nx,...
%    ymin,ymax,ny,...
%    zmin,zmax,nz,...
%    t_initial,IntegrationTime,tmax,...
%    TimeStamp,timesteps,period_flag,...
%    Euclidean_flag,CloseDomain_flag,computedirection_flag);
toc
return

