% function [t_datenum, lon_wrapped, lat, p, vel_ts, GPH, u, v, OMEGA] = waccm_nc2flowmat3d(filename, filter_flag, plev_flag, highres_flag, vel_nonzero_flag)
% Function to load WACCM nc files and convert to a mat file that can be
% used by ITALCS3D. Function is modified from waccm_nc2flowmat.m.
%
% Input string filename is the full path to the WACCM velocity data .nc file, e.g.,
% Output vel_ts is a matrix of dimension [ntimes x 3(number gridpoints)]. 
% Each row of the output matrix vel_ts is of the form
% [ts0, u0, v0, w0, u1, v1, w1, ..., ts1, u0, v0, w0, u1, v1, w1, ...]
% It is a (length of time stamps) by 1+ 3*(1 + nx*ny) matrix, where nx and ny
% are the number of gridpoints along the x-axis (longitude) and y-axis
% (latitude). Note that the first and last gridpoints will be the same in
% order to close the domain and have particles wrap around the globe
% correctly.
%
% S. Datta-Barua
% 7 Oct 2026

function [t_datenum, lon_wrapped, lat, p, vel_ts, GPH, u, v, OMEGA] = waccm_nc2flowmat3d(filename, filter_flag, plev_flag, highres_flag, vel_nonzero_flag)

% Read the .nc file and return the variables desired.
% lon is [nlons x 1]
% lat is [nlats x 1]
% u is [nlons x nlats x nlev x ntimes]
% v is [nlons x nlats x nlev x ntimes]
% GPH is [nlons x nlats x nlev x ntimes]
% p is double [nlev x 1]
% t_datenum is [1 x ntimes]
% OMEGA is [nlons x nlats x nlev x ntimes]
% T is [nlons x nlats x nlev x ntimes]
[lon, lat, u, v, GPH, p, t_datenum, OMEGA, T] = waccm3d_nc2mat(filename, highres_flag);
% SDB 2/17/26
% Test out zero velocity in the 3rd dimension with ITALCS3D.
% SDB 3/10/26 The zero vertical velocity matches 2D almost perfectly when
% the J matrix is not transformed to Jtilde.
% So we can comment out and stick to 3D.
if ~vel_nonzero_flag(3)
	OMEGA = zeros(size(OMEGA));
end
% SDB 7/8/26
% Test out vertical motion when horizontal is zeroed out.
if ~vel_nonzero_flag(1)
	u = zeros(size(u));
end
if ~vel_nonzero_flag(2)
	v = zeros(size(v));
end

% Convert omega downward pressure rates to vertical velocities.
if ~plev_flag
	wfull = omega2w(p,numel(t_datenum), numel(lon), numel(lat), T, OMEGA);
end


% Dimensions of u, v, omega are each nlon x nlat x nlev x nt

rel_s = (t_datenum - t_datenum(1))*24*3600;

nt = numel(t_datenum);
nlon = numel(lon);
nlat = numel(lat);
nlev = numel(p);
npts = nlon*nlat*nlev;


% Replace 90 lat with 89.9 to avoid singularity.
lat(find(lat == 90)) = 89.9; % deg
lat(find(lat == -90)) = -89.9; % deg


% Wrap by degrees if 360 degrees are spanned
if max(lon) - min(lon) > 355
	lon_wrapped = [lon; lon(1)+360]; 
	% Repeat the westmost speeds in order to have the domain fully close around
	% in longitude.
	% SDB 6/5/25 Uncommenting this as I test out the 3D code.
	u_wrapped = cat(1,u, u(1,:,:,:));
	v_wrapped = cat(1, v, v(1,:,:,:));
	%omega_wrapped = cat(1, OMEGA, OMEGA(1,:,:));
	% SDB 2/12/26 I said I'd uncommented it in 6/5/25, but it was commented
	% as I checked just now.  Need to work with wfull and wrap it also anyway.
    if plev_flag
    	omega_wrapped = cat(1, OMEGA, OMEGA(1,:,:,:));
    else
    	w_wrapped = cat(1,wfull, wfull(1,:,:,:));
    end
    GPH_wrapped = cat(1, GPH, GPH(1,:,:,:));

    nlon_wrapped = nlon+1;
    npts_wrapped = nlon_wrapped*nlat*nlev;
end

% Convert to angular rates with units of deg/s to match gridpoints units.
% Formulas for rate conversions are in Appendix of Wang et al., 2018.
Re = 6378e3; %m
r = Re + GPH_wrapped; %ht_wrapped*1e3; %m
latmatrad = repmat(lat',nlon_wrapped,1,nlev,nt)*pi/180; % Latitudes in rad.
phidot = u_wrapped./(r.*cos(latmatrad))*180/pi; % angular rate in deg/s.
lambdadot = (v_wrapped./r)*180/pi; % angular rate in deg/s.

% Allocate space for this matrix.
% Each gridpoint has 3 velocity components.
%vel_ts = zeros(nt,1+npts*3);% Non-wrapped version.
vel_ts = zeros(nt, 1+npts_wrapped*3); % Wrapped version.

% List the times down the first column.
vel_ts(:,1) = rel_s;
% Since each time goes on one line, that must be the outer loop.
for tindex = 1:nt
    % Wrapped version.
    vel_ts(tindex,2:3:size(vel_ts,2)) = reshape(phidot(:,:,:,tindex),1,npts_wrapped);
    vel_ts(tindex,3:3:size(vel_ts,2)) = reshape(lambdadot(:,:,:,tindex),1,npts_wrapped);
    if plev_flag
	vel_ts(tindex,4:3:size(vel_ts,2)) = reshape(OMEGA(:,:,:,tindex),1,npts);%_wrapped);
%    vel_ts(tindex,4:3:size(vel_ts,2)) = reshape(omega_wrapped(:,:,:,tindex),1,npts);%_wrapped);
    else
	%vel_ts(tindex,4:3:size(vel_ts,2)) = reshape(wfull(:,:,:,tindex),1,npts);%_wrapped);
	vel_ts(tindex,4:3:size(vel_ts,2)) = reshape(w_wrapped(:,:,:,tindex),1,npts_wrapped);
    end

end

% Restore the original latitudes, including 90 deg.
lat(find(lat == 89.9)) = 90; % deg
lat(find(lat == -89.9)) = -90; % deg

