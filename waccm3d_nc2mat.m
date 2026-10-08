% function [lon, lat, U, V, GPH, p, t_datenum, OMEGA, T] = waccm3d_nc2mat_pleiades(filename)
% Function to load WACCM nc files and convert to a mat file that can be
% used by ITALCS.
%
% Input string is the full path to the WACCM velocity data .nc file, e.g.,
% /data2/public/Model_outputs/WACCM-X+DART/WACCMX+DART_UV_200901_100km_ensmean_ensstd.nc.
% Input ht is altitude of the atmospheric shell in km.
% Output is a matrix of [ntimes x number gridpoints]. The output is also
% stored as a .mat file for future use.
%
% Outputs have changed as this function has been adapted from other functions.
% For WACCMX+DART 2010 data set,
% u is mean zonal wind in m/s and is nlons x nlats x npressures x ntimes.
%
% S. Datta-Barua
% 19 Dec 2019
% 12 Jan 2021 Return the velocities themselves also u_wrapped, v_wrapped in m/s.
% 19 May 2023 Return the vertical velocities also, for meridional plane study.
% 17 Feb 2026 Adapting to read 6-hr assimilation files.

function [lon, lat, U, V, GPH, p, t_datenum, OMEGA, T] = waccm3d_nc2mat_pleiades(filename, highres_flag)

% Commenting the read data out since we're running a new date, which has
% different file names and variables.
%if ~isempty(strfind(filename{1}, '2010'))
%	%trueht = nan;
%	try
%	% Read longitudes, converting to doubles if they were stored as singles to
%	% make mapping toolbox functions work better.
%	lon = double(ncread(filename{1}, 'LONGITUDE'));
%	lat = double(ncread(filename{1}, 'LATITUDE'));
%	p = double(ncread(filename{1}, 'PRESSURE'));
%	catch
%	% Read longitudes, converting to doubles if they were stored as singles to
%	% make mapping toolbox functions work better.
%	lon = double(ncread(filename{1}, 'lon'));
%	lat = double(ncread(filename{1}, 'lat'));
%	p = double(ncread(filename{1}, 'lev'));
%	%p = nan;
%	GPH = nan;
%	end
%
%	try
%	varstr = {'U','V','GPH'};
%	for i = 1:3
%		% u is mean zonal wind in m/s and is nlons x nlats x npressures x ntimes.
%		cmdstr = [varstr{i} ' = ncread(filename{i}, varstr{i} );'];
%		eval(cmdstr);
%
%	end
%	catch
%		% v = mean meridional wind in m/s, same dimensions as u.
%		U = ncread(filename{1}, 'UMEAN');
%		V = ncread(filename{1}, 'VMEAN');
%	end
%
%	% Read times
%	try
%	t = ncread(filename{i}, 'YYYYMMDDHH'); % numeric YYYYMMDD
%	year = floor(t/1e6);
%	month = floor((t-year*1e6)/1e4);
%	day = floor((t-year*1e6-month*1e4)/1e2);
%	HH = rem(t,1e2);
%	catch
%		yymmdd = ncread(filename{i}, 'date');
%		year = floor(t/1e4);
%	month = floor((t-year*1e4)/1e2);
%	day = floor((t-year*1e4-month*1e2)/1e0);
%
%	HH = double(ncread(filename{j}, 'datesec'))/3600;
%	end
%	t_datenum = datenum(double([year, month, day, HH, ...
%		zeros(numel(year),2)]));
%%return
%%else

% SDB 2/17/26 The 6-hr assimilations have an array of datesec, so initialize
% an empty array.
% SDB 3/10/26 It seems the hourly assimilations need this initialized also.
t_datenum = [];

%	% For working with Nick's waccmx-dart runs directly.
	for j = 1:numel(filename)
		% Read longitudes, converting to doubles if they were stored as singles to
		% make mapping toolbox functions work better.
%		try
		if j == 1
			lon = double(ncread(filename{j}, 'lon'));%'LONGITUDE'));
			lat = double(ncread(filename{j}, 'lat'));%'LATITUDE'));
			%trueht = nan;
			p = double(ncread(filename{j}, 'lev'));%'PRESSURE'));
%		catch
%			lon = double(ncread(filename{j}, 'LONGITUDE'));
%			lat = double(ncread(filename{j}, 'LATITUDE'));
			%trueht = nan;
			%p(:,j) = double(ncread(filename{j}, 'PRESSURE'));
		end
		% Read times
		try
		    % SDB 2/16/26 This reads the 6-hr assimilation files, where
		    % sec is a 6x1 array.
		    t = ncread(filename{j}, 'date'); % numeric YYYYMMDD
		    year = floor(t/1e4);
		    month = floor((t-year*1e4)/1e2);
		    %day = floor((t-year*1e4-month*1e2);%/1e2);
		    day = rem(t,1e2);
		    sec = double(ncread(filename{j}, 'datesec'));
		    %t_datenum(j) = datenum(double([year, month, day, 0, 0, sec]));
		    % Append timestamps
		    t_datenum = [t_datenum; 
		    datenum(double([year, month, day, ...
			zeros(numel(sec),2), sec]))];
		catch
		    t = ncread(filename{j}, 'YYYYMMDDHH');
		    year = floor(t/1e6);
		    month = floor((t-year*1e6)/1e4);
		    day = floor((t-year*1e4-month*1e2));%/1e2);
		    hour = rem(t,1e2);
		    t_datenum(j) = datenum(double([year, month, day, hour, 0, 0]));
		end

		varstr = {'U','V','Z3', 'OMEGA', 'T'};%'GPH'};
		
		for i = 1:numel(varstr)
		    % SDB 2/17/26 In the hourly assimilation files,
		    % u is mean zonal wind in m/s and is nlons x nlats x npressures x ntimes.
		    try
		    	cmdstr = [varstr{i} '(:,:,:,j) = double(ncread(filename{j}, varstr{i} ));'];
			
			eval(cmdstr);

		    % u is mean zonal wind in m/s and is nlons x nlats x npressures x (assim_rate/model_rate) x ntimes (of which there may be multiple in a single file).
		%cmdstr = [varstr{i} ' = [' varstr{i} '; ncread(fullname, varstr{i} )];'];
		    catch
		    	try
			cmdstr = [varstr{i} '5d(:,:,:,:,j) = ncread(filename{j}, varstr{i} );'];
			eval(cmdstr);
		    	catch
			cmdstr = [varstr{i} '5d(:,:,:,:,j) = nan'];
		    	end
		    end
		end % for i
	end % for j
%end % ~isempty(strfind(filename{1}, '2010'))

% In Nick's file the geopotential height is variable Z3.
if highres_flag == 2 %exist('Z3','var')
	GPH = Z3;
%	clear Z3
elseif highres_flag == 0 %exist('Z35d', 'var')
    GPH = reshape(Z35d,numel(lon),numel(lat),size(Z35d,3),[]);
% In Nick's 6-hrly assim file the geopotential height is variable Z3.
U = reshape(U5d, numel(lon), numel(lat), size(U5d,3),[]);
V = reshape(V5d, numel(lon), numel(lat), size(V5d,3),[]);
OMEGA = reshape(OMEGA5d, numel(lon), numel(lat), size(OMEGA5d,3),[]);
T = reshape(T5d, numel(lon), numel(lat), size(T5d,3),[]);
%clear U5d, V5d, Z35d

end

return

