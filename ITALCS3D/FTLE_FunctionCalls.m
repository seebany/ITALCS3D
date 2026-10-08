% This function computes the FTLE calculation and tracer positions if 
% tracer_flag is set.
%
% Created by Andrew Rukujzo merging work of Ningchao Wang, Zachary Bonson,
% Marta Cabos-Valle.
% Commented and modified by Seebany Datta-Barua
% 8 Sept 2025
% Trying to correct the 3D tracer calculation for WACCMX+DART flows.
% 25 Sept 2025 Integrating over plev can make pressures negative.  Try to work
% with height since we're now passing the whole array anyway. May as well pass
% the altitude grid.

function [FTLE,tracer_pos] = FTLE_FunctionCalls(flow, ...
	xmin,xmax,Xres, ...
	ymin,ymax,Yres, ...%	zmin,zmax,Zres, ...
	zarray, ...
	tmax,TimeStamp,timesteps, ...
	t_initial,IntegrationTime, ...
	CloseDomain_flag,period_flag,computedirection_flag,Euclidean_flag, ...
	dim_active_flag,FTLE_flag, ...
	tracer_flag,tracer_start, ...
	FDFTLE_index,Delta_xyz,fdmin_xyz,fdmax_xyz, ...
	vel_nonzero_flag, J_in_Cartesian_flag)
% Inputs:
%       flow: matrix containing u,v,w velocities at every gridpoint and timestamp.
%       xmin/ymin: minimum x,y distances.
%       xmax/ymax: maximum x,y distances.
%       Xres/Yres: number of gridpoints in each (x,y,z) direction.
% 	zarray: an array of the third coordinate. Passed because with atmospheric model, this is gph and changes at each time.
%       TimeStamp: time difference between data points.
%       timesteps: number of iterations at each timestamp (integration).
%       t_initial: time where the user wants to start the computation.
%       tmax: maximum time the integration can go.
%       IntegrationTime: total time of integration.
%       CloseDomain_flag: defines if particles do(1) or do not(0) wrap around
%       and on the x, y, or z axis(array input: [x,y,z])
%       period_flag: defines if flow is time-periodic or not (1,0).
%       computedirection_flag: forward or backward computation (1,-1).
%       dim_active_flag: a 1x3 array of 0 (do not use) or 1 (use).
%		The sum of the elements should be 2 for 2D or 3 for 3D.
%       FTLE_flag: No FTLE, regular FTLE, HWM14 mean velocity FDFTLE, or 
%       homoclinic linearized velocity FDFTLE (0,1,2,3).
%       tracer_flag: tracer computation on or off (1,0).
%       tracer_start: [n x 2] or [n x 3] array of initial particle positions.
%       FDFTLE_index: finite domain edge smoothing on or off (1,0) for x,y,z 
%       boundaries (array input: [x,y,z])
%       Delta_xyz: distance from each boundary that finite domain will
%       smooth velocity(array input: [x,y,z])
% Outputs:
%       FTLE: FTLE matrix
%       tracer_pos: traced particle(s) position data
%
% Created by Andrew Rukujzo consolidating work of Marta Cabos Valle, Zachary Bonson,
%       Ningchao Wang. 
% Modified by Seebany Datta-Barua
% 16 May 2024 Tracer calculation doesn't pass back multiple initial tracers' positions correctly.
% 23 May 2025 Have to be careful with meshgrids and permuting.  This code creates the spatial grids X, Y, Z using meshgrid. But VelInDirection.m seems to do some permuting.
tracer_pos = [];
x=linspace(xmin,xmax,Xres); 
y=linspace(ymin,ymax,Yres);
%z=linspace(zmin,zmax,Zres);
%z = zarray;
%Zres = numel(z);
Zres = size(zarray, 3);
% SDB 2/12/26 Check if longitudes have been wrapped around, to get dimensions to match.
if xmax - xmin > 357.5
	% Need to augment zarray to repeat the values at the 0 longitude.
	zarray = cat(1,zarray, zarray(1,:,:,:));
end


% 2D Computation
if dim_active_flag == [1,1,0]
    [X,Y]=meshgrid(x,y);
    Z=0;
elseif dim_active_flag == [1, 0, 1]
	    [X, Z] = meshgrid(x, z);
	    Y = 0;
elseif dim_active_flag == [1, 1, 1]
%	% meshgrid will make each of these nlat x nlon x nlev.
%    [X,Y,Z]=meshgrid(x,y,z);
	[X, Y] = meshgrid(x, y);
end
% FD calculations for FTLE_flag > 1, which are finite domain options.
if FTLE_flag > 1
    radis_lat = [];
    re = 6378100;
    ht = 250;
    lat_deg = linspace(ymin,ymax,Yres);
    lat_deg_rad = lat_deg'*pi/180;
    cos_lat_deg = cos(lat_deg_rad);
    for radis_id = 1:Xres,
        radis_lat_line = (ht*1000+re)*cos_lat_deg;
        radis_lat = [radis_lat;radis_lat_line'];
    end
    [t, ~] = size(flow);
    for i = 1:t
        
        u = flow(i,2:3:end-2);
        v = flow(i,3:3:end-1);
        U = reshape(u,Xres,Yres)';
        V = reshape(v,Xres,Yres)';
        
        U = U/180*pi.*radis_lat'; %deg to cart
        V = V*pi/180*(ht*1000+re); %deg to cart
        
        %                 [U, V] = vel_modifiedTopBottomBoundEorW(X,Y,Delta_xyz,fdmin_xyz,fdmax_xyz,U,V);
        
        [U, V] = vel_modified(X,Y,Delta_xyz,fdmin_xyz,fdmax_xyz,U,V,FDFTLE_index,FTLE_flag);
        
        %                 Deltax=Delta_xyz(1);
        %                 Deltay=Delta_xyz(2);
        %                 [U,V] = vel_modifiedv3(X,Y,xmax,xmin,ymax,ymin,Deltax,Deltay,U,V);
        
        U = U*180/pi./radis_lat'; %cart to deg
        V = V*180/pi/(ht*1000+re); %cart to deg
        
        U=U';
        V=V';
        
        u = reshape(U,1,Xres*Yres);
        v = reshape(V,1,Xres*Yres);
        flow(i,2:3:end-2) = u;
        flow(i,3:3:end-1) = v;
    end
    
end % if FTLE_flag > 1 % execute finite-domain options.


% Grid and time spacing. Needed for integration.
dt = TimeStamp/timesteps;
dx = (xmax-xmin)/(Xres-1);
dy = (ymax-ymin)/(Yres-1);
%dz = (zmax-zmin)/(Zres-1);

%% Tracer Calculation
if tracer_flag == 1

    t_end = t_initial+IntegrationTime;
    [num_tracers,~] = size(tracer_start);
    % Option for computing tracers in 2D.
    if dim_active_flag == [1,1,0]
        % Loop over each tracer.
        for i = 1 : num_tracers
            tracer = [];  
            x = tracer_start(i,1);
            y = tracer_start(i,2);
            z = 0;
            tracer = [tracer;x,y];
            % Loop over Delta t.
            for t = t_initial : TimeStamp : t_end
		    [U1,V1,W1,U2,V2,W2]=VelInDirection(flow,Xres,Yres,Zres,t,TimeStamp,computedirection_flag,dim_active_flag);
                % Loop over dt.
		for ind_timesteps = 1:timesteps
                    [x,y,z] = EndPointPosition(t,x,y,z,U1,V1,W1,U2,V2,W2,dt,Xres,Yres,Zres,dx,dy,dz,xmin,ymin,zmin,computedirection_flag,ind_timesteps,TimeStamp,dim_active_flag, spherical_coords_flag);
                end
                tracer = [tracer;x,y];
            end
            % Store in array that is [ntimes x ndims x ntracers].
            tracer_pos(:,:,i) = tracer;
        end
    elseif dim_active_flag == [1,0,1]
        % Loop over each tracer.
        for i = 1 : num_tracers
            tracer = [];  
            x = tracer_start(i,1);
            z = tracer_start(i,3);
            y = 0;
            tracer = [tracer;x,z];
            % Loop over Delta t.
            for t = t_initial : TimeStamp : t_end
                [U1,V1,W1,U2,V2,W2]=VelInDirection(flow,Xres,Yres,Zres,t,TimeStamp,computedirection_flag,dim_active_flag);
                % Loop over dt.
		for ind_timesteps = 1:timesteps
                    [x,y,z] = EndPointPosition(t,x,y,z,U1,V1,W1,U2,V2,W2,dt,Xres,Yres,Zres,dx,dy,dz,xmin,ymin,zmin,computedirection_flag,ind_timesteps,TimeStamp,dim_active_flag, spherical_coords_flag);
                end
                tracer = [tracer;x,z];
            end
            % Store in array that is [ntimes x ndims x ntracers].
            tracer_pos(:,:,i) = tracer;
        end
    % Compute the 3D tracer positions.
    elseif dim_active_flag == [1, 1, 1]

   	for i = 1 : num_tracers
		%disp(['Tracer number ' num2str(i)])
            % Initialize an array of the tracer intermediate positions with the initial position, to be concatenated after the outer integration loop.
	    clear tracer
	    tracer = tracer_start(i,:);

%    	    % Try to use the full height variable.
%	    lonind = find(abs(tracer(1) - x) == min(abs(tracer(1) - x)));
%	    latind = find(abs(tracer(2) - y) == min(abs(tracer(2) - y)));
%
%    	    % If using WACCMX-based model, z is log(plev), and so dz should be approximately
%    	    % regularly spaced around the mesopause altitudes.
%    	    %zrow = find(abs(zarray - z) == min(abs(zarray - z)));
%    	    htind = find(abs(zarray(lonind, latind, :,1) - tracer(3)) == ...
%	    	min(abs(tracer(3) - zarray(lonind, latind, :,1))));
%    	    dz = zarray(lonind, latind, htind + 1,1) ...
%	        - zarray(lonind, latind, htind, 1);


		% Initialize tracer positions to be overwritten at each integration step, both inner and outer.
	 	clear x y z
            	x = tracer(1);
            	y = tracer(2);
            	z = tracer(3);

            for t = t_initial : TimeStamp : t_end-TimeStamp
                % For each integration at data cadence find the flows before and after (or after and before, for backward integration).
		[U1,V1,W1,U2,V2,W2]=VelInDirection(flow,Xres,Yres,Zres,t,TimeStamp,computedirection_flag,dim_active_flag);
% ts: index generated based on the time input
ts = fix(t/TimeStamp) + 1;
	
		% Within the successive data cadence, subdivide the integration.
                for ind_timesteps = 1:timesteps
    			[x,y,z] = EndPointPosition(t,x,y,z,U1,V1,W1,U2,V2,W2,dt,...
		    	Xres,Yres,Zres,dx,dy, ...dz, ...
			xmin, xmax, ymin, ymax, ...zmin, zmax, 
			zarray(:,:,:,ts), ...
			computedirection_flag,ind_timesteps,TimeStamp, CloseDomain_flag, dim_active_flag);
                end
                tracer = [tracer;x,y,z];
            end
            tracer_pos(:,:,i) = tracer;
        end
    end
end % if tracer_flag == 1

% FTLE Calculation
% Initial gridpoints, needed for computing Jacobian.
% 10/20/25 When X, Y get created by FTLE_FunctionCalls.m using meshgrid, they appear to be nlat x nlon. So I need to transpose them before using in EndPosition.m. They don't seem to actually be used in EndPointPosition.m.
% The arrays of X, Y, Z, U, V, W each have dimensions [q x p x l] since rows correspond to latitude and columns to longitude, I think. SDB 1/26/26. No, U, V, W are [q, p, l] because they were created by VelInDirection.m.  We appear to have transposed X, Y, to match Z to be [nlon, nlat, nht] before passing in. So what I need to do is keep X, Y to match U, V, W.  Instead I need to permute Z.
%X = X';
%Y = Y';
%X0 = X;
%Y0 = Y;
%Z0 = zarray(:,:,:,1);
X0 = X;
Y0 = Y;
% Before this step, zarray is [nlon, nlat, nht, nt] = [144, 96, 126, 2]
Z = permute(zarray, [2, 1, 3, 4]);
%Z0 = Z;


if FTLE_flag > 0
    % Forward time
    if computedirection_flag == 1
        t_end = t_initial+IntegrationTime;
        % Time-periodic
	if period_flag == 1
	    % Set limits of integration.
            for t = t_initial : TimeStamp : t_end
                % Modulo to make the flow time-periodic.
		if t >= tmax
                    t = t - tmax;
                end
                % Extract the flow at time 1 and time 2.
                [U1,V1,W1,U2,V2,W2]=VelInDirection(flow,Xres,Yres,Zres,t,TimeStamp,computedirection_flag,dim_active_flag);
		% Integrate over dt's within cadence of the flow fields provided.
		for ind_timesteps = 1:timesteps
                    [X,Y,Z] = EndPosition(t,X,Y,Z,U1,V1,W1,U2,V2,W2,Xres,Yres,Zres,dx,dy,dz,dt,xmin,xmax,ymin,ymax,zmin,zmax,ind_timesteps,TimeStamp,CloseDomain_flag,computedirection_flag,dim_active_flag);
                end
            end % for t = t_initial : TimeStamp : t_end
        else % Non-time-periodic flow.
	    % Convert the 2D horizontal grid into a 3D repeated grid.
	    X3d = repmat(X,1,1,Zres);
	    Y3d = repmat(Y,1,1,Zres);
	    X03d = X3d;
	    Y03d = Y3d;
	    ts = fix(t_initial/TimeStamp) + 1;
	    Z03d = Z(:,:,:,ts); % The point locations Z3D will change over time.
	    %Z03d = Z(:,:,:,1);
	    % 8/20/26 Initialize the z positions that will keep being updated.
	    Z3d = Z03d;

	    % Set limits of integration, in elapsed sec, based on data cadence.
            for t = t_initial : TimeStamp : t_end
                if t >= tmax
                    t = tmax;
                end
		[t, Xres, Yres, Zres]
                % Extract the flow at time 1 and time 2.
		[U1,V1,W1,U2,V2,W2]=VelInDirection(flow,Xres,Yres,Zres,t,TimeStamp,computedirection_flag,dim_active_flag);

		% When the whole vertical matrix is passed in, use it as Z.
		% ts: index generated based on the time input
		ts = fix(t/TimeStamp) + 1;
		% 8/20/26 Commenting out the next line so it doesn't overwrite positions obtained after prior dt integrations, at each TimeStamp.
		%Z3d = Z(:,:,:,ts); % The point locations Z3D will change over time.
		% SDB 8/19/26 Not sure why zarray was reassigned here. May be contributing to errors in FTLE calculation for fig (b) of paper.
		%zarray = Z3d; % The height array is fixed at each time.
		% I changed the input arg of EndPosition() to be zarray instead of passing Z3d twice.
		% Integrate over dt's within cadence of the flow fields provided.
		% Update the positions at each small dt.
                for ind_timesteps = 1:timesteps
			%dz = diff(Z3d,1,3);%zarray(lonind, latind, htind + 1,1) ...
	        %- zarray(lonind, latind, htind, 1);

%    			[x,y,z] = EndPointPosition(t,x,y,z,U1,V1,W1,U2,V2,W2,dt,...
%		    	Xres,Yres,Zres,dx,dy, ...dz, ...
%			xmin, xmax, ymin, ymax, ...zmin, zmax, 
%			zarray(:,:,:,ts), ...
%			computedirection_flag,ind_timesteps,TimeStamp, CloseDomain_flag, dim_active_flag);
		% 8/20/26 Here, args 2-4 are the current positions of tracers at original gridpoints. The input arg after ymax, is the array of heights relevant at this TimeStamp for interpolating speeds.
		[X3d, Y3d, Z3d] = EndPosition(t, X3d, Y3d, Z3d, ...
		    U1,V1,W1,U2,V2,W2,...
		    Xres,Yres,Zres,dx,dy, ...dz,
		    dt,xmin,xmax,ymin,ymax, ...zmin,zmax,
		    Z(:,:,:,ts), ...Z3d, 
		    ind_timesteps,TimeStamp,CloseDomain_flag, ...
		    computedirection_flag,dim_active_flag);
                end %for ind_timesteps = 1:timesteps
            %end %for gridpoint_ind = 1:numel(X3d)
        end %for t = t_initial : TimeStamp : t_end
	%X3d(gridpoint_ind) = x;
	%Y3d(gridpoint_ind) = y;
	%Z(gridpoint_ind) = z;
	end %else % Non-time-periodic flow.
    else % computedirection_flag ~= 1
        t_end = t_initial;
        t_start = t_initial+IntegrationTime;
        if period_flag == 1
            for t = t_start : TimeStamp : t_end
                tc = t_start-t;
                if tc >= tmax
                    t = t + tmax;
                end
                [U1,V1,W1,U2,V2,W2]=VelInDirection(flow,Xres,Yres,Zres,t,TimeStamp,computedirection_flag,dim_active_flag);
                for ind_timesteps = 1:timesteps
                    [X,Y,Z] = EndPosition(t,X,Y,Z,U1,V1,W1,U2,V2,W2,Xres,Yres,Zres,dx,dy,dz,dt,xmin,xmax,ymin,ymax,zmin,zmax,ind_timesteps,TimeStamp ,CloseDomain_flag,computedirection_flag,dim_active_flag);
                end
            end
        else
            for t = t_start : TimeStamp : t_end - TimeStamp
                tc = t_start-t;
                if tc >= tmax
                    t = tmax;
                end
                [U1,V1,W1,U2,V2,W2]=VelInDirection(flow,Xres,Yres,Zres,t,TimeStamp,computedirection_flag,dim_active_flag);
                for ind_timesteps = 1:timesteps
                    [X,Y,Z] = EndPosition(t,X,Y,Z,U1,V1,W1,U2,V2,W2,Xres,Yres,Zres,dx,dy,dz,dt,xmin,xmax,ymin,ymax,zmin,zmax,ind_timesteps,TimeStamp ,CloseDomain_flag,computedirection_flag,dim_active_flag);
                end
            end
        end
    end
Jmatrix = JacobianGradient3datmo(X03d,Y03d,Z03d,X3d,Y3d,Z3d,Xres,Yres,Zres,dim_active_flag,Euclidean_flag);

% 1/30/26 SDB Use a function that transforms the flow map Jacobian from units of lat, lon, ht to m, m, m. Makes dimensions consistent for 3D.
if J_in_Cartesian_flag
	Jtilde = transform_units(X03d, Y03d, Z03d, Jmatrix);

	% Switch to compute FTLE in units of m.
	clear Jmatrix;
	Jmatrix = Jtilde;
end
% 2/20/26 SDB Testing if I can reproduce the 2D FTLE map by selecting the 2x2 submatrix.
% 8/22/26 SDB Having difficulty reproducing the w=0 simulation.  I think seeing the equations I do need to choose the submatrix.
if ~vel_nonzero_flag(3) & ~J_in_Cartesian_flag
temp = Jmatrix;
clear Jmatrix;
Jmatrix = temp(1:2, 1:2, :,:,:);
clear temp
end


Cauchy = zeros(size(Jmatrix));
if dim_active_flag==0
    % Preallocating for speed
    FTLE=zeros(Yres,Xres);
     for i=1:Yres
        for j=1:Xres
            Cauchy(:,:,i,j)=Jmatrix(:,:,i,j)'*Jmatrix(:,:,i,j);
            FTLE(i,j)=(1/abs(IntegrationTime))*(0.5*log(max(eig(Cauchy(:,:,i,j)))));
        end
     end
elseif dim_active_flag == [1, 1, 1]
    % Preallocating for speed
    FTLE=zeros(Yres,Xres,Zres);
    for i = 1:Yres
        for j = 1:Xres
            for k = 1:Zres
                Cauchy(:,:,i,j,k) = Jmatrix(:,:,i,j,k)'*Jmatrix(:,:,i,j,k);
                FTLE(i,j,k) = (1/abs(IntegrationTime))*(0.5*log(max(eig(Cauchy(:,:,i,j,k)))));
            end
        end
    end
end
else
    FTLE=[];
end
end
