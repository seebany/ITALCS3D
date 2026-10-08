% EndPointPostion: This function computes x,y,z final position of a single
% point at a given time-step estimated based on the velocity of the point and
% surrounding points. In FTLE_FunctionCalls, this function is run at time 
% intervals (t_loc) between the given time-step and the next, in order to
% increase the accuracy of the final postion values.

function [xf,yf,zf] = EndPointPosition(t,xi,yi,zi, ...
	U1,V1,W1,U2,V2,W2,dt, ...
	Xres,Yres,Zres, ...
	dx,dy,... dz, ...
	xmin,xmax, ...
	ymin,ymax, ...zmin,zmax, 
	zarray, ...
	computedirection_flag,ind_timesteps,TimeStamp,CloseDomain_flag, dimension_flag)
% Inputs:
%       t: time within the integration slot(not an array its a specific time)
%       x/y/z: initial x,y,z positions of tracked particle.
%       U1/V1/W1: x,y,z-direction velocity field at time t (calculated 
%       using VelInDirection)
%       U2/V2/W2: x,y,z-direction velocity field at time t + TimeStamp
%       (calculated using VelInDirection).
%       Xres/Yres/Zres: number of gridpoints in each (x,y,z) direction.
%       dx/dy/dz: distance between consecutive gridpoints.
%       dt: time difference between timesteps at each integration.
%       xmin/ymin/zmin: minimum x,y,z distances.
%       xmax/ymax/zmax: maximum x,y,z values.
%       ind_timesteps: timestep in the inner integration.
%       TimeStamp: time difference between data points.
%       CloseDomain_flag: defines if particles do(1) or do not(0) wrap around
%       and on the x, y, or z axis(array input=[x,y,z])
%       computedirection_flag: forward or backward computation (1,-1).
%       dimension_flag: a 1x3 [x, y, z] truth array depending on whether that dimension has active motion (1) or not (0).
% Outputs:
%       x/y/z: initial, intermediate, and final x,y,z positions after integration,
%		each variable arrayed as [ntimes x ?? x ntracers?]
% Modified by Seebany Datta-Barua
% 7 Oct 2025
%
% Passing the height array and using 3-d atmosphere. 
% Note that x is longitude, y is latitude, z is height.
% Indices used internally are p for longitude, q for latitude, l for height.

% 9/9/25 For 3D with the 3rd dimension altitude, trying to pass zarray, so need zmin and Zres extracted from it.
zmin = min(zarray);
Zres = size(zarray,3);

%-------------Step 1: linear time-interpolation of the velocity field.
% t_loc: time interval between time-steps based on the variable "timesteps"
% Example: if timesteps=10, then t_loc=0.1,0.2,0.2,...,0.9 so between t=1
% and t=2, we get t=1.1,1.2,1.3,...,1.9
t_loc = (mod(t,TimeStamp) + ind_timesteps*dt)/TimeStamp;

% p: column index generated based on a given X coordinate
% q: row index generated based on a given Y coordinate
% l: third index generated based on a given coordinate
% Example: if minimum x value is 5, then x=5 returns a p value of 1
%p = fix((x-xmin)/dx)+1;
%q = fix((y-ymin)/dy)+1;
%l = fix((z-zmin)/dz)+1;

% U,V,W: Velocities estimated based on velocities at 2 consecutive
% time-steps and the t_loc between those time-steps.
% Example: if U=0 at t=1, U=1 at t=2 and t_loc=0.5 (halfway between
% timesteps) U=0.5 (halfway between the velocities)
U = U2*t_loc + U1*(1-t_loc);
V = V2*t_loc + V1*(1-t_loc);
W = W2*t_loc + W1*(1-t_loc);


% 2D Computation along the x, y dimensions with motion u, v.
if dimension_flag == [1,1,0]
    if (p<=0)||(p>Xres) || (q<=0) || (q>Yres)
        % Checks if the point is outside of the domain otherwise it adds a velocity to the point
        x = x;
        y = y;
    else
        if (p<Xres)
            % x_loc: weighting factor for points between gridpoints
            % Example: if grid points are x=7.5 and x=7.4 then x=7.45
            % returns x_loc=0.5
            x_loc = (x-xmin-(p-1)*dx)/dx;
        else
            x_loc = 1;
        end
        
        if (q<Yres)
            % y_loc: weighting factor for points between gridpoints
            y_loc = (y-ymin-(q-1)*dy)/dy;
        else
            y_loc = 1;
        end
        
        v1 = x_loc * y_loc;
        v2 = x_loc * (1-y_loc);
        v3 = (1-x_loc) * y_loc;
        v4 = (1-x_loc) * (1-y_loc);
        
        if (p<Xres)&&(q<Yres)
            % f1x-f4x/f1y-f4y: the chosen point and neighboring points based on
            % position in the grid(Default: the point itself, the point to the
            % right of it, the point beneath it, and the point right one
            % and down one)
            % Example: If the given point is the top right corner of the grid, 
            % then there are no points to the right of it so only the given
            % point, and the point beneath it are used for interpolation.
            f1x = U(q,p);
            f2x = U(q+1,p);
            f3x = U(q,p+1);
            f4x = U(q+1,p+1);
            
            f1y = V(q,p);
            f2y = V(q+1,p);
            f3y = V(q,p+1);
            f4y = V(q+1,p+1);

         else if (p<Xres)&&(q == Yres)
                  f1x = U(q,p);
                  f2x = U(q,p);
                  f3x = U(q,p+1);
                  f4x = U(q,p+1);
                  
                  f1y = V(q,p);
                  f2y = V(q,p);
                  f3y = V(q,p+1);
                  f4y = V(q,p+1);
                  else if (p == Xres)&&(q<Yres)
                           f1x = U(q,p);
                           f2x = U(q+1,p);
                           f3x = U(q,p);
                           f4x = U(q+1,p);
                           
                           f1y = V(q,p);
                           f2y = V(q+1,p);
                           f3y = V(q,p);
                           f4y = V(q+1,p);
                       else
                           f1x = U(q,p);
                           f2x = U(q,p);
                           f3x = U(q,p);
                           f4x = U(q,p);
                           
                           f1y = V(q,p);
                           f2y = V(q,p);
                           f3y = V(q,p);
                           f4y = V(q,p);
                       end
              end
        end
        dxdt = f1x*v4 + f2x*v3 + f3x*v2 + f4x*v1;
        dydt = f1y*v4 + f2y*v3 + f3y*v2 + f4y*v1;
        
        x = x + dxdt*dt*computedirection_flag;
        y = y + dydt*dt*computedirection_flag;
    end
% 2D Computation along the x, z dimensions with motion u, w. Untested.
elseif dimension_flag == [1,0,1]
	% To make this work, I need to change all q's to l's, y->z, v->w. SDB 7/12/24.
    if (p<=0)||(p>Xres) || (q<=0) || (q>Yres)
        % Checks if the point is outside of the domain otherwise it adds a velocity to the point
        x = x;
        y = y;
    else
        if (p<Xres)
            % x_loc: weighting factor for points between gridpoints
            % Example: if grid points are x=7.5 and x=7.4 then x=7.45
            % returns x_loc=0.5
            x_loc = (x-xmin-(p-1)*dx)/dx;
        else
            x_loc = 1;
        end
        
        if (q<Yres)
            % y_loc: weighting factor for points between gridpoints
            y_loc = (y-ymin-(q-1)*dy)/dy;
        else
            y_loc = 1;
        end
        
        v1 = x_loc * y_loc;
        v2 = x_loc * (1-y_loc);
        v3 = (1-x_loc) * y_loc;
        v4 = (1-x_loc) * (1-y_loc);
        
        if (p<Xres)&&(q<Yres)
            % f1x-f4x/f1y-f4y: the chosen point and neighboring points based on
            % position in the grid(Default: the point itself, the point to the
            % right of it, the point beneath it, and the point right one
            % and down one)
            % Example: If the given point is the top right corner of the grid, 
            % then there are no points to the right of it so only the given
            % point, and the point beneath it are used for interpolation.
            f1x = U(q,p);
            f2x = U(q+1,p);
            f3x = U(q,p+1);
            f4x = U(q+1,p+1);
            
            f1y = V(q,p);
            f2y = V(q+1,p);
            f3y = V(q,p+1);
            f4y = V(q+1,p+1);

         elseif (p<Xres)&&(q == Yres)
              f1x = U(q,p);
              f2x = U(q,p);
              f3x = U(q,p+1);
              f4x = U(q,p+1);
                  
              f1y = V(q,p);
              f2y = V(q,p);
              f3y = V(q,p+1);
              f4y = V(q,p+1);
         elseif (p == Xres)&&(q<Yres)
              f1x = U(q,p);
              f2x = U(q+1,p);
              f3x = U(q,p);
              f4x = U(q+1,p);
                           
              f1y = V(q,p);
              f2y = V(q+1,p);
              f3y = V(q,p);
              f4y = V(q+1,p);
         else
              f1x = U(q,p);
              f2x = U(q,p);
              f3x = U(q,p);
              f4x = U(q,p);
                           
              f1y = V(q,p);
              f2y = V(q,p);
              f3y = V(q,p);
              f4y = V(q,p);
%                       end
%              end
        end
        dxdt = f1x*v4 + f2x*v3 + f3x*v2 + f4x*v1;
        dydt = f1y*v4 + f2y*v3 + f3y*v2 + f4y*v1;
        
        x = x + dxdt*dt*computedirection_flag;
        y = y + dydt*dt*computedirection_flag;
    end
%------------------- 3D Computation
elseif dimension_flag == [1, 1, 1]
% Check if the point needs to be wrapped around the domain before integration, or is outside of the domain otherwise it adds a velocity to the point
%if (p<=0) || (p>Xres) || (q<=0) || (q>Yres) || (l<=0) || (l>Zres) 
% Copying the 30 lines from EndPosition.m to wrap points around correctly.

    % Particles wrap around in the x (longitude if CloseDomain_flag(1)==1 ) direction
    % SDB 6/6/25 This seems to have to come after if CloseDomain_flag(2) otherwise the indices end up negative.  Not sure why it works to have this first in EndPosition.m.
    if CloseDomain_flag(1) == 1
	if xi < xmin
             xi = xi + dx*(Xres-1)*(1+abs(fix((xi - xmin)/((Xres-1)*dx)))); %X(i,j) could be much smaller than xmin
        elseif xi > xmax
             xi = xi - dx*(Xres-1)*(1+abs(fix((xi - xmax)/((Xres-1)*dx)))); %X(i,j) could be much bigger than xmax
%        else
%             x = x;
        end
    end
    % Particles wrap around in the y (latitude if CloseDomain_flag(2) == 1) direction
    if CloseDomain_flag(2) == 1
	% Check if longitude is negative and select the power of -1 to be used.
%	if x < xmin || x > xmax
            sign_id = (xi < xmin); %1;
%        else
%            sign_id = 0;
%        end
%             x = x - (xmax-xmin)/2*(-1)^sign_id; %v2 updated here for the latitude boundary
%	end
	% If below -90 deg lat in the southern hemisphere. 
        if yi < ymin
             % Bring to the opposite hemisphere by shifting 180 deg
%		x = x - (xmax-xmin)/2*(-1)^sign_id; %v2 updated here for the latitude boundary
		% Try an alternative function to shift by 180 deg.
		xi = mod(xi, 180) + 180*(xi == mod(xi, 180));
		yi = 2*ymin - yi; % if Y(i,j) is much smallerthan ymin, need to fix!!!!!! but for this flow, need not right now
        elseif yi > ymax
             % Bring to the opposite hemisphere by shifting 180 deg
%             x = x - (xmax-xmin)/2*(-1)^sign_id; %v2 updated here for the latitude boundary;
		% Try an alternative function to shift by 180 deg.
		xi = mod(xi, 180) + 180*(xi == mod(xi, 180));
             yi = 2 * ymax - yi;
%        else
%             x = x;
%             y = y;
        end
%                    end
    end
%        x = x;
%        y = y;
%        z = z;
%else
%                 end

% Now that the positions have been wrapped, recompute the indices.
p = fix((xi - xmin)/dx)+1;
q = fix((yi - ymin)/dy)+1;
% Not enough to find nearest.  Need to find nearest lower, to recreate fix.
l = find(abs(zi - zarray(p,q,:)) == min(abs(zi - zarray(p,q,:))));%53;%fix((z-zmin)/dz)+1;
%disp([q p l zarray(p,q,l) zarray(p,q,l+1)])

%l = max(find(zi - zarray(q,p,:) >= 0));%53;%fix((z-zmin)/dz)+1;
    % Checks if the point is outside of the domain otherwise it adds a velocity to the point
if isempty(l)
	xf = xi;
	yf = yi;
	zf = zi;
	return
end
if (p<=0) || (p>Xres) || (q<=0) || (q>Yres) || (l<=0) || (l>Zres)
	xf = xi;
	yf = yi;
	zf = zi;
else
    % Set the fractions for tri-linear spatial interpolation.
    if p>=Xres
         x_loc = 1;
    else
        % x_loc: weighting factor for points between gridpoints
        x_loc = (xi - xmin-(p-1)*dx)/dx;
    end
    if q>=Yres
        y_loc = 1;
    else
        % y_loc: weighting factor for points between gridpoints
        y_loc = (yi - ymin-(q-1)*dy)/dy;
    end
    if l>=Zres
        z_loc = 1;
    else
        % z_loc: weighting factor for points between gridpoints
        %z_loc = (z-zmin-(l-1)*dz)/dz;
	% Modify for passing the whole irregular zarray.
	%z_loc = (zi - zarray(l))/(zarray(l+1) - zarray(l));
	% Modify for passing the whole 3d grid.
	z_loc = (zi - zarray(p,q,l))/(zarray(p,q, l+1) - zarray(l));
    end

    v1 = (1-x_loc)*(1-y_loc)*(1-z_loc);
    v2 = x_loc*(1-y_loc)*(1-z_loc);
    v3 = (1-x_loc)*y_loc*(1-z_loc);
    v4 = x_loc*y_loc*(1-z_loc);
    v5 = (1-x_loc)*(1-y_loc)*z_loc;
    v6 = x_loc*(1-y_loc)*z_loc;
    v7 = (1-x_loc)*y_loc*z_loc;
    v8 = x_loc*y_loc*z_loc;

    if p<Xres && q<Yres && l<Zres            
        
        %A111_x-A222_x/A111_7-A222_y/A111_z-A222_z: the chosen point and
        % neighboring points based on its position in the grid (for 3D these 
        % are the point itself and 7 other points which are all 
        % combinations of adding 1 to the x, y, and z position)
        % Example: If the point is (0,0,0), then the neighboring points are 
        % (0,0,1), (0,1,0), (0,1,1), (1,0,0), (1,0,1),(1,1,0), and (1,1,1).  
        % Exemptions are made the same way they are in the 2D case (if
        % neighboring points are off the grid, they are not used for the
        % interpolation).
	% The subscripts on the "A" variables are 1 plus the zero and one indices
	% listed above. So point (0,0,0) is for point (q, p, l) and A111. SDB 6/6/25
	A111_x = U(q,p,l);
        A112_x = U(q,p,l+1);
        A211_x = U(q+1,p,l);
        A212_x = U(q+1,p,l+1);
        A121_x = U(q,p+1,l);
        A122_x = U(q,p+1,l+1);
        A221_x = U(q+1,p+1,l);
        A222_x = U(q+1,p+1,l+1);

        A111_y = V(q,p,l);
        A112_y = V(q,p,l+1);
        A211_y = V(q+1,p,l);
        A212_y = V(q+1,p,l+1);
        A121_y = V(q,p+1,l);
        A122_y = V(q,p+1,l+1);
        A221_y = V(q+1,p+1,l);
        A222_y = V(q+1,p+1,l+1);

        A111_z = W(q,p,l);
        A112_z = W(q,p,l+1);
        A211_z = W(q+1,p,l);
        A212_z = W(q+1,p,l+1);
        A121_z = W(q,p+1,l);
        A122_z = W(q,p+1,l+1);
        A221_z = W(q+1,p+1,l);
        A222_z = W(q+1,p+1,l+1);
    
    elseif p<Xres && q<Yres && l==Zres
        A111_x = U(q,p,l);
        A112_x = U(q,p,l);
        A211_x = U(q+1,p,l);
        A212_x = U(q+1,p,l);
        A121_x = U(q,p+1,l);
        A122_x = U(q,p+1,l);
        A221_x = U(q+1,p+1,l);
        A222_x = U(q+1,p+1,l);

        A111_y = V(q,p,l);
        A112_y = V(q,p,l);
        A211_y = V(q+1,p,l);
        A212_y = V(q+1,p,l);
        A121_y = V(q,p+1,l);
        A122_y = V(q,p+1,l);
        A221_y = V(q+1,p+1,l);
        A222_y = V(q+1,p+1,l);

        A111_z = W(q,p,l);
        A112_z = W(q,p,l);
        A211_z = W(q+1,p,l);
        A212_z = W(q+1,p,l);
        A121_z = W(q,p+1,l);
        A122_z = W(q,p+1,l);
        A221_z = W(q+1,p+1,l);
        A222_z = W(q+1,p+1,l);

    elseif p<Xres && q==Yres && l<Zres
        A111_x = U(q,p,l);
        A112_x = U(q,p,l+1);
        A211_x = U(q,p,l);
        A212_x = U(q,p,l+1);
        A121_x = U(q,p+1,l);
        A122_x = U(q,p+1,l+1);
        A221_x = U(q,p+1,l);
        A222_x = U(q,p+1,l+1);

        A111_y = V(q,p,l);
        A112_y = V(q,p,l+1);
        A211_y = V(q,p,l);
        A212_y = V(q,p,l+1);
        A121_y = V(q,p+1,l);
        A122_y = V(q,p+1,l+1);
        A221_y = V(q,p+1,l);
        A222_y = V(q,p+1,l+1);

        A111_z = W(q,p,l);
        A112_z = W(q,p,l+1);
        A211_z = W(q,p,l);
        A212_z = W(q,p,l+1);
        A121_z = W(q,p+1,l);
        A122_z = W(q,p+1,l+1);
        A221_z = W(q,p+1,l);
        A222_z = W(q,p+1,l+1);

    elseif p==Xres && q<Yres && l<Zres
        A111_x = U(q,p,l);
        A112_x = U(q,p,l+1);
        A211_x = U(q+1,p,l);
        A212_x = U(q+1,p,l+1);
        A121_x = U(q,p,l);
        A122_x = U(q,p,l+1);
        A221_x = U(q+1,p,l);
        A222_x = U(q+1,p,l+1);

        A111_y = V(q,p,l);
        A112_y = V(q,p,l+1);
        A211_y = V(q+1,p,l);
        A212_y = V(q+1,p,l+1);
        A121_y = V(q,p,l);
        A122_y = V(q,p,l+1);
        A221_y = V(q+1,p,l);
        A222_y = V(q+1,p,l+1);

        A111_z = W(q,p,l);
        A112_z = W(q,p,l+1);
        A211_z = W(q+1,p,l);
        A212_z = W(q+1,p,l+1);
        A121_z = W(q,p,l);
        A122_z = W(q,p,l+1);
        A221_z = W(q+1,p,l);
        A222_z = W(q+1,p,l+1);

    elseif p<Xres && q==Yres && l==Zres
        A111_x = U(q,p,l);
        A112_x = U(q,p,l);
        A211_x = U(q,p,l);
        A212_x = U(q,p,l);
        A121_x = U(q,p+1,l);
        A122_x = U(q,p+1,l);
        A221_x = U(q,p+1,l);
        A222_x = U(q,p+1,l);

        A111_y = V(q,p,l);
        A112_y = V(q,p,l);
        A211_y = V(q,p,l);
        A212_y = V(q,p,l);
        A121_y = V(q,p+1,l);
        A122_y = V(q,p+1,l);
        A221_y = V(q,p+1,l);
        A222_y = V(q,p+1,l);

        A111_z = W(q,p,l);
        A112_z = W(q,p,l);
        A211_z = W(q,p,l);
        A212_z = W(q,p,l);
        A121_z = W(q,p+1,l);
        A122_z = W(q,p+1,l);
        A221_z = W(q,p+1,l);
        A222_z = W(q,p+1,l);

    elseif p==Xres && q<Yres && l==Zres
        A111_x = U(q,p,l);
        A112_x = U(q,p,l);
        A211_x = U(q+1,p,l);
        A212_x = U(q+1,p,l);
        A121_x = U(q,p,l);
        A122_x = U(q,p,l);
        A221_x = U(q+1,p,l);
        A222_x = U(q+1,p,l);

        A111_y = V(q,p,l);
        A112_y = V(q,p,l);
        A211_y = V(q+1,p,l);
        A212_y = V(q+1,p,l);
        A121_y = V(q,p,l);
        A122_y = V(q,p,l);
        A221_y = V(q+1,p,l);
        A222_y = V(q+1,p,l);

        A111_z = W(q,p,l);
        A112_z = W(q,p,l);
        A211_z = W(q+1,p,l);
        A212_z = W(q+1,p,l);
        A121_z = W(q,p,l);
        A122_z = W(q,p,l);
        A221_z = W(q+1,p,l);
        A222_z = W(q+1,p,l);


    elseif p==Xres && q==Yres && l<Zres
        A111_x = U(q,p,l);
        A112_x = U(q,p,l+1);
        A211_x = U(q,p,l);
        A212_x = U(q,p,l+1);
        A121_x = U(q,p,l);
        A122_x = U(q,p,l+1);
        A221_x = U(q,p,l);
        A222_x = U(q,p,l+1);

        A111_y = V(q,p,l);
        A112_y = V(q,p,l+1);
        A211_y = V(q,p,l);
        A212_y = V(q,p,l+1);
        A121_y = V(q,p,l);
        A122_y = V(q,p,l+1);
        A221_y = V(q,p,l);
        A222_y = V(q,p,l+1);

        A111_z = W(q,p,l);
        A112_z = W(q,p,l+1);
        A211_z = W(q,p,l);
        A212_z = W(q,p,l+1);
        A121_z = W(q,p,l);
        A122_z = W(q,p,l+1);
        A221_z = W(q,p,l);
        A222_z = W(q,p,l+1);

    else
        A111_x = U(q,p,l);
        A112_x = U(q,p,l);
        A211_x = U(q,p,l);
        A212_x = U(q,p,l);
        A121_x = U(q,p,l);
        A122_x = U(q,p,l);
        A221_x = U(q,p,l);
        A222_x = U(q,p,l);

        A111_y = V(q,p,l);
        A112_y = V(q,p,l);
        A211_y = V(q,p,l);
        A212_y = V(q,p,l);
        A121_y = V(q,p,l);
        A122_y = V(q,p,l);
        A221_y = V(q,p,l);
        A222_y = V(q,p,l);

        A111_z = W(q,p,l);
        A112_z = W(q,p,l);
        A211_z = W(q,p,l);
        A212_z = W(q,p,l);
        A121_z = W(q,p,l);
        A122_z = W(q,p,l);
        A221_z = W(q,p,l);
        A222_z = W(q,p,l);
    end

                        dxdt = A111_x*v1 + A121_x*v2 + A211_x*v3 + A221_x*v4 + A112_x*v5 + A122_x*v6 + A212_x*v7 + A222_x*v8;
                        dydt = A111_y*v1 + A121_y*v2 + A211_y*v3 + A221_y*v4 + A112_y*v5 + A122_y*v6 + A212_y*v7 + A222_y*v8;
                        dzdt = A111_z*v1 + A121_z*v2 + A211_z*v3 + A221_z*v4 + A112_z*v5 + A122_z*v6 + A212_z*v7 + A222_z*v8;
% SDB 6/6/25 These were the ones originally in this code, but commented out in the equivalent spot in EndPosition.m.  I think these might be incorrect and have pasted above the ones that would make the tracers be computed consistently with the FTLEs.
%    dxdt = A111_x*v1 + A211_x*v2 + A121_x*v3 + A221_x*v4 + A112_x*v5 * A212_x*v6 + A122_x*v7 + A222_x*v8;
%    dydt = A111_y*v1 + A211_y*v2 + A121_y*v3 + A221_y*v4 + A112_y*v5 * A212_y*v6 + A122_y*v7 + A222_y*v8;
%    dzdt = A111_z*v1 + A211_z*v2 + A121_z*v3 + A221_z*v4 + A112_z*v5 * A212_z*v6 + A122_z*v7 + A222_z*v8;
    xf = xi + dxdt*dt*computedirection_flag;
    yf = yi + dydt*dt*computedirection_flag;
    zf = zi + dzdt*dt*computedirection_flag;
    %disp([dzdt W(q,p,l)])
	       end	       

end
%end
