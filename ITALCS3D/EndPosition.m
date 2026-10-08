% End Postion: This function computes x,y,z final positions at each gridpoint 
% at a given timestep estimated based on the velocity of the point and
% surrounding points. In FTLE_FunctionCalls, this function is run at time 
% intervals (t_loc) between the given time-step and the next, in order to
% increase the accuracy of the final postion values.

function [X_end,Y_end,Z_end] = EndPosition(t,X,Y,Z, ...
	U1,V1,W1,U2,V2,W2,Xres,Yres,Zres,dx,dy, ...dz,
	dt,xmin,xmax,ymin,ymax, ...zmin,zmax, ...
	zarray, ...
	ind_timesteps,TimeStamp,CloseDomain_flag,computedirection_flag,dimension_flag)

% Inputs:
%       t: time within the integration slot(not an array its a specific time)
%       X/Y/Z: initial x,y,z positions(the grid).
%       U1/V1/W1: x,y,z-direction velocity field at time t (calculated
%       using VelInDirection)
%       U2/V2/W2: x,y,z-direction velocity field at time t + TimeStamp
%       (calculated using VelInDirection).
%       Xres/Yres/Zres: number of gridpoints in each (x,y,z) direction.
%       dx/dy/dz: distance between consecutive gridpoints.
%       dt: time difference between timesteps at each integration.
%       xmin/ymin/zmin: minimum x,y,z values.
%       xmax/ymax/zmax: maximum x,y,z values.
%       ind_timesteps: timestep within the integration time where you are at.
%       TimeStamp: time difference between data points.
%       CloseDomain_flag: defines if particles do(1) or do not(0) wrap around
%       and on the x, y, or z axis(array input=[x,y,z])
%       computedirection_flag: forward or backward computation (1,-1).
%       dimension_flag: a 1x3 [x, y, z] truth array depending on whether that dimension has active motion (1) or not (0).
% Outputs:
%       Xend/Yend/Zend: final x,y,z positions after integration
% Modified by Seebany Datta-Barua
% 20 Oct 2025
%
% Passing the height array and using 3-d atmosphere. 
% Note that x is longitude, y is latitude, z is height.
% Indices used internally are p for longitude, q for latitude, l for height.
% The arrays of X, Y, Z, U, V, W each have dimensions [q x p x l] since rows correspond to latitude and columns to longitude, I think. SDB 1/26/26. No, U, V, W are [q, p, l] because they were created by VelInDirection.m.  We appear to have transposed X, Y, to match Z to be [nlon, nlat, nht] before passing in. So what I need to do is keep X, Y to match U, V, W.  Instead I need to permute Z.

% 10/20/25 For 3D with the 3rd dimension altitude, use initial points Z as zarray, so need zmin and Zres extracted from it. Need to keep Z (updated gridpoints) separate from zarray (static 3D grid).
% 08/19/26 Since zarray is now passed as an input arg, don't need to assign it here.
%zarray = Z;
zmin = min(Z);
Zres = size(Z,3);
%------------Step 1: linear time-interpolation of the velocity field.
% t_loc: time interval between time-steps based on the variable "timesteps"
% Example: if timesteps=10, then t_loc=0.1,0.2,0.3,...,0.9 so between t=1
% and t=2, we get t=1.1,1.2,1.3,...,1.9
t_loc=(mod(t,TimeStamp)+ind_timesteps*dt)/TimeStamp;

% U,V,W: Velocities estimated based on velocities at 2 consecutive
% time-steps and the t_loc between those time-steps.
% Example: if U=0 at t=1, U=1 at t=2 and t_loc=0.5 (halfway between
% timesteps) U=0.5 (halfway between the velocities)
U = U2*t_loc + U1*(1-t_loc);
V = V2*t_loc + V1*(1-t_loc);
W = W2*t_loc + W1*(1-t_loc);
% 2D Calculation
if any(dimension_flag == 0)
    % Preallocating for speed
    X_end = zeros(Yres,Xres);
    Y_end = zeros(Yres,Xres);
    Z_end = [];
    for i=1:Yres
        for j=1:Xres
            if isnan(X(i,j))==1 || isnan(Y(i,j))==1
                % Checks if coordinates of particle velocity are NaN otherwise it adds a velocity to the coordinate
                X_end(i,j) = X(i,j);
                Y_end(i,j) = Y(i,j);
            else
                % Particles wrap around in the x direction
                if CloseDomain_flag(1) == 1
                    if X(i,j) < xmin
                        X(i,j) = X(i,j)+dx*(Xres-1)*(1+abs(fix((X(i,j)-xmin)/((Xres-1)*dx)))); %X(i,j) could be much smaller than xmin
                    elseif X(i,j) > xmax
                        X(i,j) = X(i,j)-dx*(Xres-1)*(1+abs(fix((X(i,j)-xmax)/((Xres-1)*dx)))); %X(i,j) could be much bigger than xmax
%                    else
%                        X(i,j) = X(i,j);
%                        end
                    end
                 end
                % Particles wrap around in the y driection
                if CloseDomain_flag(2) == 1
                    if X(i,j) < 0
                        sign_id = 1;
                     else
                        sign_id = 0;
                     end
                    if Y(i,j)<ymin
                        X(i,j) = X(i,j)-(xmax-xmin)/2*(-1)^sign_id; %v2 updated here for the latitude boundary
                        Y(i,j) = 2*ymin-Y(i,j); % if Y(i,j) is much smallerthan ymin, need to fix!!!!!! but for this flow, need not right now
                    else if Y(i,j)>ymax 
                        X(i,j) = X(i,j)-(xmax-xmin)/2*(-1)^sign_id; %v2 updated here for the latitude boundary;
                        Y(i,j) = 2*ymax-Y(i,j);
                    else
                        X(i,j) = X(i,j);
                        Y(i,j) = Y(i,j);
                        end
                    end
                end
                
                % p: column index generated based on a given X coordinate
                % q: row index generated based on a given Y coordinate
                % Example: if minimum x value is 5, then x=5 returns a p value of 1
                p = fix((X(i,j)-xmin)/dx)+1;
                q = fix((Y(i,j)-ymin)/dy)+1;
            
                if (p<=0) || (p>Xres) ||(q<=0) || (q>Yres)
                    % Checks if coordinates are outside of the domain otherwise it adds a velocity to the coordinate
                    X_end(i,j) = X(i,j);
                    Y_end(i,j) = Y(i,j);
                else
                    if (p<Xres)
                        % X_loc: weighting factor for points between gridpoints
                        % Example: if grid points are x=7.5 and x=7.4 then x=7.45
                        % returns x_loc=0.5
                        X_loc = (X(i,j)-xmin-(p-1)*dx)/dx;
                    else
                        X_loc = 1;
                    end
                        
                    if (q<Yres)
                        % Y_loc: weighting factor for points between gridpoints
                        Y_loc = (Y(i,j)-ymin-(q-1)*dy)/dy;
                    else
                        Y_loc = 1;
                    end
                        
                    v1 = X_loc*Y_loc;
                    v2 = X_loc*(1-Y_loc);
                    v3 = (1-X_loc)*Y_loc;
                    v4 = (1-X_loc)*(1-Y_loc);
                        
                    if (p<Xres) && (q<Yres)
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
                           
                    else if (p<Xres) && (q==Yres)
                            if (q<=0)
                                l=q+180;
                            else
                                l=q-180;
                            end
                            
                            f1x = U(q,p);
                            f2x = U(q,p);
                            f3x = U(q,p+1);
                            f4x = U(q,p+1);
                                
                            f1y = V(q,p);
                            f2y = V(q,p);
                            f3y = V(q,p+1);
                            f4y = V(q,p+1);
                                
                        else if (p==Xres) && (q<Yres)
                                f1x = U(q,p);
                                f2x = U(q+1,p);
                                f3x = U(q,p-(Xres-1));
                                f4x = U(q+1,p-(Xres-1));
                                    
                                f1y = V(q,p);
                                f2y = V(q+1,p);
                                f3y = V(q,p-(Xres-1));
                                f4y = V(q+1,p-(Xres-1));
                                    
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
                    
                    
                    X(i,j) = X(i,j) + dxdt*dt*computedirection_flag;
                    Y(i,j) = Y(i,j) + dydt*dt*computedirection_flag;
                    
                     % Particles wrap around in the x direction
                if CloseDomain_flag(1) == 1
                    if X(i,j)<xmin
                        X(i,j) = X(i,j)+dx*(Xres-1)*(1+abs(fix((X(i,j)-xmin)/((Xres-1)*dx)))); %X(i,j) could be much smaller than xmin
                    else if X(i,j)>xmax
                        X(i,j) = X(i,j)-dx*(Xres-1)*(1+abs(fix((X(i,j)-xmax)/((Xres-1)*dx)))); %X(i,j) could be much bigger than xmax
                    else
                        X(i,j) = X(i,j);
                        end
                    end
                 end
                % Particles wrap around in the y driection
                if CloseDomain_flag(2) == 1
                    if X(i,j) < 0
                        sign_id = 1;
                     else
                        sign_id = 0;
                     end
                    if Y(i,j)<ymin
                        X(i,j) = X(i,j)-(xmax-xmin)/2*(-1)^sign_id; %v2 updated here for the latitude boundary
                        Y(i,j) = 2*ymin-Y(i,j); % if Y(i,j) is much smallerthan ymin, need to fix!!!!!! but for this flow, need not right now
                    else if Y(i,j)>ymax 
                        X(i,j) = X(i,j)-(xmax-xmin)/2*(-1)^sign_id; %v2 updated here for the latitude boundary;
                        Y(i,j) = 2*ymax-Y(i,j);
                    else
                        X(i,j) = X(i,j);
                        Y(i,j) = Y(i,j);
                        end
                    end
                end
                
                    X_end(i,j) = X(i,j);
                    Y_end(i,j) = Y(i,j);
                
                end
            end
        end % for j = 1:Xres
    end % for i = 1:Yres
    
%-----------------------3D Computation
elseif dimension_flag == [1, 1, 1]
    % Preallocating for speed
    X_end = zeros(Xres, Yres, Zres);%Yres,Xres,Zres);
    Y_end = zeros(Xres, Yres, Zres);%Yres,Xres,Zres);
    Z_end = zeros(Xres, Yres, Zres);%Yres,Xres,Zres);

    % Should be able to vectorize the check for nans, commented out below,
    % if needed, here.
    % nanrows = find(isnan(X || Y || Z));
    % X_end(nanrows) = X(nanrows);
    % Y_end(nanrows) = Y(nanrows);
    % Z_end(nanrows) = Z(nanrows);

    % Next wrap any particles that are out of range in a closed domain.
    % Particles wrap around in the x direction
    if CloseDomain_flag(1) == 1
	x_outofrangerows = find(X < xmin);
	X(x_outofrangerows) = X(x_outofrangerows) + dx*(Xres-1) * (1 + abs(fix((X(x_outofrangerows) - xmin)/((Xres - 1) * dx))));
	clear x_outofrangerows
	x_outofrangerows = find(X > xmax);
        X(x_outofrangerows) = X(x_outofrangerows) - dx*(Xres-1) * (1 + abs(fix((X(x_outofrangerows) - xmax)/((Xres - 1) * dx)))); %X(i,j) could be much bigger than xmax
    end % if CloseDomain_flag(1) == 1

    % Particles wrap around in the y (latitude) direction
    if CloseDomain_flag(2) == 1
	    clear y_oorrows
	% Check if longitude is negative and select the power of -1 to be used.
	sign_id = (X < xmin);
	% if below -90 deg lat in the southern hemisphere.
        y_oorrows = find(Y < ymin);%if Y(i,j,k) < ymin
	if ~isempty(y_oorrows)
        % Bring longitude to the opposite hemisphere by shifting 180 deg.
	X(y_oorrows) = mod(X(y_oorrows), 180) + 180*(X(y_oorrows) == mod(X(y_oorrows), 180));
	Y(y_oorrows) = ymin - rem(Y(y_oorrows), ymin);
	end % if ~isempty(y_oorrows)
	
	clear y_oorrows
	y_oorrows = find(Y > ymax);
	if ~isempty(y_oorrows)
        % Bring longitude to the opposite hemisphere by shifting 180 deg.
	X(y_oorrows) = mod(X(y_oorrows), 180) + 180*(X(y_oorrows) == mod(X(y_oorrows), 180));
	Y(y_oorrows) = ymax - rem(Y(y_oorrows), ymax);
    	end % if ~isempty(y_oorrows)
    end % if CloseDomain_flag(2) == 1

    % Atmosphere is not closed so shouldn't ever execute this part.
    if CloseDomain_flag(3) == 1
	    keyboard
    end

    % Now that the positions have been wrapped, recompute indices.   
    % p: column index generated based on a given X coordinate
    % q: row index generated based on a given Y coordinate
    % l: third index generated based on a given coordinate
    % Example: if minimum x value is 5, then x=5 returns a p value of 1
    p = fix((X - xmin)/dx)+1;
    q = fix((Y - ymin)/dy)+1;
    for i = 1:numel(p)
    % Not enough to find nearest.  Need to find nearest lower, to recreate fix. SDB 12/12/25 I think "lower" means lower altitude, so take max index when two heights are equidistant.
    % 1/26/26 SDB Since I have to keep dimensions of X, Y, Z as [nlat, nlon, nht] to match U, V, W, the indexing of those arrays here must go as zarray(q, p, l).
    % 2/16/26 SDB Repeated the longitudes 0 and 360 in the input arrays.
    	l(i) = max(find(abs(Z(i) - zarray(q(i),p(i),:)) == min(abs(Z(i) - zarray(q(i),p(i),:)))));%53;%fix((z-zmin)/dz)+1;
    end % for i = 1:numel(p)
%    toc % Took about 12 s in interactive mode.
    l = reshape(l, size(p));
    
    % Checks if coordinate indices are outside of the domain otherwise it adds a velocity to the coordinate
%    ind_oorrows = find((p<=0)||(p>Xres)||(q<=0)||(q>Yres)||(l<=0)||(l>Zres));
%
%    X_end(ind_oorrows) = X(ind_oorrows);
%    Y_end(ind_oorrows) = Y(ind_oorrows);
%    Z_end(ind_oorrows) = Z(ind_oorrows);    
%[p, q, l]
			
    % Set the fractions for tri-linear interpolation.    
    ind_oorrows = find(p >= Xres);
    % X_loc: weighting factor for points between gridpoints
    X_loc(ind_oorrows) = 1;
    X_loc = (X - xmin - (p-1)*dx)/dx;
    clear ind_oorrows
    % Y_loc: weighting factor for points between gridpoints
    ind_oorrows = find(q >= Yres);
    
    Y_loc(ind_oorrows) = 1;
    Y_loc = (Y - ymin - (q-1)*dy)/dy;
                           
    clear ind_oorrows
    % Z_loc: weighting factor for points between gridpoints
    % For any points beyond the height grid, set the fraction to 1 so it will just stay at that edge. This is not the same as the finite domain fix, but can deal with it later. 1/26/26 SDB
    % indices of the out-of-range rows
    ind_oorrows = find(l >= Zres);
                            
    Z_loc(ind_oorrows) = 1;
    % Try this next line to see if we can have it parallel the other dims.
    % dz is not defined at this point so won't work.
    %    Z_loc = (Z - zmin - (l-1)*dz)/dz;
    % Loop through the ones that are in-range, not out of bounds.
    ind_irrows = find(l < Zres);
    % Modify for passing the whole 3d grid.
    for i = 1:numel(ind_irrows)
   
	% Z_loc represents the fraction of the distance between the two nearest heights, for linear interpolation. 
	Z_loc(ind_irrows(i)) = (Z(ind_irrows(i)) - zarray( q(ind_irrows(i)),p(ind_irrows(i)), l(ind_irrows(i))))/...
	(zarray(q(ind_irrows(i)), p(ind_irrows(i)), l(ind_irrows(i))+1) - zarray(q(ind_irrows(i)), p(ind_irrows(i)),l(ind_irrows(i))));
	%Z_loc(ind_irrows(i)) = (Z(ind_irrows(i)) - zarray(p(ind_irrows(i)),q(ind_irrows(i)),l(ind_irrows(i))))/(zarray(p(ind_irrows(i)),q(ind_irrows(i)), l(ind_irrows(i))+1) - zarray(p(ind_irrows(i), q(ind_irrows(i)), l(ind_irrows(i))));
    end
   % Or if we don't use the analogous line, reshape here.
   Z_loc = reshape(Z_loc, size(X_loc));
   
   % The fractional weights.   
   v1 = (1-X_loc) .* (1-Y_loc) .* (1-Z_loc);
   v2 = X_loc .* (1-Y_loc) .* (1-Z_loc);
   v3 = (1-X_loc) .* Y_loc .* (1-Z_loc);
   v4 = X_loc .* Y_loc .* (1-Z_loc);
   v5 = (1-X_loc) .* (1-Y_loc) .* Z_loc;
   v6 = X_loc .* (1-Y_loc) .* Z_loc;
   v7 = (1-X_loc) .* Y_loc .* Z_loc;
   v8 = X_loc .* Y_loc .* Z_loc;
%                            Z_loc = (Z(i,j,k)-zmin-(l-1)*dz)/dz;

% Initialize nearby grid speeds to have the right matrix dimensions.
A111_x = zeros(size(v1));
A121_x = zeros(size(v1));
A122_x = zeros(size(v1));
A112_x = A111_x;
A211_x = A111_x;
A212_x = A111_x;
A221_x = A111_x;
A222_x = A111_x;
A111_y = zeros(size(v1));
A121_y = zeros(size(v1));
A122_y = zeros(size(v1));
A112_y = A111_x;
A211_y = A111_x;
A212_y = A111_x;
A221_y = A111_x;
A222_y = A111_x;

A111_z = zeros(size(v1));
A121_z = zeros(size(v1));
A122_z = zeros(size(v1));
A112_z = A111_x;
A211_z = A111_x;
A212_z = A111_x;
A221_z = A111_x;
A222_z = A111_x;
   for i = 1:numel(p)                        
%if mod(i,1e3) == 0
%	i
%end
   %tic%	if p(i) < Xres && q(i) < Yres && l(i) < Zres      
	plower = p(i);
	qlower = q(i);
	llower = l(i);

	if p(i) == Xres
		pupper = p(i);
	else
		pupper = p(i)+1;
	end
	if q(i) == Yres
		qupper = q(i);
	else
		qupper = q(i)+1;
	end
	if l(i) == Zres
		lupper = l(i);
	else
		lupper = l(i)+1;
	end		
		
	%A111_x-A222_x/A111_7-A222_y/A111_z-A222_z: the chosen point and
	% neighboring points based on its position in the grid (for 3D these 
	% are the point itself and 7 other points which are all 
	% combinations of adding 1 to the x, y, and z position)
	% Example: If the point is (0,0,0), then the neighboring points are 
	% (0,0,1), (0,1,0), (0,1,1), (1,0,0), (1,0,1),(1,1,0), and (1,1,1).  
	% Exemptions are made the same way they are in the 2D case (if
	% neighboring points are off the grid, they are not used for the
	% interpolation).
	A111_x(i) = U(qlower, plower, llower);
	A112_x(i) = U(qlower, plower, lupper);
	A211_x(i) = U(qupper, plower, llower);
	A212_x(i) = U(qupper, plower, lupper);
	A121_x(i) = U(qlower, pupper, llower);
	A122_x(i) = U(qlower, pupper, lupper);
	A221_x(i) = U(qupper, pupper, llower);
	A222_x(i) = U(qupper, pupper, lupper);

	A111_y(i) = V(qlower, plower, llower);
	A112_y(i) = V(qlower, plower, lupper);
	A211_y(i) = V(qupper, plower, llower);
	A212_y(i) = V(qupper, plower, lupper);
	A121_y(i) = V(qlower, pupper, llower);
	A122_y(i) = V(qlower, pupper, lupper);
	A221_y(i) = V(qupper, pupper, llower);
	A222_y(i) = V(qupper, pupper, lupper);

	A111_z(i) = W(qlower, plower, llower);
	A112_z(i) = W(qlower, plower, lupper);
	A211_z(i) = W(qupper, plower, llower);
	A212_z(i) = W(qupper, plower, lupper);
	A121_z(i) = W(qlower, pupper, llower);
	A122_z(i) = W(qlower, pupper, lupper);
	A221_z(i) = W(qupper, pupper, llower);
	A222_z(i) = W(qupper, pupper,lupper);
    end 
	% Interpolated velocities.
	dxdt = A111_x.*v1 + A121_x.*v2 + A211_x.*v3 + A221_x.*v4 + A112_x.*v5 + A122_x.*v6 + A212_x.*v7 + A222_x.*v8;    
   
	dydt = A111_y.*v1 + A121_y.*v2 + A211_y.*v3 + A221_y.*v4 + A112_y.*v5 + A122_y.*v6 + A212_y.*v7 + A222_y.*v8;
   
	dzdt = A111_z.*v1 + A121_z.*v2 + A211_z.*v3 + A221_z.*v4 + A112_z.*v5 + A122_z.*v6 + A212_z.*v7 + A222_z.*v8;

	% Updated positions of all points.
	% 1/26/26 SDB Because U/V/W = [nlat, nlon, nht] reshaped, so X/Y/Z need to be [nlat, nlon, nht] reshaped.	
	X_end = X + dxdt*dt*computedirection_flag;
	Y_end = Y + dydt*dt*computedirection_flag;
	Z_end = Z + dzdt*dt*computedirection_flag;
end

% 2/10/26 SDB Putting in a fix in case the end points go beyond [0, 360]. May need to do same for latitude following procedure above [-90, 90].
X_end = X_end + 360 * (X_end < 0) - 360 * (X_end > 360);
%    for i=1:Xres%Yres
%        for j=1:Yres%Xres
%            for k=1:Zres
%                if isnan(X(i,j,k))||isnan(Y(i,j,k))||isnan(Z(i,j,k)) 
%                    % Checks if coordinates of particle velocity are NaN otherwise it adds a velocity to the coordinate
%                    X_end(i,j,k) = X(i,j,k);
%                    Y_end(i,j,k) = Y(i,j,k);
%                    Z_end(i,j,k) = Z(i,j,k);
%                else
                    
                % 3D
                % Particles wrap around in the x direction
%                 if CloseDomain_flag(1) == 1
%                     if X(i,j,k)<xmin
%                         X(i,j,k) = X(i,j,k)+dx*(Xres-1)*(1+abs(fix((X(i,j,k)-xmin)/((Xres-1)*dx)))); %X(i,j) could be much smaller than xmin
%                     elseif X(i,j,k)>xmax
%                         X(i,j,k) = X(i,j,k)-dx*(Xres-1)*(1+abs(fix((X(i,j,k)-xmax)/((Xres-1)*dx)))); %X(i,j) could be much bigger than xmax
%%                     else
%%                         X(i,j,k) = X(i,j,k);
%                     end
%%                     end
%                  end % if CloseDomain_flag(1) == 1

%                 % Particles wrap around in the y (latitude) direction
%                 if CloseDomain_flag(2) == 1
%			 % Check if longitude is negative and select the power of -1 to be used.
%			 sign_id = (X(i,j,k) < xmin);
%		% if below -90 deg lat in the southern hemisphere.
%                     if Y(i,j,k) < ymin
%			 % Bring to the opposite hemisphere by shifting 180 deg.
%			 X(i,j,k) = mod(X(i,j,k), 180) + 180*(X(i,j,k) == mod(X(i,j,k), 180));
%			 Y(i,j,k) = 2* ymin - Y(i,j,k);
%%                         Y(i,j,k) = 2*ymin-Y(i,j,k); % if Y(i,j) is much smallerthan ymin, need to fix!!!!!! but for this flow, need not right now
%                     elseif Y(i,j,k) > ymax
%			    X(i,j,k) = mod(X(i,j,k), 180) + 180*(X(i,j,k) == mod(X(i,j,k), 180));
%			   Y(i,j,k) = 2 * ymax - Y(i,j,k); 
%                     end
%		end % if CloseDomain_flag(2) == 1

% 3D Vertical, not closed in an atmosphere, so should never get executed.
%                if CloseDomain_flag(3) == 1
%                    if X(i,j,k) < xmin
%                        X(i,j,k) = X(i,j,k) + dx*(Xres-1)*(1+abs(fix((X(i,j,k)-xmin)/((Xres-1)*dx)))); %X(i,j) could be much smaller than xmin
%                    elseif X(i,j,k) > xmax
%                        X(i,j,k) = X(i,j,k) - dx*(Xres-1)*(1+abs(fix((X(i,j,k)-xmax)/((Xres-1)*dx)))); %X(i,j) could be much bigger than xmax
%                    %else
%                    %    X(i,j,k) = X(i,j,k);
%                    %    end
%                    end
%                 end % if CloseDomain_flag(3) == 1

%                % Particles wrap around in the y driection
%                if CloseDomain_flag(2) == 1
%                    %if X(i,j,k) < 0
%                        sign_id = (X(i,j,k) < xmin);%1;
%                    % else
%                    %    sign_id = 0;
%                    % end
%		    % If below -90 deg lat in the southern hemisphere
%		    if Y(i,j,k) < ymin
%			% Bring to the opposite hemisphere by shifting 180 deg
%			X(i,j,k) = mod(X(i,j,k), 180) + 180 * (X(i,j,k) == mod(X(i,j,k), 180));
%			Y(i,j,k) = 2 * ymin - Y(i,j,k);
%		elseif Y(i,j,k) > ymax
%			% Bring to the opposite hemisphere by shifting 180 deg
%			X(i,j,k) = mod(X(i,j,k), 180) + 180 * (X(i,j,k) == mod(X(i,j,k), 180));
%			Y(i,j,k) = 2 * ymin - Y(i,j,k);
%		end
			% 9/23/25 SDB thinks this logic is incorrect for true 3D calculations because 
% (x,z) will not be (lon, lat) but instead (lon, plev).
% 		     if Z(i,j,k)<zmin
%                        X(i,j,k) = X(i,j,k)-(xmax-xmin)/2*(-1)^sign_id; %v2 updated here for the latitude boundary
%                        Z(i,j,k) = 2*zmin-Z(i,j,k); % if Y(i,j) is much smallerthan ymin, need to fix!!!!!! but for this flow, need not right now
%                    else if Z(i,j,k)>zmax 
%                        X(i,j,k) = X(i,j,k)-(xmax-xmin)/2*(-1)^sign_id; %v2 updated here for the latitude boundary;
%                        Z(i,j,k) = 2*zmax-Z(i,j,k);
%                         else
%                        X(i,j,k) = X(i,j,k);
%                        Z(i,j,k) = Z(i,j,k);
%                        end
%                    end
%                end
%                 % Now that the positions have been wrapped, recompute indices.   
%                    % p: column index generated based on a given X coordinate
%                    % q: row index generated based on a given Y coordinate
%                    % l: third index generated based on a given coordinate
%                    % Example: if minimum x value is 5, then x=5 returns a p value of 1
%                    p = fix((X(i,j,k)-xmin)/dx)+1;
%                    q = fix((Y(i,j,k)-ymin)/dy)+1;
%% Not enough to find nearest.  Need to find nearest lower, to recreate fix.
%l = find(abs(Z(i,j,k) - zarray(p,q,:)) == min(abs(Z(i,j,k) - zarray(p,q,:))));%53;%fix((z-zmin)/dz)+1;
%[i, j, k, p, q, l]
%disp([q p l])
%l = max(find(zi - zarray(q,p,:) >= 0));%53;%fix((z-zmin)/dz)+1;
    % Checks if the point is outside of the domain otherwise it adds a velocity to the point
%                    l = fix((Z(i,j,k)-zmin)/dz)+1;
%if l > 126 || isempty(l)
%keyboard
%	X_end(i,j,k) = X(i,j,k);
%	Y_end(i,j,k) = Y(i,j,k);
%	Z_end(i,j,k) = Z(i,j,k);
%	return
%end
%                    if (p<=0)||(p>Xres)||(q<=0)||(q>Yres)||(l<=0)||(l>Zres)
%                        % Checks if coordinates are outside of the domain otherwise it adds a velocity to the coordinate
%                        X_end(i,j,k) = X(i,j,k);
%                        Y_end(i,j,k) = Y(i,j,k);
%                        Z_end(i,j,k) = Z(i,j,k);    
%                    else
%			% Set the fractions for tri-linear interpolation.
%                        if p>=Xres
%                            X_loc = 1;
%                        else
%                            % X_loc: weighting factor for points between gridpoints
%                            X_loc = (X(i,j,k)-xmin-(p-1)*dx)/dx;
%                        end
%                        if q>=Yres
%                            Y_loc = 1;
%                        else
%                            % Y_loc: weighting factor for points between gridpoints
%                            Y_loc = (Y(i,j,k)-ymin-(q-1)*dy)/dy;
%                        end
%                        if l>=Zres
%                            Z_loc = 1;
%                        else
%                            % Z_loc: weighting factor for points between gridpoints
%			    % Modify for passing the whole 3d grid.
%			    Z_loc = (Z(i,j,k) - zarray(p,q,l))/(zarray(p,q, l+1) - zarray(l));
%%                            Z_loc = (Z(i,j,k)-zmin-(l-1)*dz)/dz;
%                        end
                        
%                        v1 = (1-X_loc)*(1-Y_loc)*(1-Z_loc);
%                        v2 = X_loc*(1-Y_loc)*(1-Z_loc);
%                        v3 = (1-X_loc)*Y_loc*(1-Z_loc);
%                        v4 = X_loc*Y_loc*(1-Z_loc);
%                        v5 = (1-X_loc)*(1-Y_loc)*Z_loc;
%                        v6 = X_loc*(1-Y_loc)*Z_loc;
%                        v7 = (1-X_loc)*Y_loc*Z_loc;
%                        v8 = X_loc*Y_loc*Z_loc;
%                        if p(i) < Xres && q(i) < Yres && l(i) < Zres      
                            %A111_x-A222_x/A111_7-A222_y/A111_z-A222_z: the chosen point and
                            % neighboring points based on its position in the grid (for 3D these 
                            % are the point itself and 7 other points which are all 
                            % combinations of adding 1 to the x, y, and z position)
                            % Example: If the point is (0,0,0), then the neighboring points are 
                            % (0,0,1), (0,1,0), (0,1,1), (1,0,0), (1,0,1),(1,1,0), and (1,1,1).  
                            % Exemptions are made the same way they are in the 2D case (if
                            % neighboring points are off the grid, they are not used for the
                            % interpolation).
%                            A111_x = U(q,p,l);
%                            A112_x = U(q,p,l+1);
%                            A211_x = U(q+1,p,l);
%                            A212_x = U(q+1,p,l+1);
%                            A121_x = U(q,p+1,l);
%                            A122_x = U(q,p+1,l+1);
%                            A221_x = U(q+1,p+1,l);
%                            A222_x = U(q+1,p+1,l+1);
%    
%                            A111_y = V(q,p,l);
%                            A112_y = V(q,p,l+1);
%                            A211_y = V(q+1,p,l);
%                            A212_y = V(q+1,p,l+1);
%                            A121_y = V(q,p+1,l);
%                            A122_y = V(q,p+1,l+1);
%                            A221_y = V(q+1,p+1,l);
%                            A222_y = V(q+1,p+1,l+1);
%    
%                            A111_z = W(q,p,l);
%                            A112_z = W(q,p,l+1);
%                            A211_z = W(q+1,p,l);
%                            A212_z = W(q+1,p,l+1);
%                            A121_z = W(q,p+1,l);
%                            A122_z = W(q,p+1,l+1);
%                            A221_z = W(q+1,p+1,l);
%                            A222_z = W(q+1,p+1,l+1);
%                    
%                        elseif p<Xres && q<Yres && l==Zres
%                            A111_x = U(q,p,l);
%                            A112_x = U(q,p,l);
%                            A211_x = U(q+1,p,l);
%                            A212_x = U(q+1,p,l);
%                            A121_x = U(q,p+1,l);
%                            A122_x = U(q,p+1,l);
%                            A221_x = U(q+1,p+1,l);
%                            A222_x = U(q+1,p+1,l);
%    
%                            A111_y = V(q,p,l);
%                            A112_y = V(q,p,l);
%                            A211_y = V(q+1,p,l);
%                            A212_y = V(q+1,p,l);
%                            A121_y = V(q,p+1,l);
%                            A122_y = V(q,p+1,l);
%                            A221_y = V(q+1,p+1,l);
%                            A222_y = V(q+1,p+1,l);
%    
%                            A111_z = W(q,p,l);
%                            A112_z = W(q,p,l);
%                            A211_z = W(q+1,p,l);
%                            A212_z = W(q+1,p,l);
%                            A121_z = W(q,p+1,l);
%                            A122_z = W(q,p+1,l);
%                            A221_z = W(q+1,p+1,l);
%                            A222_z = W(q+1,p+1,l);
%    
%                        elseif p<Xres && q==Yres && l<Zres
%                            A111_x = U(q,p,l);
%                            A112_x = U(q,p,l+1);
%                            A211_x = U(q,p,l);
%                            A212_x = U(q,p,l+1);
%                            A121_x = U(q,p+1,l);
%                            A122_x = U(q,p+1,l+1);
%                            A221_x = U(q,p+1,l);
%                            A222_x = U(q,p+1,l+1);
%    
%                            A111_y = V(q,p,l);
%                            A112_y = V(q,p,l+1);
%                            A211_y = V(q,p,l);
%                            A212_y = V(q,p,l+1);
%                            A121_y = V(q,p+1,l);
%                            A122_y = V(q,p+1,l+1);
%                            A221_y = V(q,p+1,l);
%                            A222_y = V(q,p+1,l+1);
%    
%                            A111_z = W(q,p,l);
%                            A112_z = W(q,p,l+1);
%                            A211_z = W(q,p,l);
%                            A212_z = W(q,p,l+1);
%                            A121_z = W(q,p+1,l);
%                            A122_z = W(q,p+1,l+1);
%                            A221_z = W(q,p+1,l);
%                            A222_z = W(q,p+1,l+1);
%    
%                        elseif p==Xres && q<Yres && l<Zres
%                            A111_x = U(q,p,l);
%                            A112_x = U(q,p,l+1);
%                            A211_x = U(q+1,p,l);
%                            A212_x = U(q+1,p,l+1);
%                            A121_x = U(q,p,l);
%                            A122_x = U(q,p,l+1);
%                            A221_x = U(q+1,p,l);
%                            A222_x = U(q+1,p,l+1);
%    
%                            A111_y = V(q,p,l);
%                            A112_y = V(q,p,l+1);
%                            A211_y = V(q+1,p,l);
%                            A212_y = V(q+1,p,l+1);
%                            A121_y = V(q,p,l);
%                            A122_y = V(q,p,l+1);
%                            A221_y = V(q+1,p,l);
%                            A222_y = V(q+1,p,l+1);
%    
%                            A111_z = W(q,p,l);
%                            A112_z = W(q,p,l+1);
%                            A211_z = W(q+1,p,l);
%                            A212_z = W(q+1,p,l+1);
%                            A121_z = W(q,p,l);
%                            A122_z = W(q,p,l+1);
%                            A221_z = W(q+1,p,l);
%                            A222_z = W(q+1,p,l+1);
%    
%                        elseif p<Xres && q==Yres && l==Zres
%                            A111_x = U(q,p,l);
%                            A112_x = U(q,p,l);
%                            A211_x = U(q,p,l);
%                            A212_x = U(q,p,l);
%                            A121_x = U(q,p+1,l);
%                            A122_x = U(q,p+1,l);
%                            A221_x = U(q,p+1,l);
%                            A222_x = U(q,p+1,l);
%    
%                            A111_y = V(q,p,l);
%                            A112_y = V(q,p,l);
%                            A211_y = V(q,p,l);
%                            A212_y = V(q,p,l);
%                            A121_y = V(q,p+1,l);
%                            A122_y = V(q,p+1,l);
%                            A221_y = V(q,p+1,l);
%                            A222_y = V(q,p+1,l);
%    
%                            A111_z = W(q,p,l);
%                            A112_z = W(q,p,l);
%                            A211_z = W(q,p,l);
%                            A212_z = W(q,p,l);
%                            A121_z = W(q,p+1,l);
%                            A122_z = W(q,p+1,l);
%                            A221_z = W(q,p+1,l);
%                            A222_z = W(q,p+1,l);
%    
%                        elseif p==Xres && q<Yres && l==Zres
%                            A111_x = U(q,p,l);
%                            A112_x = U(q,p,l);
%                            A211_x = U(q+1,p,l);
%                            A212_x = U(q+1,p,l);
%                            A121_x = U(q,p,l);
%                            A122_x = U(q,p,l);
%                            A221_x = U(q+1,p,l);
%                            A222_x = U(q+1,p,l);
%    
%                            A111_y = V(q,p,l);
%                            A112_y = V(q,p,l);
%                            A211_y = V(q+1,p,l);
%                            A212_y = V(q+1,p,l);
%                            A121_y = V(q,p,l);
%                            A122_y = V(q,p,l);
%                            A221_y = V(q+1,p,l);
%                            A222_y = V(q+1,p,l);
%    
%                            A111_z = W(q,p,l);
%                            A112_z = W(q,p,l);
%                            A211_z = W(q+1,p,l);
%                            A212_z = W(q+1,p,l);
%                            A121_z = W(q,p,l);
%                            A122_z = W(q,p,l);
%                            A221_z = W(q+1,p,l);
%                            A222_z = W(q+1,p,l);
%     
%                        elseif p==Xres && q==Yres && l<Zres
%                            A111_x = U(q,p,l);
%                            A112_x = U(q,p,l+1);
%                            A211_x = U(q,p,l);
%                            A212_x = U(q,p,l+1);
%                            A121_x = U(q,p,l);
%                            A122_x = U(q,p,l+1);
%                            A221_x = U(q,p,l);
%                            A222_x = U(q,p,l+1);
%    
%                            A111_y = V(q,p,l);
%                            A112_y = V(q,p,l+1);
%                            A211_y = V(q,p,l);
%                            A212_y = V(q,p,l+1);
%                            A121_y = V(q,p,l);
%                            A122_y = V(q,p,l+1);
%                            A221_y = V(q,p,l);
%                            A222_y = V(q,p,l+1);
%    
%                            A111_z = W(q,p,l);
%                            A112_z = W(q,p,l+1);
%                            A211_z = W(q,p,l);
%                            A212_z = W(q,p,l+1);
%                            A121_z = W(q,p,l);
%                            A122_z = W(q,p,l+1);
%                            A221_z = W(q,p,l);
%                            A222_z = W(q,p,l+1);
%    
%                        else
%                            A111_x = U(q,p,l);
%                            A112_x = U(q,p,l);
%                            A211_x = U(q,p,l);
%                            A212_x = U(q,p,l);
%                            A121_x = U(q,p,l);
%                            A122_x = U(q,p,l);
%                            A221_x = U(q,p,l);
%                            A222_x = U(q,p,l);
%    
%                            A111_y = V(q,p,l);
%                            A112_y = V(q,p,l);
%                            A211_y = V(q,p,l);
%                            A212_y = V(q,p,l);
%                            A121_y = V(q,p,l);
%                            A122_y = V(q,p,l);
%                            A221_y = V(q,p,l);
%                            A222_y = V(q,p,l);
%    
%                            A111_z = W(q,p,l);
%                            A112_z = W(q,p,l);
%                            A211_z = W(q,p,l);
%                            A212_z = W(q,p,l);
%                            A121_z = W(q,p,l);
%                            A122_z = W(q,p,l);
%                            A221_z = W(q,p,l);
%                            A222_z = W(q,p,l);
%                        end
%                        
%                        dxdt = A111_x*v1 + A121_x*v2 + A211_x*v3 + A221_x*v4 + A112_x*v5 + A122_x*v6 + A212_x*v7 + A222_x*v8;
%                        dydt = A111_y*v1 + A121_y*v2 + A211_y*v3 + A221_y*v4 + A112_y*v5 + A122_y*v6 + A212_y*v7 + A222_y*v8;
%                        dzdt = A111_z*v1 + A121_z*v2 + A211_z*v3 + A221_z*v4 + A112_z*v5 + A122_z*v6 + A212_z*v7 + A222_z*v8;

%                         dxdt = A111_x*v1 + A211_x*v2 + A121_x*v3 + A221_x*v4 + A112_x*v5 * A212_x*v6 + A122_x*v7 + A222_x*v8;
%                         dydt = A111_y*v1 + A211_y*v2 + A121_y*v3 + A221_y*v4 + A112_y*v5 * A212_y*v6 + A122_y*v7 + A222_y*v8;
%                         dzdt = A111_z*v1 + A211_z*v2 + A121_z*v3 + A221_z*v4 + A112_z*v5 * A212_z*v6 + A122_z*v7 + A222_z*v8;
                        
%                        X_end(i,j,k) = X(i,j,k) + dxdt*dt*computedirection_flag;
%                        Y_end(i,j,k) = Y(i,j,k) + dydt*dt*computedirection_flag;
%                        Z_end(i,j,k) = Z(i,j,k) + dzdt*dt*computedirection_flag;
%                    end
%                end
%            end
    
%	    end
%    end
%end
%end
