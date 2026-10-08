% Function to create the soomth boundary layer of the local domain flow
% fields
function [ U_updated,V_updated ] = vel_modified(X,Y,Delta_xyz,fdmin_xyz,fdmax_xyz,U,V,FDFTLE_index,FTLE_flag)
%Note on function inputs:
%       X: the horizontal coordinate of the point on the grid
%       Y: the vertical coordinate of the point on the grid
%       Delta_xyz: distance from each boundary that finite domain will
%       smooth velocity(array input: [x,y,z])
%       fdmin_xyz: minimum index value of the x,y,z coordinates
%       fdmax_xyz: maximum index value of the x,y,z coordinates
%       U: the vertical velocity at the point (X,Y) on the grid
%       V: the horizontal velocity at the point (X,Y) on the grid
%       FDFTLE_index: finite domain edge smoothing on or off (1,0) for x,y,z 
%       boundaries (array input: [x,y,z])
%       FTLE_flag: No FTLE, regular FTLE, HWM14 mean velocity FDFTLE, or 
%       homoclinic linearized velocity FDFTLE (0,1,2,3).

% Max and Min values extracted from fdmax_xyz and fdmin_xyz
xmax=fdmax_xyz(1);
ymax=fdmax_xyz(2);
zmax=fdmax_xyz(3);
xmin=fdmin_xyz(1);
ymin=fdmin_xyz(2);
zmin=fdmin_xyz(3);

% Delta Values extracted from Delta_xyz
Deltax=Delta_xyz(1);
Deltay=Delta_xyz(2);
Deltaz=Delta_xyz(3);

% Generating the smoothing functions for all four boundaries
[v1,v2] = size(X);
g = linspace(-1,1,101);
f = g.^3 - 3*g + 2;

% Bottom Boundary
gp = linspace(ymin,ymin+Deltay,101);
polymin = polyfit(gp,flip(f/4),3);

% Top Boundary
gp = linspace(ymax-Deltay,ymax,101);
polymax = polyfit(gp,(f/4),3);

% Left Boundary
gp = linspace(xmin,xmin+Deltax,101);
Polymin = polyfit(gp,flip(f/4),3);

% Right Boundary
gp = linspace(xmax-Deltax,xmax,101);
Polymax = polyfit(gp,(f/4),3);

% Preallocation for speed
u_l = zeros(v1,v2);
v_l = zeros(v1,v2);

% Mean velocity FDFTLE or linearized velocity FDFTLE
if FTLE_flag == 2
    % For zonal/meridianal wind mean(U) = meridianal, mean (U') = zonal
%     U_l = zeros(v1,v2);
%     V_l = zeros(v1,v2);
    U_updated=U;
    V_updated=V;
%     G = mean(U);% Meridianal U mean
%     H = mean(V');% Zonal V mean
%     K = mean(U');% Zonal U mean
    [ U_l, V_l ] = vel_linear(X, Y, U, V);
else
    [ U_l, V_l ] = vel_linear(X, Y, U, V);% Gets linearized U and V velocities
end

for i = 1:v1
    for j = 1:v2
        x= X(i,j);
        y = Y(i,j);
        u = U(i,j);
        v = V(i,j);
        if FTLE_flag == 2
%             u_l(i,j) = G(j);% For zonal/meridianal wind mean(U) = meridianal, mean (U') = zonal
%             v_l(i,j) = H(i);% For zonal/meridianal wind mean(V) = meridianal, mean (V') = zonal
%             U_l = K(j);
        else
            u_l = U_l(i,j);
            v_l = V_l(i,j);
        end
        
        if y > (ymin+Deltay) && y < (ymax-Deltay)
                delta_y = Deltay^3;
        elseif y >= ymax || y <= ymin
                delta_y = 0;
        elseif y <= (ymin+Deltay)
                delta_y = (Deltay^3)*(polymin(1)*(y.^3) + polymin(2)*(y.^2) + polymin(3)*(y) + polymin(4));     
        elseif y >= (ymax-Deltay)
                delta_y = (Deltay^3)*(polymax(1)*(y.^3) + polymax(2)*(y.^2) + polymax(3)*(y) + polymax(4));  
        end
     
        if x > (xmin+Deltax) && x < (xmax-Deltax)
                delta_x = Deltax^3;
        elseif x >= xmax || x <= xmin
                delta_x = 0;
        elseif x <= (xmin+Deltax)
                delta_x = (Deltax^3)*(Polymin(1)*(x.^3) + Polymin(2)*(x.^2) + Polymin(3)*(x) + Polymin(4));     
        elseif x >= (xmax-Deltax)
                delta_x = (Deltax^3)*(Polymax(1)*(x.^3) + Polymax(2)*(x.^2) + Polymax(3)*(x) + Polymax(4));  
        end
        
        if FTLE_flag == 2
            if FDFTLE_index(1) == 1
                % Updates the x boundary velocities
                U_updated(i,j) = U_updated(i,j) +(U_l(i,j)-u)*(1-delta_x/(Deltax^3));
                V_updated(i,j) = V_updated(i,j) +(V_l(i,j)-v)*(1-delta_x/(Deltax^3));
            end
        
            if FDFTLE_index(2) == 1
                % Updates the y boundary velocities
                U_updated(i,j) = U_updated(i,j) +(U_l(i,j)-u)*(1-delta_y/(Deltay^3));
                V_updated(i,j) = V_updated(i,j) +(V_l(i,j)-v)*(1-delta_y/(Deltay^3));
            end
        else
            % Updates all 4 boundaries at the same time(specific to linearized case)
            U_updated(i,j) = u_l +(u-u_l)*delta_x*delta_y/(Deltax^3 * Deltay^3);
            V_updated(i,j) = v_l +(v-v_l)*delta_x*delta_y/(Deltax^3 * Deltay^3);
        end
    end
    end
    
end
