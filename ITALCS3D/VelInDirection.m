% Extracts U,V,W velocities at a given time(U1,V1,W1) and at one timestep 
% forward or backwards (U2,V2,W2) depending on computation direction. 

function [U1,V1,W1,U2,V2,W2] = VelInDirection(flow,Xres,Yres,Zres,t,TimeStamp,computedirection_flag,dimension_flag)

% Inputs:
%       flow: velocity field matrix.
%       Xres/Yres/Zres: number of gridpoints in x,y,z-direction.
%       t: time within the integration slot(not an array its a specific time)
%       TimeStamp: increment of time in within an integration loop.
%       computedirection_flag: forward or backward computation (1,-1).
%       dimension_flag: 2D or 3D (0,1).
% Outputs
%       U1/U2/U3: x,y,z-direction velocities for one timestamp.
%       U2/V2/W2: x,y,z-direction velocities at the next time stamp.

% ts: index generated based on the time input
ts = fix(t/TimeStamp) + 1;

% flow_ts: velocities at the inputted time
flow_ts = flow(ts,:);
u = flow_ts(:,2:3:end-2); 
v = flow_ts(:,3:3:end-1);
w = flow_ts(:,4:3:end);

if dimension_flag == 0
    U1 = reshape(u,Xres,Yres)'; 
    V1 = reshape(v,Xres,Yres)';
    W1 = reshape(w,Xres,Yres)';
else
    
% 3D
%     U1 = reshape(u,Xres,Yres,Zres);
%     U1 = permute(U1,[2,1,3]);
%     V1 = reshape(v,Xres,Yres,Zres);
%     V1 = permute(V1,[2,1,3]);
%     W1 = reshape(w,Xres,Yres,Zres);
%     W1 = permute(W1,[2,1,3]);
    
% 3D Vert
%    % This first line extracts them to be nlon x nlev x nlat.
%    U1 = reshape(u,Xres,Zres,Yres);
%    % Next they seem to be made nlev x nlon x nlat.
%    U1 = permute(U1,[2,1,3]);
%    V1 = reshape(v,Xres,Zres,Yres);
%    V1 = permute(V1,[2,1,3]);
%    W1 = reshape(w,Xres,Zres,Yres);
%    W1 = permute(W1,[2,1,3]);
%   
%    % Then changed again to be nlat x nlon x nlev.
%    % These dimensions match the X, Y, Z spatial grid meshgrid dimensions.
%    % But I'm not sure they've been extracted in the correct order.
%    % SDB 5/23/25. 
%    U1 = permute(U1,[3,2,1]);
%    V1 = permute(V1,[3,2,1]);
%    W1 = permute(W1,[3,2,1]);

    % This first line extracts them to be nlon x nlat x nlev,
    % which is the order in which they were placed in the flow matrix.
    % (nlon changes most rapidly, nlat is the middle loop, and nlev is the outer loop).
    U1 = reshape(u,Xres,Yres,Zres);
    % Next they need to be made so that rows correspond to latitudes and
    % columns are longitudes, in order to match the meshgrid arrangement of
    % the X, Y, and Z matrices.
    U1 = permute(U1,[2,1,3]);
    V1 = reshape(v,Xres,Yres,Zres);
    V1 = permute(V1,[2,1,3]);
    W1 = reshape(w,Xres,Yres,Zres);
    W1 = permute(W1,[2,1,3]);
   
end

% Checks if we are at the end of a forward computation
if computedirection_flag == 1 && ts == size(flow,1)
    U2 = U1;
    V2 = V1;
    W2 = W1;
% Checks if we are at the end of a backwards computation
else if computedirection_flag == -1 && ts == 1
        U2 = U1;
        V2 = V1;
        W2 = W1;
    else
        % Gets the velocities at each gridpoint forward or backwards one time step
        flow_ts = flow(ts+1*computedirection_flag,:); 
        u = flow_ts(:,2:3:end-2); 
        v = flow_ts(:,3:3:end-1); 
        w = flow_ts(:,4:3:end); 
        if dimension_flag == 0
            U2 = reshape(u,Xres,Yres)'; 
            V2 = reshape(v,Xres,Yres)';
            W2 = reshape(w,Xres,Yres)';
        else
% 3D
             U2 = reshape(u,Xres,Yres,Zres);
             U2 = permute(U2,[2,1,3]);
             V2 = reshape(v,Xres,Yres,Zres);
             V2 = permute(V2,[2,1,3]);
             W2 = reshape(w,Xres,Yres,Zres);
             W2 = permute(W2,[2,1,3]);
            
% 3D Vertical
%            U2 = reshape(u,Xres,Zres,Yres);
%            U2 = permute(U2,[2,1,3]);
%            V2 = reshape(v,Xres,Zres,Yres);
%            V2 = permute(V2,[2,1,3]);
%            W2 = reshape(w,Xres,Zres,Yres);
%            W2 = permute(W2,[2,1,3]);
%            
%            U2 = permute(U2,[3,2,1]);
%            V2 = permute(V2,[3,2,1]);
%            W2 = permute(W2,[3,2,1]);
        end
    end
end
