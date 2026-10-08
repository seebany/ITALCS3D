%function [Jmatrix] = JacobianGradient3datmo(X0,Y0,Z0,X,Y,Z,Xres,Yres,Zres,dimension_flag,Euclidean_flag)
%
% This function computes the Jacobian matrix (deformation matrix) in a 2D 
% or 3D domain in either euclidean or non-euclidean space.
%
% Inputs:
%       XO/Y0/Z0: initial x,y,z positions for every gridpoint and timestamp.
%       X/Y/Z: final x,y,z positions for every gridpoint and timestamp.
%       dimension_flag: 2D or 3D computation (0,1).
%       Euclidean_flag: Euclidean or Non-Euclidean computation (1,0)
% Outputs:
%       Jmatrix: Jacobian matrix (3x3xjxixk). Each particle has a unique position in the matrix, where an individual 3x3 matrix is saved.

function [Jmatrix] = JacobianGradient3datmo(X0,Y0,Z0,X,Y,Z,Xres,Yres,Zres,dimension_flag,Euclidean_flag)

if Euclidean_flag == 1
    if dimension_flag==0
        [dx_before,fy]=gradient(X0);
        [fx,dy_before]=gradient(Y0);
        [dx_after11,dx_after12]=gradient(X);
        [dy_after21,dy_after22]=gradient(Y);
    
        Jmatrix=[];
        Jmatrix(1,1,:,:)=dx_after11./dx_before;
        Jmatrix(1,2,:,:)=dx_after12./dy_before;
        Jmatrix(2,1,:,:)=dy_after21./dx_before;
        Jmatrix(2,2,:,:)=dy_after22./dy_before;
    else
        [dx_before,~,~] = gradient(X0);
        [~,dy_before,~] = gradient(Y0);
        [~,~,dz_before] = gradient(Z0);
        [dx_after11,dx_after12,dx_after13] = gradient(X);
        [dy_after21,dy_after22,dy_after23] = gradient(Y);
        [dz_after31,dz_after32,dz_after33] = gradient(Z);
    
        Jmatrix = [];
        Jmatrix(1,1,:,:,:) = dx_after11./dx_before;
        Jmatrix(1,2,:,:,:) = dx_after12./dy_before;
        Jmatrix(1,3,:,:,:) = dx_after13./dz_before;
        Jmatrix(2,1,:,:,:) = dy_after21./dx_before;
        Jmatrix(2,2,:,:,:) = dy_after22./dy_before;
        Jmatrix(2,3,:,:,:) = dy_after23./dz_before;
        Jmatrix(3,1,:,:,:) = dz_after31./dx_before;
        Jmatrix(3,2,:,:,:) = dz_after32./dy_before;
        Jmatrix(3,3,:,:,:) = dz_after33./dz_before;
    end
else % then Euclidean_flag ~= 1
	if dimension_flag==0
    [dx_before,~]=gradient(X0);
    [~,dy_before]=gradient(Y0);
    
%     for i=1:Yres
%         for j=1:Xres-1
%             if abs(X(i,j)-X(i,j+1)) > 300
%                 X(i,j+1)=X(i,j+1)+360;
%             end
%         end
%     end
    
    %  augment matrices so you don't have to use first/last row/column gradients
    % for x direction, the latitude boundary
    X = [X(1,round((Xres+1)/2):end),X(1,1:round((Xres-1)/2));X;X(end,round((Xres+1)/2):end),X(end,1:round((Xres-1)/2))];
    % for x direction, the longitude boundary
    X = [X(:,end),X,X(:,1)];
    % for y direction, the latitude boundary
    Y = [Y(1,round((Xres+1)/2):end),Y(1,1:round((Xres-1)/2));Y;Y(end,round((Xres+1)/2):end),Y(end,1:round((Xres-1)/2))];
    % for y direction, the longitude boundary
    Y = [Y(:,end),Y,Y(:,1)];
    [dx_after11,dx_after12]=gradient(X);
    [dy_after21,dy_after22]=gradient(Y);
    
    % pick the inner part of the gradient matrix
    dx_after11 = dx_after11(2:end-1,2:end-1);
    dx_after12 = dx_after12(2:end-1,2:end-1);
    dy_after21 = dy_after21(2:end-1,2:end-1);
    dy_after22 = dy_after22(2:end-1,2:end-1);

%     dx_after11 = dx_after11(:,2:end-1);
%     dx_after12 = dx_after12(:,2:end-1);
%     dy_after21 = dy_after21(:,2:end-1);
%     dy_after22 = dy_after22(:,2:end-1);
    %Defines the deformation gradient tensor (or the D matrix)
    
        %restrict the maximum distance
    for i = 1:Yres,
        for j = 1:Xres,
            if (abs(dx_after11(i,j)) > 90)
                dx_after11(i,j)=(sign(dx_after11(i,j)))*(179.75-(abs(dx_after11(i,j))));
            else
                dx_after11(i,j)=dx_after11(i,j);
            end
        end
    end
 
    
    for i = 1:Yres,
        for j = 1:Xres,
            if (abs(dx_after12(i,j)) > 90)
                dx_after12(i,j)=(sign(dx_after12(i,j)))*(179.75-(abs(dx_after12(i,j))));
            else
                dx_after12(i,j)=dx_after12(i,j);
            end
        end
    end
                
    
%     for i = 1:Yres,
%         for j = 1:Xres,
%             if (abs(dy_after21(i,j)) > 90)
%                 dy_after21(i,j)=(sign(dy_after21(i,j)))*(90-(abs(dy_after21(i,j))));
%             else
%                 dy_after21(i,j)=dy_after21(i,j);
%             end
%         end
%     end
% 
%     
%     for i = 1:Yres,
%         for j = 1:Xres,
%             if (abs(dy_after22(i,j)) > 90)
%                 dy_after22(i,j)=(sign(dy_after22(i,j)))*(180-(abs(dy_after22(i,j))));
%             else
%                 dy_after22(i,j)=dy_after22(i,j);
%             end
%         end
%     end 
    
    Jmatrix=[];
    %Calculations for the Deformation Gradient Tensor
    %In a sense, this is now 2x2xjxi matrix, same position of particle on the
    %grid, except now we tacked a matrix at each point
    %so when we want to describe D matrix of a particle at x=0 and y=0, we
    %would say (:,:,1,1) or at x=1 and y=1 would say (:,:,J+1,I+1)
    Jmatrix(1,1,:,:)=dx_after11./dx_before;
    Jmatrix(1,2,:,:)=dx_after12./dy_before;
    Jmatrix(2,1,:,:)=dy_after21./dx_before;
    Jmatrix(2,2,:,:)=dy_after22./dy_before;

    % 3D non-Euclidean (atmospheric calculation)
    elseif dimension_flag == [1 1 1]
	% SDB 12/1/25: In the 3d atmospheric version, 
	% I had made X0 = [nlon x nlat x nht]
	% so [144 x 96 x 126]. The gradient function permutes
	% the differencing of the first two dimensions.
	%[~,dx_before,~] = gradient(X0);
        %[dy_before,~,~] = gradient(Y0);
        %[~,~,dz_before] = gradient(Z0);
        %keyboard 
	% SDB 1/26/26: I had to revert X, Y, Z as
	% [nlat x nlon x nht] to match U, V, W. So changing this here.
	[dx_before,~,~] = gradient(X0);
        [~,dy_before,~] = gradient(Y0);
        [~,~,dz_before] = gradient(Z0);

	% 3D Vertical
	% SDB 12/1/25 Not sure what this was for, so comment out.
	% SDB 2/10/26 I think this is to augment the matrices so you don't
	% have to use the edge row/col gradients in the gradient calculation.
%	for i = 1:Yres
%	    %  update final matrix
%	    % for x direction, the latitude boundary
%	    P = [X(i,round((Xres+1)/2):end,1),X(i,1:round((Xres-1)/2),1);permute(X(i,:,:),[3,2,1]);X(i,round((Xres+1)/2):end,end),X(i,1:round((Xres-1)/2),end)];
%	    % for x direction, the longitude boundary
%	    XP(i,:,:) = [P(:,end),P,P(:,1)];
%	    % for y direction, the latitude boundary
%	    Q = [Z(i,round((Xres+1)/2):end,1),Z(i,1:round((Xres-1)/2),1);permute(Z(i,:,:),[3,2,1]);Z(i,round((Xres+1)/2):end,end),Z(i,1:round((Xres-1)/2),end)];
%	    % for y direction, the longitude boundary
%	    YQ(i,:,:) = [Q(:,end),Q,Q(:,1)];
%	end
%    
%	X=permute(XP,[1,3,2]);
%	Y=permute(YQ,[1,3,2]);
	% SDB 2/11/26 Attempting to close the FTLE gap between 357-360 deg.
	% 2/12/26 It doesn't appear to work, so comment out and repeat
	% velocities in the flow fields instead.
	% Repeat the matrix elements of the 0 longitude on the 360 side.
	%Xperm = permute(X, [2, 3, 1]);
	%Xaug = [Xperm; Xperm(1,:,:)];
	%X = permute(Xaug, [3, 1, 2]);
	%Yperm = permute(Y, [2, 3, 1]);
	%Yaug = [Yperm; Yperm(1,:,:)];
	%Y = permute(Yaug, [3, 1, 2]);
	%Zperm = permute(Z, [2, 3, 1]);
	%Zaug = [Zperm; Zperm(1,:,:)];
	%Z = permute(Zaug, [3, 1, 2]);
%----
    % These are [nlat x nlon x nht]
    [dx_after11,dx_after12,dx_after13]=gradient(X);
    [dy_after21,dy_after22,dy_after23]=gradient(Y);
    [dz_after31,dz_after32,dz_after33]=gradient(Z);
    
    % pick the inner part of the gradient matrix
    % SDB 12/1/25 Not sure why we need to do that. Maybe an edge effect?
    % Commenting out.
    % SDB 10/3/26 I now think this was to be used if the arrays were augmented so that edge gradients
    % could be computed the same as interior gradients with the gradient.m 
    % function, following the augmentation above that I commented on 2/10/26.
    %dx_after11 = dx_after11(:,2:end-1,2:end-1);
    %dx_after12 = dx_after12(:,2:end-1,2:end-1);
    %dx_after13 = dx_after13(:,2:end-1,2:end-1);
    %
    %dy_after21 = dy_after21(:,2:end-1,2:end-1);
    %dy_after22 = dy_after22(:,2:end-1,2:end-1);
    %dy_after23 = dy_after23(:,2:end-1,2:end-1);
    % SDB 1/28/26 Uncommenting, as it may help with the pole problem, which gives nans and infs. You only uncomment this after augmenting the spherical grid. That way the dimensions stay [nlat x nlon x nht]
    % SDB 2/12/26 This did not appear to close the gap in the plot
    %dx_after11 = dx_after11(:,1:end-1,:);
    %dx_after12 = dx_after12(:,1:end-1,:);
    %dx_after13 = dx_after13(:,1:end-1,:);
    %
    %dy_after21 = dy_after21(:,1:end-1,:);
    %dy_after22 = dy_after22(:,1:end-1,:);
    %dy_after23 = dy_after23(:,1:end-1,:);%
    
    %dz_after31 = dz_after31(:,1:end-1,:);
    %dz_after32 = dz_after32(:,1:end-1,:);
    %dz_after33 = dz_after33(:,1:end-1,:);

    %Defines the deformation gradient tensor (or the D matrix)
    
        %restrict the maximum distance
    % SDB 12/1/25 Not sure why we need to do that. Maybe an edge effect?
    % Commenting out.
    % SDB 1/26/26 I see now why it is needed: Ningchao and I discussed that going around a sphere, we would need to choose the shorter latitudinal and longitudinal difference.
% SDB 10/3/26 I'm trying to figure out why the test is for > 90 deg instead of 180
% SDB 10/3/26 Original implementation is > 90 deg.
% First try vectorizing to remove loops.
% Using bad_idx = find(abs(dx_after11) > 90); gives identical results to original implementation.
%bad_idx = find(abs(dx_after11) > 90);
%dx_after11(bad_idx) = sign(dx_after11(bad_idx)).*(180- abs(dx_after11(bad_idx)));
%clear bad_idx
%bad_idx = find(abs(dx_after12) > 90);
%dx_after12(bad_idx) = sign(dx_after12(bad_idx)).*(180- abs(dx_after12(bad_idx)));
%clear bad_idx
%bad_idx = find(abs(dx_after13) > 90);
%dx_after13(bad_idx) = sign(dx_after13(bad_idx)).*(180- abs(dx_after13(bad_idx)));
%clear bad_idx
% Next try if abs(dx) >180 deg, compute 360-abs(dx)
% This gives spurious yellow ridges. See in file 261003_WACCMX+DART_090125_0000_24hr_90km_ftle2d-vs-3d.  I still can't tell why.
bad_idx = find(abs(dx_after11) > 180);
dx_after11(bad_idx) = sign(dx_after11(bad_idx)).*(360- abs(dx_after11(bad_idx)));
clear bad_idx
bad_idx = find(abs(dx_after12) > 180);
dx_after12(bad_idx) = sign(dx_after12(bad_idx)).*(360- abs(dx_after12(bad_idx)));
clear bad_idx
bad_idx = find(abs(dx_after13) > 180);
dx_after13(bad_idx) = sign(dx_after13(bad_idx)).*(360- abs(dx_after13(bad_idx)));
clear bad_idx

%for k = 1:Zres,
%    for i = 1:Yres,
%        for j = 1:Xres,
%            if (abs(dx_after11(i,j,k)) > 90)
%%		    [i,j,k]
%%		    dx_after11(i,j,k)
%%            if (abs(dx_after12(i,j,k)) > 180)
%                dx_after11(i,j,k)=(sign(dx_after11(i,j,k)))*(180-(abs(dx_after11(i,j,k))));
%%            end
%            end
%            if (abs(dx_after12(i,j,k)) > 90)
%%		    [i,j,k]
%%		    dx_after12(i,j,k)
%%            if (abs(dx_after12(i,j,k)) > 180)
%                dx_after12(i,j,k)=(sign(dx_after12(i,j,k)))*(180-(abs(dx_after12(i,j,k))));
%%            end
%            end
%            if (abs(dx_after13(i,j,k)) > 90)
%%		    [i,j,k]
%%		    dx_after13(i,j,k)
%%            if (abs(dx_after13(i,j,k)) > 180)
%                dx_after13(i,j,k)=(sign(dx_after13(i,j,k)))*(180-(abs(dx_after13(i,j,k))));
%%            end
%            end
%        end
%    end
%end                

% 3D Vertical
%for i = 1:Yres,
%        for j = 1:Xres,
%            if (abs(dz_after31(i,j,k)) > 90)
%                dz_after31(i,j,k)=(sign(dz_after31(i,j,k)))*(90-(abs(dz_after31(i,j,k))));
%            else
%                dz_after31(i,j,k)=dz_after31(i,j,k);
%            end
%        end
%    end
%
%    
%    for i = 1:Yres,
%        for j = 1:Xres,
%            if (abs(dz_after32(i,j,k)) > 90)
%                dz_after32(i,j,k)=(sign(dz_after32(i,j,k)))*(90-(abs(dz_after32(i,j,k))));
%            else
%                dz_after32(i,j,k)=dz_after32(i,j,k);
%            end
%        end
%    end
%end
    
    Jmatrix=[];
    %Calculations for the Deformation Gradient Tensor
    %In a sense, this is now 2x2xjxi matrix, same position of particle on the
    %grid, except now we tacked a matrix at each point
    %so when we want to describe D matrix of a particle at x=0 and y=0, we
    %would say (:,:,1,1) or at x=1 and y=1 would say (:,:,J+1,I+1)
    Jmatrix(1,1,:,:,:) = dx_after11./dx_before;
    Jmatrix(1,2,:,:,:) = dx_after12./dy_before;
    Jmatrix(1,3,:,:,:) = dx_after13./dz_before;
    Jmatrix(2,1,:,:,:) = dy_after21./dx_before;
    Jmatrix(2,2,:,:,:) = dy_after22./dy_before;
    Jmatrix(2,3,:,:,:) = dy_after23./dz_before;
    Jmatrix(3,1,:,:,:) = dz_after31./dx_before;
    Jmatrix(3,2,:,:,:) = dz_after32./dy_before;
    Jmatrix(3,3,:,:,:) = dz_after33./dz_before;
end
end
% 1/26/26 SDB When I first updated this in Dec 2025, I permuted it to be [mxnx nlon x nlat x nht], but
% now that I've transposed to have X,Y,Z,U,V,W all be [nlat, nlon, nht], I should not permute.
%Jmatrix = permute(Jmatrix, [1, 2, 4, 3, 5]);

