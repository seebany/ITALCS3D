% This function computes the Jacobian matrix (deformation matrix) in a 2D 
% or 3D domain in either euclidean or non-euclidean space.

function [Jmatrix] = JacobianGradient(X0,Y0,Z0,X,Y,Z,Xres,Yres,Zres,dimension_flag,Euclidean_flag)

% Inputs:
%       XO/Y0/Z0: initial x,y,z positions for every gridpoint and timestamp.
%       X/Y/Z: final x,y,z positions for every gridpoint and timestamp.
%       dimension_flag: 2D or 3D computation (0,1).
%       Euclidean_flag: Euclidean or Non-Euclidean computation (1,0)
% Outputs:
%       Jmatrix: Jacobian matrix (3x3xjxixk). Each particle has a unique position in the matrix, where an individual 3x3 matrix is saved.
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
else if dimension_flag==0
    [dx_before,~]=gradient(X0);
    [~,dy_before]=gradient(Y0);
    
%     for i=1:Yres
%         for j=1:Xres-1
%             if abs(X(i,j)-X(i,j+1)) > 300
%                 X(i,j+1)=X(i,j+1)+360;
%             end
%         end
%     end
    
    %  update final matrix
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
    
        %restract the maximum distance
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
    else
        [dx_before,~,~] = gradient(X0);
        [~,dy_before,~] = gradient(Y0);
        [~,~,dz_before] = gradient(Z0);
        
    % 3D
%     for i = 1:Zres
%     %  update final matrix
%     % for x direction, the latitude boundary
%     P = [X(1,round((Xres+1)/2):end,i),X(1,1:round((Xres-1)/2),i);X(:,:,i);X(end,round((Xres+1)/2):end,i),X(end,1:round((Xres-1)/2),i)];
%     % for x direction, the longitude boundary
%     XP(:,:,i) = [P(:,end),P,P(:,1)];
%     % for y direction, the latitude boundary
%     Q = [Y(1,round((Xres+1)/2):end,i),Y(1,1:round((Xres-1)/2),i);Y(:,:,i);Y(end,round((Xres+1)/2):end,i),Y(end,1:round((Xres-1)/2),i)];
%     % for y direction, the longitude boundary
%     YQ(:,:,i) = [Q(:,end),Q,Q(:,1)];
%     end
%     
%     X=XP;
%     Y=YQ;

% 3D Vertical
for i = 1:Yres
    %  update final matrix
    % for x direction, the latitude boundary
    P = [X(i,round((Xres+1)/2):end,1),X(i,1:round((Xres-1)/2),1);permute(X(i,:,:),[3,2,1]);X(i,round((Xres+1)/2):end,end),X(i,1:round((Xres-1)/2),end)];
    % for x direction, the longitude boundary
    XP(i,:,:) = [P(:,end),P,P(:,1)];
    % for y direction, the latitude boundary
    Q = [Z(i,round((Xres+1)/2):end,1),Z(i,1:round((Xres-1)/2),1);permute(Z(i,:,:),[3,2,1]);Z(i,round((Xres+1)/2):end,end),Z(i,1:round((Xres-1)/2),end)];
    % for y direction, the longitude boundary
    YQ(i,:,:) = [Q(:,end),Q,Q(:,1)];
end
    
X=permute(XP,[1,3,2]);
Z=permute(YQ,[1,3,2]);

    [dx_after11,dx_after12,dx_after13]=gradient(X);
    [dy_after21,dy_after22,dy_after23]=gradient(Y);
    [dz_after31,dz_after32,dz_after33]=gradient(Z);
    
    % pick the inner part of the gradient matrix
    
% 3D
%     dx_after11 = dx_after11(2:end-1,2:end-1,:);
%     dx_after12 = dx_after12(2:end-1,2:end-1,:);
%     dx_after13 = dx_after13(2:end-1,2:end-1,:);
%     
%     dy_after21 = dy_after21(2:end-1,2:end-1,:);
%     dy_after22 = dy_after22(2:end-1,2:end-1,:);
%     dy_after23 = dy_after23(2:end-1,2:end-1,:);
%     
% %     dz_after31 = dz_after31(2:end-1,2:end-1,:);
% %     dz_after32 = dz_after32(2:end-1,2:end-1,:);
% %     dz_after33 = dz_after33(2:end-1,2:end-1,:);


% 3D Vertical
    dx_after11 = dx_after11(:,2:end-1,2:end-1);
    dx_after12 = dx_after12(:,2:end-1,2:end-1);
    dx_after13 = dx_after13(:,2:end-1,2:end-1);
    
    dz_after31 = dz_after31(:,2:end-1,2:end-1);
    dz_after32 = dz_after32(:,2:end-1,2:end-1);
    dz_after33 = dz_after33(:,2:end-1,2:end-1);

%     dx_after11 = dx_after11(:,2:end-1);
%     dx_after12 = dx_after12(:,2:end-1);
%     dy_after21 = dy_after21(:,2:end-1);
%     dy_after22 = dy_after22(:,2:end-1);
    %Defines the deformation gradient tensor (or the D matrix)
    
        %restrict the maximum distance
keyboard
for k = 1:Zres,
    for i = 1:Yres,
        for j = 1:Xres,
            if (abs(dx_after11(i,j,k)) > 90)
                dx_after11(i,j,k)=(sign(dx_after11(i,j,k)))*(180-(abs(dx_after11(i,j,k))));
            else
                dx_after11(i,j,k)=dx_after11(i,j,k);
            end
        end
    end
    
    for i = 1:Yres,
        for j = 1:Xres,
            if (abs(dx_after12(i,j,k)) > 90)
                dx_after12(i,j,k)=(sign(dx_after12(i,j,k)))*(180-(abs(dx_after12(i,j,k))));
            else
                dx_after12(i,j,k)=dx_after12(i,j,k);
            end
        end
    end
    
    for i = 1:Yres,
        for j = 1:Xres,
            if (abs(dx_after13(i,j,k)) > 90)
                dx_after13(i,j,k)=(sign(dx_after13(i,j,k)))*(180-(abs(dx_after13(i,j,k))));
            else
                dx_after13(i,j,k)=dx_after13(i,j,k);
            end
        end
    end
                
% 3D
%     for i = 1:Yres,
%         for j = 1:Xres,
%             if (abs(dy_after21(i,j,k)) > 90)
%                 dy_after21(i,j,k)=(sign(dy_after21(i,j,k)))*(180-(abs(dy_after21(i,j,k))));
%             else
%                 dy_after21(i,j,k)=dy_after21(i,j,k);
%             end
%         end
%     end
% 
%     
%     for i = 1:Yres,
%         for j = 1:Xres,
%             if (abs(dy_after22(i,j,k)) > 90)
%                 dy_after22(i,j,k)=(sign(dy_after22(i,j,k)))*(180-(abs(dy_after22(i,j,k))));
%             else
%                 dy_after22(i,j,k)=dy_after22(i,j,k);
%             end
%         end
%     end
% end

% 3D Vertical
for i = 1:Yres,
        for j = 1:Xres,
            if (abs(dz_after31(i,j,k)) > 90)
                dz_after31(i,j,k)=(sign(dz_after31(i,j,k)))*(90-(abs(dz_after31(i,j,k))));
            else
                dz_after31(i,j,k)=dz_after31(i,j,k);
            end
        end
    end

    
    for i = 1:Yres,
        for j = 1:Xres,
            if (abs(dz_after32(i,j,k)) > 90)
                dz_after32(i,j,k)=(sign(dz_after32(i,j,k)))*(90-(abs(dz_after32(i,j,k))));
            else
                dz_after32(i,j,k)=dz_after32(i,j,k);
            end
        end
    end
end
    
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
