% Function for linearize the velocity
% By Ningchao Wang 2019.2.3

function [ U_l, V_l ] = vel_linear( X, Y, U, V)
%Notes on the inputs of this function:
%X is the horizontal coordinate of the point on the grid
%Y is the vertical coordinate of the point on the grid
%U is the vertical velocity at the point (X,Y) on the grid
%V is the horizontal velocity at the point (X,Y) on the grid

%Notes on the outputs of this function:
%U_l is the horizontal linear velocity at the point (X,Y) on the grid
%V_l is the vertical linear velocity at the point (X,Y) on the grid
%"linear velocity" means that U_l and V_l are the closest linear,
%divergence free velocity field to the inputs, U and V respectively, in the
%C0 norm over the grid

%V_l_U = []; %linear zonal velocity | ZJB - value unused in code so
%commented out
%V_l_V = []; %linear meriadional velocity | ZJB - value unused in code so commented
%out

%Nomenclature tendancy is as follows: N11 is the pre-averaged numerator for
%location (1,1) (top right) of the A matrix in v_L = Au + b where v_L is the 
%boundary layer linear velocity. n11 is then the post-average numerator,
%and d11 is the post-average denmoninator. The numerator is then divided by
%the denominator to give the value of A(1,1). This process is repeated for
%each of the positions in the A matrix.

% Changed average() function to mean(mean()) --> Speeds up computation
% Marta CV. 09/28/2021

N11 = X.*U-Y.*V;
n11 = mean(mean(N11));
D11 = X.^2 + Y.^2;
d11 = mean(mean(D11));
a11 = n11/d11;
N12 = Y.*U;
n12 = mean(mean(N12));
D12 = Y.^2;
d12 = mean(mean(D12));
a12 = n12/d12;
N21 = X.*V; %Changed from X.*U on 3/3/20 by ZB
n21 = mean(mean(N21));
D21 = X.^2;
d21 = mean(mean(D21));
a21 = n21/d21;
N22 = Y.*V-X.*U;
n22 = mean(mean(N22));
D22 = X.^2 + Y.^2;
d22 = mean(mean(D22));
a22 = n22/d22;

AA = [a11,a12;a21,a22]; %creates the A matrix

% U_l = zeros(v1,v2); %added by Zack Bonson to increase speed of loop by preallocating
% V_l = zeros(v1,v2); %added by Zack Bonson to increase speed of loop by preallocating

%The following loop fills in each location of the linear velocity output
%commented out due to note below
% for i = 1:v1
%     for j = 1:v2
%         vel_L = AA*[U(i,j);V(i,j)]+[average(U);average(V)];
%         U_l(i,j) = vel_L(1,1);
%         V_l(i,j) = vel_L(2,1);
%     end
% end

%proposed alternate loop to see if there is any change to the error - ZJB
%I have changed the U and V in the "vel_L = AA*..." line to X and Y as that
%appears to be more consistent with the source formula.

U_ave = mean(mean(U));
V_ave = mean(mean(V));
% Computing V and U means outside the 'for' loop reduces the computation time.
% Marta CV. 09/08/2021


[m n] = size(X);
[V_aux] = solve_linear_eq(AA,X,Y,U_ave,V_ave);
U_l = reshape(V_aux(1,:),m,n);
V_l = reshape(V_aux(2,:),m,n);
% Alternative code to solve the linear equation. Speeds up computation time. 
% This equation is solved calling 'SolvingLinearEq.m' function.
% Substitutes the double 'for' loop commented out just below.
% Marta CV. 10/06/2021


% [v1 v2] = size(U); %gives us the number of rows (v1) and the number of columns (v2) for linear velocity matrix
% for i = 1:v1
%     for j = 1:v2
%         vel_L = AA*[X(i,j);Y(i,j)]+[U_ave;V_ave];
%         U_l(i,j) = vel_L(1,1);
%         V_l(i,j) = vel_L(2,1);
%     end
% end
end





