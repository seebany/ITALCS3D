function [V_aux] = solve_linear_eq(AA,X,Y,U_ave,V_ave)

% Gets U and V velocicities after linearization. It solves the following
% matrix equation:
% V = Ax + b
% where:
% A: (2x2) coefficient matrix.
% x: position vector x = [X;Y].
% b: average velocities (U_ave, V_ave).
%
% INPUTS:
% AA: (2x2) coefficient matrix.
% X, Y: position vectors for 2D.
% U_ave, V_ave: average velocities.
%
% OUTPUTS:
% V_aux: matrix (2 x number of gridpoints).
%
%  Marta CV. 10/06/2021

[m n] = size(X);

X = reshape(X,1,m*n);
Y = reshape(Y,1,m*n);

XY = [X;Y];

V_aux(1,:) = AA(1,1)*XY(1,:) + AA(1,2)*XY(2,:) + U_ave;
V_aux(2,:) = AA(2,1)*XY(1,:) + AA(2,2)*XY(2,:) + V_ave;


return

