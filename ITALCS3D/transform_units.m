% function Jtilde = transform_units(X03d, Y03d, Z03d, Jmatrix)
% takes a 5d array of [3 x 3 x nlat x nlon x nht]
% of the Jacobian of the flow map as expressed in
% changes in (lat, lon, ht) and outputs the 
% Jacobian in units of (m, m, m).
%
% Seebany Datta-Barua
% 30 Jan 2026

function Jtilde = transform_units(X03d, Y03d, Z03d, Jmatrix)

% Initialize output to match dimensions of input.
Jtilde = zeros(size(Jmatrix));
[Yres, Xres, Zres] = size(X03d);
tic
% Need transformation matrices.
% The 3x3 matrix that transforms from llh to r, \theta, \phi.
rTl = [0 0 1; -1 0 0; 0 1 0];
% Since rTl is orthogonal, the inverse is the transpose.
lTr = rTl';

% Write llh points as (r,theta, phi) points.
Re = 6378e3; % Earth radius in m.
r = Re + Z03d;
theta = 90 - Y03d;
phi = X03d;

% Adjust theta so that matrix inversion won't be singular.
polerows = find(theta == 0);
theta(polerows) = 0.05;
clear polerows
polerows = find(theta == 180);
theta(polerows) = 179.95;
clear polerows

% Construct each 3x3 matrix that transforms deformations in
% (r, theta, phi) to Cartesian (x, y, z).
% This matrix is \begin{matrix}
%	\partial{x}/\partial{r}, \partial{x}/\partial{\theta}, \partial{x}{\partial \phi};
%	\partial{y}/\partial{r}, ...
%
xTr = zeros(size(Jmatrix));
xTr(1,1,:,:,:) = sind(theta).*cosd(phi);
xTr(1,2,:,:,:) = r.*cosd(theta).*cosd(phi);
xTr(1,3,:,:,:) = -r.*sind(theta).*sind(phi);
xTr(2,1,:,:,:) = sind(theta).*sind(phi);
xTr(2,2,:,:,:) = r.*cosd(theta).*sind(phi);
xTr(2,3,:,:,:) = r.*sind(theta).*cosd(phi);
xTr(3,1,:,:,:) = cosd(theta);
xTr(3,2,:,:,:) = -r.*sind(theta);
%xTr(3,3,:,:,:) = 0;

% Jtilde is constructed for each gridpoint by pre- and
% post-multiplying by these transformation matrices.
for i = 1:Yres
    for j = 1:Xres
        for k = 1:Zres
		% Warnings about singular matrices, 
		% presumably at lat = +/-90 deg.
                rTx(:,:,i,j,k) = inv(xTr(:,:,i,j,k));
		Jtilde(:,:,i,j,k) = xTr(:,:,i,j,k)*rTl...
			* Jmatrix(:,:,i,j,k) ...
			* lTr * rTx(:,:,i,j,k); 
        end
    end
end
toc
