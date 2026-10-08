function w = omega2w(pressurelevels, nt, nlon, nlat, T, omegafull)
% function w = omega2w(pressurelevels, nt, nlon, nlat, T, omegafull)
% computes upward velocity in m/s from pressure rate omega.

% Notes from Nick Pedatella's email 3/16/23
% S. Datta-Barua
% 23 May 2025 Checking the units.

% # get vertical velocity
rgas = 287.058;
g = 9.80665;
ptemp = repmat(pressurelevels, 1, nt, nlon, nlat);
rho = shiftdim(ptemp, 2) .* 100 ./ (rgas*T); %# lev hPa -> Pa

% Nick's output was in cm/s, but we want m/s.
%w = -omegafull./(rho*g)*100; %# m/s -> cm/s
w = -omegafull./(rho*g); % # m/s
% where lev is the pressure level and om is omega.

