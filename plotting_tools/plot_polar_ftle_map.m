function h = plot_polar_ftle_map(modelstr, t0, lat, lon, ht, FTLE_Value, regionstr, integration_hours, J_in_Cartesian_flag)
% plot polar FTLE map 
% Input modelstr should be 'navegem', 'SD-WACCM', or 'WACCM-X_DART'
% regionstr should be 'North' or 'South'
% S. Datta-Barua
% 24 Mar 2020 Cropped from waccm_lcs_script.m.
% 7 Oct 2025 Plot map axes without FTLE colors.

%% polar plot // uncomment when need to plot the polar plot
%glat = -88.75:2.5:88.75;
%glon = -180:2.5:180;
%figure
% subplot(1,2,1)

fontsz = 8;
utc = datevec(t0);
t_local = 0;%180-(utc(4)+utc(5)/60+utc(6)/3600)/15;

switch regionstr
    case {'North Pole'}
        origin = [90, t_local];
    case 'South Pole'
        origin = [-90, t_local];
end

h = axesm('MapProjection','ortho','origin',origin,'meridianlabel','on','parallellabel','on','grid','on','mlabelparallel',0,'labelrotation','on','fontsize',fontsz,'fontweight','bold','labelformat','none'); % change the font size
[gridglat, gridglon] = ndgrid(lat, lon);
% [c,h_p] = contourm(gridglat, gridglon, B','LevelStep',5e-6,'LineWidth',1.5); %for 3hour integration
if ~isempty(FTLE_Value)
pcolorm(gridglat, gridglon, FTLE_Value)
end
framem

% From R2020b onward, coast.mat is not in the Mapping Toolbox.  Use coastlines.mat.
try
	load coast
	DDD=plotm(lat,long,'w');
catch
	load coastlines
	DDD=plotm(coastlat,coastlon,'w');
end

hc = colorbar;
ylabel(hc, 'FTLE Value [1/s]');
% Check for backward (negative) or forward (positive) values of FTLE.
%if any(FTLE_Value < 0)
%	set(gca, 'CLim', [-4e5, 0]);
%else
if J_in_Cartesian_flag == 2
	set(gca, 'CLim', [4e-5 6e-5])% 
elseif J_in_Cartesian_flag == 1
	set(gca, 'CLim', [1e-4 1.5e-4])% 
else
	set(gca, 'CLim', [0.5e-5 4e-5]);%[1e-4 1.5e-4]);
end

titleofplot = [modelstr ' FTLE map \newline ', ...
    datestr(t0, 'mm/dd/yyyy hh:MM') ' UT, ', ...
    num2str(ht * ((ht>1e4)/1e3)) ' km, ', '\tau = ' num2str(integration_hours), ' hrs'];
title(titleofplot,'FontSize',fontsz);%12)
    

