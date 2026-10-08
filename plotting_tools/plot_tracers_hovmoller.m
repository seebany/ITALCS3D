% function h = plot_tracers_unwrapped(tracerX, tracerY)
% plots tracers onto an existing map projection.  
%
% S. Datta-Barua
% 24 Feb 2020
% 24 Mar 2020 Include integration time tau on the plot title.
% 3 Apr 2020 Optional input args tracerX and tracerY will have tracers
% plotted also.
% 5 Aug 2020 Making this map to more closely match the Meier et al., 2011 paper.
% Requires Mapping Toolbox.
% 12 Jan 2021 Putting in plotting of GUVI data.
% 26 Jul 2021 Removing GUVI plotting again. Plotting multiple tracers.
% 16 May 2024 To plot on other 2D slices, make this version not use mapping toolbox.
% 9 Feb 2026 Assuming tracerX is longitude, work to separate continuous segments.
%   tracerX = [ntimes x numtracers], x coordinate only.
% tracerY = [ntimes x numtracers], the vertical coordinate.
% tracerZ = [ntimes x numtracers], the coordinate to be assigned color.
% 25 Aug 2026 Adapting to make a plot of tracers in lon-time axes (hovmoller). For this tracerX is the longitudes, tracerY is datetime. tracerZ is not used.
% 22 Sept 2026 Now setting tracerZ to be the color traced on the line (corresponding to height or latitude as desired).

function hA = plot_tracers_hovmoller(tracerX, tracerY, tracerZ)%, guvi_flag)

%plot tracer position
hold on
delta = 2; %deg
set(gca, 'YDir', 'reverse');

% There are as many columns in tracerX as there are particles being traced.
numtracers = size(tracerX,2);
	
% Put in logic to check if the x dimension needs to be unwrapped, 
% e.g., if it is longitude.
% First check if there are any that fall outside of longitude limits
% and enforce [0, 360].
tracerX = tracerX + 360* (tracerX < 0) - 360*(tracerX > 360);

for i = 1:numtracers
	%tracerstr = num2str(i);


	% discrows will mark the index of the last point of a continuous
	% segment. discrows+1 will mark the start of the next continuous part.
	discrows = find(abs(diff(tracerX(:,i))) > 30);

	% There were no discontinuities for this tracer.
	if isempty(discrows)
	    % Plot a line for the tracer path.
	    hB = plot(tracerX(:,i),tracerY(:,i), 'k-', 'LineWidth', 1.5);%...
	    %ones(size(tracerX(:,i))), 
	    set(hB, 'Color',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers]);
	    % Plot a starting marker.
	    hA = plot(tracerX([1,1],i),tracerY([1,1],i), 'ok');
	    %ones(size(tracerX([1, end],i))), 'om');%...
	    set(hA,	'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers], ...
		'MarkerSize', 6, 'LineWidth',1);
	    % Set the endpoint symbol to be a different shape.
	    hA = plot(tracerX([end,end],i),tracerY([end,end],i), '^k');%...
%	    	ones(size(tracerX([end, end],i))), '^m');%...
	    set(hA,	'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers], ...
		'MarkerSize', 6, 'LineWidth',1);

	    % Plot the line with a color map to show the height.
	    clear x y z c
	    x = tracerX([2:end-1], i)';
	    y = tracerY([2:end-1], i)';
	    z = zeros(size(tracerX([2:end-1],i)))';
	    c = tracerZ([2:end-1], i)';
	    hA = surf([x;x], [y;y], [z;z], [c;c], ...
	    'facecolor', 'none', 'edgecolor', 'interp', 'linewidth', 2);
	    %map = colormap;
	    %set(hA, 'MarkerEdgeColor', tracerZ([2:end], i);
	else % There are discontinuities due to crossing a boundary.

	    % Plot the segment from element 1:discrows(1).
	    hB = plot(tracerX(1:discrows(1),i),tracerY(1:discrows(1),i),...
	   'k-', 'LineWidth', 1.5);
	  % ones(size(tracerX(1:discrows(1),i))), 'k-', 'LineWidth', 1.5);%...
	    set(hB, 'Color',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers]);
	    hA = plot(tracerX([1,discrows(1)],i),tracerY([1,discrows(1)],i), 'om');%ones(size(tracerX([1, discrows(1)],i))), 'om');%...
	    set(hA,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers], ...
		'MarkerSize', 6, 'LineWidth',1);
	    % Plot the line with a color map to show the height.
	    clear x y z c
	    x = tracerX([1:discrows(1)], i)';
	    y = tracerY([1:discrows(1)], i)';
	    z = zeros(size(tracerX([1:discrows(1)],i)))';
	    c = tracerZ([1:discrows(1)], i)';
	    %keyboard
	    hA = surface([x;x], [y;y], [z;z], [c;c], ...
	    'facecolor', 'none', 'edgecolor', 'interp', 'linewidth', 2);

	    % Plot all intermediate segments.
	    for j = 2:numel(discrows)%-1
		startind = discrows(j-1)+1;
		endind = discrows(j);
		hB = plot(tracerX(startind:endind,i), ...
		    tracerY(startind:endind,i), ...
		    'k-', 'LineWidth', 1.5);%...
		    %ones(size(tracerX(startind:endind,i))), 
		set(hB, 'Color',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers]);
	    % Plot the line with a color map to show the height.
	    clear x y z c
	    x = tracerX(startind:endind, i)';
	    y = tracerY(startind:endind, i)';
	    z = zeros(size(tracerX(startind:endind,i)))';
	    c = tracerZ(startind:endind, i)';
	    %keyboard
	    surface([x;x], [y;y], [z;z], [c;c], ...
	    'facecolor', 'none', 'edgecolor', 'interp', 'linewidth', 2);

		hA = plot(tracerX([startind,endind],i), ...
		    tracerY([startind,endind],i), 'om');%...
		    %ones(size(tracerX([startind, endind],i))), 'om');%...
		set(hA,	'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers], ...
			'MarkerSize', 6, 'LineWidth',1);
	    end

	    % Plot the last segment.
	    startind = discrows(end)+1;
	    hB = plot(tracerX(startind:end,i),tracerY(startind:end,i), ... ones(size(tracerX(startind:end,i))), 
	    'k-', 'LineWidth', 1.5);%...
	    set(hB, 'Color',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers]);
	    % Plot the line with a color map to show the height.
	    clear x y z c
	    x = tracerX(startind:end, i)';
	    y = tracerY(startind:end, i)';
	    z = zeros(size(tracerX(startind:end,i)))';
	    c = tracerZ(startind:end, i)';
	    %keyboard
	    surface([x;x], [y;y], [z;z], [c;c], ...
	    'facecolor', 'none', 'edgecolor', 'interp', 'linewidth', 2);
	    hA = plot(tracerX([startind, startind],i),tracerY([startind,startind],i), 'om');%ones(size(tracerX([startind, startind],i))), 'om');%...
		set(hA,	'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers], ...
			'MarkerSize', 6, 'LineWidth',1);

	    % Plot the endpoint with a different symbol.
	    hA = plot(tracerX([end, end],i),tracerY([end,end],i), '^m');%ones(size(tracerX([end, end],i))), '^m');%...
		set(hA,	'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers], ...
			'MarkerSize', 6, 'LineWidth',1);

	end
end
ytickformat('MM-dd');
axis tight
