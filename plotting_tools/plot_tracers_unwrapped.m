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

function hA = plot_tracers_unwrapped(tracerX, tracerY)%, guvi_flag)

%plot tracer position
hold on
delta = 2; %deg

% There are as many columns in tracerX as there are particles being traced.
numtracers = size(tracerX,2);
for i = 1:numtracers
	tracerstr = num2str(i);
	%eval(['tracer_pos_' tracerstr '(:,1) = tracerY(:,' tracerstr ')']);
	%eval(['tracer_pos_' tracerstr '(:,2) = tracerX(:,' tracerstr ')']);

	% Put in logic to check if the x dimension needs to be unwrapped, 
	% e.g., if it is longitude.
	% First check if there are any that fall outside of longitude limits
	% and enforce [0, 360].
	tracerX = tracerX + 360* (tracerX < 0) - 360*(tracerX > 360);

	% discrows will mark the index of the last point of a continuous
	% segment. discrows+1 will mark the start of the next continuous part.
	discrows = find(abs(diff(tracerX(:,i))) > 30);

	% There were no discontinuities for this tracer.
	if isempty(discrows)
	    hB = plot3(tracerX(:,i),tracerY(:,i), ones(size(tracerX(:,i))), 'k-', 'LineWidth', 1.5);%...
	    set(hB, 'Color',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers]);
	    hA = plot3(tracerX([1,end],i),tracerY([1,end],i), ones(size(tracerX([1, end],i))), 'om');%...
	    set(hA,	'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers], ...
		'MarkerSize', 6, 'LineWidth',1);
	    % Set the endpoint symbol to be a different shape.
	    hA = plot3(tracerX([end,end],i),tracerY([end,end],i), ...
	    	ones(size(tracerX([end, end],i))), '^m');%...
set(hA,	'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers], ...
	'MarkerSize', 6, 'LineWidth',1);
	else

	    % Plot the segment from element 1:discrows(1).
	    hB = plot3(tracerX(1:discrows(1),i),tracerY(1:discrows(1),i), ones(size(tracerX(1:discrows(1),i))), 'k-', 'LineWidth', 1.5);%...
	    set(hB, 'Color',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers]);
		hA = plot3(tracerX([1,discrows(1)],i),tracerY([1,discrows(1)],i), ones(size(tracerX([1, discrows(1)],i))), 'om');%...
		set(hA,	'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers], ...
			'MarkerSize', 6, 'LineWidth',1);

	    % Plot all intermediate segments.
	    for j = 2:numel(discrows)%-1
		startind = discrows(j-1)+1;
		endind = discrows(j);
		hB = plot3(tracerX(startind:endind,i), ...
		    tracerY(startind:endind,i), ...
		    ones(size(tracerX(startind:endind,i))), 'k-', 'LineWidth', 1.5);%...
		set(hB, 'Color',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers]);
		hA = plot3(tracerX([startind,endind],i), ...
		    tracerY([startind,endind],i), ...
		    ones(size(tracerX([startind, endind],i))), 'om');%...
		set(hA,	'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers], ...
			'MarkerSize', 6, 'LineWidth',1);
	    end

	    % Plot the last segment.
	    startind = discrows(end)+1;
	    hB = plot3(tracerX(startind:end,i),tracerY(startind:end,i), ones(size(tracerX(startind:end,i))), 'k-', 'LineWidth', 1.5);%...
	    set(hB, 'Color',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers]);
	    hA = plot3(tracerX([startind, startind],i),tracerY([startind,startind],i), ones(size(tracerX([startind, startind],i))), 'om');%...
		set(hA,	'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers], ...
			'MarkerSize', 6, 'LineWidth',1);

	    % Plot the endpoint with a different symbol.
	    hA = plot3(tracerX([end, end],i),tracerY([end,end],i), ones(size(tracerX([end, end],i))), '^m');%...
		set(hA,	'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers], ...
			'MarkerSize', 6, 'LineWidth',1);

	end
end


%hB = plotm(tracer_pos_2(1:1:end,2),tracer_pos_2(1:1:end,1), 'o-r');
%set(hB, 'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0], ...
%	'MarkerSize', 10, 'LineWidth',2);
%str1 = {'A_0'};
%strf = {'A_f'};
%if guvi_flag
%	Atextclrstr = 'm';
%	Btextclrstr = 'r';
%else
%	Atextclrstr = 'k';
%	Btextclrstr = 'k';
%end
%textm(tracer_pos_1(1,2)+delta*0,tracer_pos_1(1,1)-delta*6.5,str1, ...
%	'FontWeight', 'bold', 'FontSize', 24, 'Color', Atextclrstr);
%textm(tracer_pos_1(end,2)+delta*0,tracer_pos_1(end,1)+delta*2.5,strf, ...
%	'FontWeight', 'bold', 'FontSize', 24, 'Color', Atextclrstr);
%str1 = {'B_0'};
%strf = {'B_f'};
%textm(tracer_pos_2(1,2)+delta*0,tracer_pos_2(1,1)-delta*6.5,str1, ...
%	'FontWeight', 'bold', 'FontSize', 24, 'Color', Btextclrstr);
%textm(tracer_pos_2(end,2)+delta*0,tracer_pos_2(end,1)+delta*2.5,strf, ...
%	'FontWeight', 'bold', 'FontSize', 24, 'Color', Btextclrstr);
%end
