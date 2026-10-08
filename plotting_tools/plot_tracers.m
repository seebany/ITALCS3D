% function h = plot_tracers(tracerX, tracerY)
% plots start and endpoint tracers connected by a line
% onto an existing map projection.  
%   tracerX = [ntimes x numtracers], x coordinate only.
%
% S. Datta-Barua
% 7 Oct 2026 Preparing for publication to github.

function hA = plot_tracers(tracerX, tracerY)%, guvi_flag)

%plot tracer position
hold on
delta = 2; %deg

% There are as many columns in tracerX as there are particles being traced.
numtracers = size(tracerX,2);
for i = 1:numtracers
	tracerstr = num2str(i);
	%eval(['tracer_pos_' tracerstr '(:,1) = tracerY(:,' tracerstr ')']);
	%eval(['tracer_pos_' tracerstr '(:,2) = tracerX(:,' tracerstr ')']);
end

for i = 1:numtracers
hB = plot3(tracerX(:,i),tracerY(:,i), ones(size(tracerX(:,i))), 'k-', 'LineWidth', 1.5);%...
set(hB, 'Color', [(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers]);
% Plot the startpoints with one symbol.
hA = plot3(tracerX([1,1],i),tracerY([1,1],i), ones(size(tracerX([1, 1],i))), 'om');%...
set(hA,	'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers], ...
	'MarkerSize', 6, 'LineWidth',1);
	% Plot the endpoints with a different symbol.

	hA = plot3(tracerX([end,end],i),tracerY([end,end],i), ones(size(tracerX([end, end],i))), '^m');%...
set(hA,	'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers], ...
	'MarkerSize', 6, 'LineWidth',1);
end
