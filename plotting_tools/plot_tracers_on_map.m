% function h = plot_tracers_on_map(tracerlat, tracerlon)
% plots tracers onto an existing map projection.  
%
% S. Datta-Barua
% 24 Feb 2020
% 24 Mar 2020 Include integration time tau on the plot title.
% 3 Apr 2020 Optional input args tracerlat and tracerlon will have tracers
% plotted also.
% 5 Aug 2020 Making this map to more closely match the Meier et al., 2011 paper.
% Requires Mapping Toolbox.
% 12 Jan 2021 Putting in plotting of GUVI data.
% 26 Jul 2021 Removing GUVI plotting again. Plotting multiple tracers.

function hA = plot_tracers_on_map(tracerlat, tracerlon)%, guvi_flag)

%plot tracer position
hold on
delta = 2; %deg

% There are as many columns in tracerlat as there are particles being traced.
numtracers = size(tracerlat,2);
for i = 1:numtracers
	tracerstr = num2str(i);
	%eval(['tracer_pos_' tracerstr '(:,1) = tracerlon(:,' tracerstr ')']);
	%eval(['tracer_pos_' tracerstr '(:,2) = tracerlat(:,' tracerstr ')']);
end

for i = 1:numtracers
hB = plotm(tracerlat(:,i),tracerlon(:,i), 'k-', 'LineWidth', 1.5);%...
set(hB, 'Color',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers]);
hA = plotm(tracerlat([1,1],i),tracerlon([1,1],i), 'om');%...
set(hA,	'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers], ...
	'MarkerSize', 6, 'LineWidth',1);
	hA = plotm(tracerlat([end,end],i),tracerlon([end,end],i), '^m');%...
set(hA,	'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[(i-1)/numtracers, (i-1)/numtracers, (i-1)/numtracers], ...
	'MarkerSize', 6, 'LineWidth',1);
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
