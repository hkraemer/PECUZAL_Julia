clear, clc

% plotting the results from script
% `analyze_standard_ecoscillator_try2.jl.jl`


% h = plot(1:10,1:10,1:10,2:11,1:10,3:12,1:10,4:13);
% c = get(h,'Color');
% 
% c1 = c{1};
% c2 = c{4};
% c3 = c{3};
% c4 = c{2};

% legend fontsize
fsl = 30;
% label fontsize
fsla = 30;
% title fontsize
fst = 40;

% columns indicate:
% ENTR, LAM, RTE, TRANS

% Experiment 1
% Initialize data points
TDE = [.9997 .9997 .8708 .7307];
GA = [.952 .9976 .9778 .9744 ];
MDOP = [0.938 .9955 .9959 .9509];
PECUZAL = [.9557 .9997 .9822 .9869];
P = [TDE; GA; MDOP; PECUZAL];

figure('Units','normalized','Position',[.3 .3 .8 .8])
% Spider plot
%'Color', [c1; c2; c3; c4],...
spider_plot(P,...
    'AxesLabels', {'ENTR', 'LAM', 'RTE', 'TRANS'},...
    'AxesLimits', [.85, .99, .8, .7; 1, 1, 1, 1],...
    'FillOption', 'on',...
    'AxesPrecision', 3,...
    'LabelFontSize', fsla,...
    'FillTransparency', 0.2);
title('Experiment I','FontSize',fst,'FontName', 'Times New Roman')

legend('standard TDE', 'Garcia & Almeida', 'MDOP', 'PECUZAL', 'Location', 'southoutside','Fontsize',fsl);

% Experiment 2
% Initialize data points

% columns indicate:
% ENTR, LAM, RTE, TRANS
TDE = [.9957 .9996 .8792 .9561];
GA = [.935 .9994 .9212 .9312 ];
MDOP = [0.8776 .999 .9567 .9421];
PECUZAL = [.9588 .9996 .8892 .9227];
P = [TDE; GA; MDOP; PECUZAL];

figure('Units','normalized','Position',[.3 .3 .8 .8])
% Spider plot
spider_plot(P,...
    'AxesLabels', {'ENTR', 'LAM', 'RTE', 'TRANS'},...
    'AxesLimits', [.85, .99, .8, .7; 1, 1, 1, 1],...
    'FillOption', 'on',...
    'AxesPrecision', 3,...
    'LabelFontSize', fsla,...
    'FillTransparency', 0.2);
title('Experiment II','FontSize',fst,'FontName', 'Times New Roman')

legend('standard TDE', 'Garcia & Almeida', 'MDOP', 'PECUZAL', 'Location', 'southoutside','Fontsize',fsl);
