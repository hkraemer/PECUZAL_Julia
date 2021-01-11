clear, clc

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
% JRRF, MFNN, L, DET, ENTR, RTE

% Rössler uni
% Initialize data points
TDE = [.71 .85 .9 .9999 .95 .97];
GA = [.48 -.64 .73 .9992 .68 .93];
MDOP = [0.64 .79 .88 .9999 .82 .93];
PECUZAL = [.73 1 1 .9999 .92 .96];
P = [TDE; GA; MDOP; PECUZAL];

figure('Units','normalized','Position',[.3 .3 .8 .8])
% Spider plot
%'Color', [c1; c2; c3; c4],...
spider_plot(P,...
    'AxesLabels', {'JRRF', 'MFNN', '\DeltaL', 'DET', 'ENTR', 'RTE'},...
    'AxesLimits', [.4, -.64, .6, .6, .6, .6; 1, 1, 1, 1, 1, 1],...
    'FillOption', 'on',...
    'LabelFontSize', fsla,...
    'FillTransparency', 0.2);
title('Rössler (y-component)','FontSize',fst,'FontName', 'Times New Roman')

legend('standard TDE', 'Garcia & Almeida', 'MDOP', 'PECUZAL', 'Location', 'southoutside','Fontsize',fsl);

% Rössler multi
% Initialize data points
TDE = [.71 -.98 0.8 .9999 .95 .97];
GA = [.83 .52 1 .999 .84 .94];
MDOP = [0.86 .9 .8 .9998 .91 .98];
PECUZAL = [.87 1 .94 .9999 .92 .99];
P = [TDE; GA; MDOP; PECUZAL];

figure('Units','normalized','Position',[.3 .3 .8 .8])
% Spider plot
spider_plot(P,...
    'AxesLabels', {'JRRF', 'MFNN', '\DeltaL', 'DET', 'ENTR', 'RTE'},...
    'AxesLimits', [.4, -1, .6, .6, .6, .6; 1, 1, 1, 1, 1, 1],...
    'FillOption', 'on',...
    'LabelFontSize', fsla,...
    'FillTransparency', 0.2);
title('Rössler (x- & y-component)','FontSize',fst,'FontName', 'Times New Roman')

legend('standard TDE', 'Garcia & Almeida', 'MDOP', 'PECUZAL', 'Location', 'southoutside','Fontsize',fsl);


% Duffing uni
% Initialize data points
TDE = [.82 .94 .95 .9961 .56 .97];
GA = [.8 .71 .96 0.9964 .64 .98];
MDOP = [.83 .95 .95 0.9959 .58 .99];
PECUZAL = [.84 1 1 0.9958 .55 .98];
P = [TDE; GA; MDOP; PECUZAL];

figure('Units','normalized','Position',[.3 .3 .8 .8])
% Spider plot
spider_plot(P,...
    'AxesLabels', {'JRRF', 'MFNN', '\DeltaL', 'DET', 'ENTR', 'RTE'},...
    'AxesLimits', [.5, .7, .7, .7, .4, .7; 1, 1, 1, 1, 1, 1],...
    'FillOption', 'on',...
    'LabelFontSize', fsla,...
    'FillTransparency', 0.2);
title('Duffing (x-component)','FontSize',fst,'FontName', 'Times New Roman')

legend('standard TDE', 'Garcia & Almeida', 'MDOP', 'PECUZAL', 'Location', 'southoutside','Fontsize',fsl);

% Duffing multi
% Initialize data points
TDE = [.82 .73 .94 .9961 .56 .97];
GA = [.8 .76 .88 .9917 0.51 .84];
MDOP = [.83 1 .81 0.9958 0.49 .98];
PECUZAL = [.84 .93 1 0.9958 0.56 .98];
P = [TDE; GA; MDOP; PECUZAL];

figure('Units','normalized','Position',[.3 .3 .8 .8])
% Spider plot
spider_plot(P,...
    'AxesLabels', {'JRRF', 'MFNN', '\DeltaL', 'DET', 'ENTR', 'RTE'},...
    'AxesLimits', [.5, .7, .7, .7, .4, .7; 1, 1, 1, 1, 1, 1],...
    'FillOption', 'on',...
    'LabelFontSize', fsla,...
    'FillTransparency', 0.2);
title('Duffing (x- & y-component)','FontSize',fst,'FontName', 'Times New Roman')

legend('standard TDE', 'Garcia & Almeida', 'MDOP', 'PECUZAL', 'Location', 'southoutside','Fontsize',fsl);