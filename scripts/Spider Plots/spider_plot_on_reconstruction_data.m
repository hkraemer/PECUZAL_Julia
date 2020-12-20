clear, clc

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
TDE = [.78 .9 .97 .9999 .93 .95];
GA = [.58 -1.05 .82 .9992 .710 .81];
MDOP = [0.75 .9 1 .9999 .9 .9];
PECUZAL = [.78 1 1 .9999 .92 .95];
P = [TDE; GA; MDOP; PECUZAL];

figure('Units','normalized','Position',[.3 .3 .8 .8])
% Spider plot
spider_plot(P,...
    'AxesLabels', {'JRRF', 'MFNN', 'L', 'DET', 'ENTR', 'RTE'},...
    'AxesLimits', [.5, -1.05, .7, .7, .7, .7; 1, 1, 1, 1, 1, 1],...
    'FillOption', 'on',...
    'LabelFontSize', fsla,...
    'FillTransparency', 0.2);
title('Rössler (y-component)','FontSize',fst,'FontName', 'Times New Roman')

legend('standard TDE', 'Garcia & Almeida', 'MDOP', 'PECUZAL', 'Location', 'southoutside','Fontsize',fsl);

% Rössler multi
% Initialize data points
TDE = [.78 .31 0.96 .9999 .93 .95];
GA = [.8 -0.06 0.93 .9997 .84 .93];
MDOP = [0.8 .68 .99 .9999 .94 .97];
PECUZAL = [.865 1 1 1 .94 .97];
P = [TDE; GA; MDOP; PECUZAL];

figure('Units','normalized','Position',[.3 .3 .8 .8])
% Spider plot
spider_plot(P,...
    'AxesLabels', {'JRRF', 'MFNN', 'L', 'DET', 'ENTR', 'RTE'},...
    'AxesLimits', [.7, -0.1, .7, .7, .7, .7; 1, 1, 1, 1, 1, 1],...
    'FillOption', 'on',...
    'LabelFontSize', fsla,...
    'FillTransparency', 0.2);
title('Rössler (x- & y-component)','FontSize',fst,'FontName', 'Times New Roman')

legend('standard TDE', 'Garcia & Almeida', 'MDOP', 'PECUZAL', 'Location', 'southoutside','Fontsize',fsl);


% Duffing uni
% Initialize data points
TDE = [.78 .9 .97 .9999 .93 .95];
GA = [.58 -1.1 .82 0.9999 .71 .81];
MDOP = [.75 .9 1 0.9999 .9 .9];
PECUZAL = [.78 1 1 0.9999 .92 .95];
P = [TDE; GA; MDOP; PECUZAL];

figure('Units','normalized','Position',[.3 .3 .8 .8])
% Spider plot
spider_plot(P,...
    'AxesLabels', {'JRRF', 'MFNN', 'L', 'DET', 'ENTR', 'RTE'},...
    'AxesLimits', [.5, -1.1, .7, .7, .4, .7; 1, 1, 1, 1, 1, 1],...
    'FillOption', 'on',...
    'LabelFontSize', fsla,...
    'FillTransparency', 0.2);
title('Duffing (x-component)','FontSize',fst,'FontName', 'Times New Roman')

legend('standard TDE', 'Garcia & Almeida', 'MDOP', 'PECUZAL', 'Location', 'southoutside','Fontsize',fsl);

% Duffing multi
% Initialize data points
TDE = [.78 1 .78 .9999 .93 .95];
GA = [.79 .17 .98 .9999 0.52 .98];
MDOP = [.84 .74 .95 0.999 0.45 .97];
PECUZAL = [.84 .43 1 0.999 0.56 .98];
P = [TDE; GA; MDOP; PECUZAL];

figure('Units','normalized','Position',[.3 .3 .8 .8])
% Spider plot
spider_plot(P,...
    'AxesLabels', {'JRRF', 'MFNN', 'L', 'DET', 'ENTR', 'RTE'},...
    'AxesLimits', [.5, -4.1, .7, .7, .4, .7; 1, 1, 1, 1, 1, 1],...
    'FillOption', 'on',...
    'LabelFontSize', fsla,...
    'FillTransparency', 0.2);
title('Duffing (x- & y-component)','FontSize',fst,'FontName', 'Times New Roman')

legend('standard TDE', 'Garcia & Almeida', 'MDOP', 'PECUZAL', 'Location', 'southoutside','Fontsize',fsl);