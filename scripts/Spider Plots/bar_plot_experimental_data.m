% We plot the results from script `analyze_standard_ecoscillator_try2.jl`

clear, clc

% get the same color sequence for the different methods as in the Spider
% plots

h = plot(1:10,1:10,1:10,2:11,1:10,3:12,1:10,4:13);
c = get(h,'Color');

c1 = c{1};
c2 = c{2};
c3 = c{3};
c4 = c{4};

close all


data1 = load("../Experimental data/Standard electrochemical chaos/s3ch4.txt");
data2 = load("../Experimental data/Standard electrochemical chaos/s1ch3.txt");
% length of the whole time series
M1 = length(data1);
M2 = length(data2);
M = min([M1, M2]);

% downsampling
data1 = data1(1:2:M);
M = length(data1);

% timestamp
t = 0:0.02:0.02*(M-1);

% regularize time series and add minimal noise component
data1 = data1 + .00000001*randn(length(data1),1);
data1 = (data1 - mean(data1))/std(data1);

% draw some subsamples
idxs = randsample(M-2500, 1);

% construct RP
idx1 = idxs(1);

x1 = data1(idx1:idx1+2500);
t_sub = t(idx1:idx1+2500)/10;


%% 3D Bar Plot

legend_labels = ["Cao's TDE", "Garcia & Almeida", "MDOP", "PECUZAL"];
group_labels = {'ENTR', 'LAM', 'RTE', "{T}"};

results_exp_1 = load("relative_dev_RQA_1.csv");
results = reshape(results_exp_1, [4,4]);

[best, best_idx] = min(results, [], 2)

fs = 26; % fontsize
ls = 34; % legendfontsize

width = .8; % bar width
 
% figure('Units','normalized','Position',[.2 .2 .8 .8])
% h = bar3(results,width);
% set(h(1), 'FaceColor', c1);
% set(h(2), 'FaceColor', c2);
% set(h(3), 'FaceColor', c3);
% set(h(4), 'FaceColor', c4);
% xticklabels([])
% yticklabels(group_labels)
% set(gca,'TickLabelInterpreter', 'tex');
% zlabel("rel. deviation from reference")
% lgd = legend(legend_labels,'box','off','Location','best');
% lgd.FontSize = ls;
% set(gca,'LineWidth',2)
% set(gca,'FontSize',fs)
% view(60,30)

figure('Units','normalized','Position',[.01 .01 .9 .9])

subplot(2,2,1)
plot(data1, 'LineWidth', 1); hold on
xlabel('time')
ylabel('current [mA]')
xticklabels([])
grid on
set(gca,'LineWidth',2)
set(gca,'Fontsize', fs)

subplot(2,2,3)
plot(t_sub,x1, '.-')
xlabel('time [ms]')
ylabel('current [mA]')
xlim([t_sub(1) t_sub(end)])
grid on
set(gca,'LineWidth',2)
set(gca,'Fontsize', fs)

subplot(2,2,[2,4])
h = bar3(results,width);
set(h(1), 'FaceColor', c1);
set(h(2), 'FaceColor', c2);
set(h(3), 'FaceColor', c3);
set(h(4), 'FaceColor', c4);
xticklabels([])
yticklabels(group_labels)
set(gca,'TickLabelInterpreter', 'tex');
zlabel("rel. deviation from reference")
lgd = legend(legend_labels,'box','off','Location','best');
lgd.FontSize = ls;
set(gca,'LineWidth',2)
set(gca,'FontSize',fs)
view(60,30)

%%
% 
% %% Plot time series and sub-RPs
% 
% clear, clc
% 
% data1 = load("../Experimental data/Standard electrochemical chaos/s3ch4.txt");
% data2 = load("../Experimental data/Standard electrochemical chaos/s1ch3.txt");
% % length of the whole time series
% M1 = length(data1);
% M2 = length(data2);
% M = min([M1, M2]);
% 
% % downsampling
% data1 = data1(1:2:M);
% data2 = data2(1:2:M);
% M = length(data1);
% 
% % regularize time series and add minimal noise component
% data1 = data1 + .00000001*randn(length(data1),1);
% data2 = data2 + .00000001*randn(length(data2),1);
% 
% data1 = (data1 - mean(data1))/std(data1);
% data2 = (data2 - mean(data2))/std(data2);
% 
% % draw some subsamples
% idxs = randsample(M-5000, 3);
% 
% % construct RP
% idx1 = idxs(1);
% idx2 = idxs(2);
% idx3 = idxs(3);
% 
% x1 = data1(idx1:idx1+5000);
% x2 = data1(idx2:idx2+2500);
% x3 = data1(idx3:idx3+2500);
% 
% % compute pecuzal embedding (toolbox need to be installed from Matlab Central)
% [Y1, tau_vals1, ts_vals1, epsilon_mins1, Ls1] = ...
%                 pecuzal_embedding(x1, 0:150, 'theiler', 20);
% [Y2, tau_vals2, ts_vals2, epsilon_mins2, Ls2] = ...
%                 pecuzal_embedding(x2, 0:150, 'theiler', 20);
% [Y3, tau_vals3, ts_vals3, epsilon_mins3, Ls3] = ...
%                 pecuzal_embedding(x3, 0:150, 'theiler', 20);
% 
% %% create sub-RPs
% 
% RP1 = rp(Y1,.08,'var');
% RP2 = rp(Y2,.08,'var');
% RP3 = rp(Y3,.08,'var');
% 
% %% Plot stuff
% 
% t = 0:0.02:0.02*M;
% 
% t_sub = t(idx1:idx1+5000)/10;
% figure('Units','normalized','Position',[.2 .2 .4 .3])
% plot(t_sub,x1, '.-')
% xlabel('time [ms]')
% ylabel('Voltage [\muV]')
% xlim([t_sub(1) t_sub(end)])
% grid on
% set(gca,'LineWidth',2)
% set(gca,'Fontsize', 14)
% 
% figure('Units','normalized','Position',[.2 .2 .4 .5])
% plot3(Y1(:,1),Y1(:,2),Y1(:,4))
% title('Reference Reconstruction')
% set(gca,'XTickLabel',[], 'YTickLabel', [], 'ZTickLabel', [])
% grid on
% set(gca,'LineWidth',2)
% set(gca,'Fontsize', 14)
% 
% figure('Units','normalized','Position',[.2 .2 .4 .5])
% plot3(Y2(:,1),Y2(:,2),Y2(:,3))
% title('Subsample Reconstruction')
% set(gca,'XTickLabel',[], 'YTickLabel', [], 'ZTickLabel', [])
% grid on
% set(gca,'LineWidth',2)
% set(gca,'Fontsize', 14)
% 
% figure('Units','normalized','Position',[.2 .2 .4 .5])
% plot3(Y3(:,1),Y3(:,2),Y3(:,3))
% title('Subsample Reconstruction')
% set(gca,'XTickLabel',[], 'YTickLabel', [], 'ZTickLabel', [])
% grid on
% set(gca,'LineWidth',2)
% set(gca,'Fontsize', 14)
% %%
% 
% figure('Units','normalized','Position',[.2 .2 .4 .5])
% imagesc(RP1), colormap([1 1 1; 0 0 0]), axis xy square
% title('Reference RP')
% set(gca,'XTickLabel',[], 'YTickLabel', [])
% set(gca,'LineWidth',2)
% set(gca,'Fontsize', 14)
% 
% 
% figure('Units','normalized','Position',[.2 .2 .4 .5])
% imagesc(RP2), colormap([1 1 1; 0 0 0]), axis xy square
% title('Subsample RP')
% set(gca,'XTickLabel',[], 'YTickLabel', [])
% set(gca,'LineWidth',2)
% set(gca,'Fontsize', 14)
% 
% figure('Units','normalized','Position',[.2 .2 .4 .5])
% imagesc(RP3), colormap([1 1 1; 0 0 0]), axis xy square
% title('Subsample RP')
% set(gca,'XTickLabel',[], 'YTickLabel', [])
% set(gca,'LineWidth',2)
% set(gca,'Fontsize', 14)
% 
% %%
% RQA_distr = load("RQA_distribution.csv");
% 
% fs = 10;
% 
% figure('Units','normalized','Position',[.2 .2 .4 .5])
% subplot(2,2,1)
% histogram(RQA_distr(1,:))
% title('ENTR')
% grid on
% set(gca,'XTickLabel',[], 'YTickLabel', [])
% set(gca,'LineWidth',1.5)
% set(gca,'Fontsize', fs)
% 
% subplot(2,2,2)
% histogram(RQA_distr(2,:))
% title('LAM')
% grid on
% set(gca,'XTickLabel',[], 'YTickLabel', [])
% set(gca,'LineWidth',1.5)
% set(gca,'Fontsize', fs)
% 
% subplot(2,2,3)
% histogram(RQA_distr(3,:))
% title('RTE')
% grid on
% set(gca,'XTickLabel',[], 'YTickLabel', [])
% set(gca,'LineWidth',1.5)
% set(gca,'Fontsize', fs)
% 
% subplot(2,2,4)
% histogram(RQA_distr(4,:))
% title('TRANS')
% grid on
% set(gca,'XTickLabel',[], 'YTickLabel', [])
% set(gca,'LineWidth',1.5)
% set(gca,'Fontsize', fs)
