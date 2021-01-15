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

width2 = .9;

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
h = bar(results,width2);
xticklabels(group_labels)
ylabel("rel. deviation from reference")
set(gca,'TickLabelInterpreter', 'tex');
lgd = legend(legend_labels,'box','off','Location','best');
lgd.FontSize = ls;
grid on
set(gca,'LineWidth',2)
set(gca,'FontSize',fs)

%%
