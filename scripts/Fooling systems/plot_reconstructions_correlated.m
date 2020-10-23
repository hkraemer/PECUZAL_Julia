
%% Plot PS-reconstructions 

clear, clc

Y_GA = load("./correlated results/Y_GA.csv");
Y_mdop = load("./correlated results/Y_mdop.csv");
Y_pec = load("./correlated results/Y_pec.csv");
tau_vals_GA = load("./correlated results/taus_GA.csv");
tau_vals_mdop = load("./correlated results/taus_mdop.csv");
tau_vals_pec = load("./correlated results/taus_pec.csv");
ts_vals_GA = load("./correlated results/ts_GA.csv");
ts_vals_mdop = load("./correlated results/ts_mdop.csv");
ts_vals_pec = load("./correlated results/ts_pec.csv");


fs = 20;

figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(131)
plot3(Y_GA(:,1),Y_GA(:,2),Y_GA(:,3))
str1 = sprintf('Garcia & Almeida',1.4);
str2 = sprintf(strcat("delays: [",num2str(tau_vals_GA'),"]"));
str3 = sprintf(strcat("ts: [",num2str(ts_vals_GA'),"]"));
title({str1,str2,str3});
xticklabels([])
yticklabels([])
zticklabels([])
set(gca,'LineWidth',2)
set(gca,'FontSize',fs)
grid on


subplot(132)
plot3(Y_mdop(:,2),Y_mdop(:,3),Y_mdop(:,1))
str1 = sprintf('MDOP',1.4);
str2 = sprintf(strcat("delays: [",num2str(tau_vals_mdop'),"]"));
str3 = sprintf(strcat("ts: [",num2str(ts_vals_mdop'),"]"));
title({str1,str2,str3});
xticklabels([])
yticklabels([])
zticklabels([])
set(gca,'LineWidth',2)
set(gca,'FontSize',fs)
grid on


subplot(133)
plot3(Y_pec(:,3),Y_pec(:,1),Y_pec(:,2))
str1 = sprintf('PECUZAL',1.4);
str2 = sprintf(strcat("delays: [",num2str(tau_vals_pec'),"]"));
str3 = sprintf(strcat("ts: [",num2str(ts_vals_pec'),"]"));
title({str1,str2,str3});
xticklabels([])
yticklabels([])
zticklabels([])
set(gca,'LineWidth',2)
set(gca,'FontSize',fs)
grid on