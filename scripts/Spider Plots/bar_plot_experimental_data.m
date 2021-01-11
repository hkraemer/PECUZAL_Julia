% We plot the results from script `analyze_standard_ecoscillator_try2.jl`

clear, clc

legend_labels = ["Cao's TDE", "Garcia & Almeida", "MDOP", "PECUZAL"];
group_labels = categorical({'ENTR', 'LAM', 'RTE', 'TRANS'});

results_exp_1 = load("relative_dev_RQA_1.csv");
results = reshape(results_exp_1, [4,4]);

fs = 26; % fontsize
ls = 30; % legendfontsize

figure('Units','normalized','Position',[.2 .2 .8 .8])
bar(group_labels,results, 1, 'EdgeColor', 'none')
lgd = legend(legend_labels, 'box', 'off', 'Location','northwest');
ylabel("relative deviation to reference")
lgd.FontSize = ls;
ylim([-0.05 0.3])
set(gca,'LineWidth',2)
set(gca,'FontSize',fs)
grid on
