% Farmer's relation of the delay in the Mackey-Glass-equation and the
% dimensionality of the corresponding attractor.
% Farmer, Physica 4D (1982, 366-393)

clear, clc

% fill delays and dimension from Table 1 
Delay = [17 23 23.8 30];
Dimension = [2.1 2.82 3.04 3.58];
dy = [0.02 0.03 0.03 0.04];

% fit = fitlm(Delay,Dimension,'Weights',1./(dy.^2));
fit = fitlm(Delay,Dimension);

ci = coefCI(fit)

%% Linear model
figure
plot(fit)
% fit.Coefficients(1,1)
% fit.Coefficients(2,1)

% b = 0.092759;
% m = 0.11923;

b = 0.20385;
b_u = -0.8719;
b_o = 1.2796;
m = 0.11433;
m_u = 0.0693;
m_o = 0.1593;

%% Adapt our data
dim_for_delay_44 = b + m*44
under_ci = b_u + m_u*44
upper_ci = b_o + m_o*44

for i = 1:50
    y(i) = b + m*i;
    y_u(i) = b_u + m_u*i;
    y_o(i) = b_o + m_o*i;
end

figure
plot(1:50, y, 'LineWidth',3), hold on
plot(1:50, y_u, 'r-', 'LineWidth',2), hold on
plot(1:50, y_o, 'r-', 'LineWidth',2), hold on
scatter(Delay,Dimension,80,'filled')
xline(44)
title('Delay-attractor-dimension-relation after Farmer (Kaplan-Yorke-Conjecture)')
legend('least square linear fit', 'lower 95% CI','upper 95% CI','data')
xlabel('delay')
ylabel('dimension')
set(gca,'LineWidth',2)
set(gca,'FontSize',14)
grid on

