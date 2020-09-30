% Farmer's relation of the delay in the Mackey-Glass-equation and the
% dimensionality of the corresponding attractor.
% Farmer, Physica 4D (1982, 366-393)

clear, clc

% fill delays and dimension from Table 1 
Delay = [17 23 23.8 30];
Dimension = [2.1 2.82 3.04 3.58];
dy = [0.02 0.03 0.03 0.04];

fit = fitlm(Delay,Dimension,'Weights',1./(dy.^2));

%% Linear model
figure
plot(fit)
fit.Coefficients(1,1)
fit.Coefficients(2,1)

b = 0.092759;
m = 0.11923;

%% Adapt our data
dim_for_delay_44 = b + m*44;