clear
clc
close all

%% Parameters

m.N = 3;
m.sigma = 4;
m.A = [2;1;1];
m.L = [1;2;4];
m.tau = 2*(1-eye(m.N))+eye(m.N);
tar = zeros(m.N,m.N);

%% Baseline equilibrium

[w,X,P,welfare,lambda_mat] = func_eqm_iter(tar,m);

welfare_baseline = welfare;
RI_baseline = w./P;

%% Unilateral tariff by Country 1

tar = zeros(m.N,m.N);
tar(2,1) = 0.2;
tar(3,1) = 0.2;

[w,X,P,welfare,lambda_mat] = func_eqm_iter(tar,m);

welfare_uni = welfare;
RI_uni = w./P;

disp([welfare_baseline, welfare_uni])
disp([RI_baseline, RI_uni])

%}
