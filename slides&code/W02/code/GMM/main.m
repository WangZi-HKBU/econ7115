clear
clc
close all

%% Data generation
m.N = 1000; % sample size
m.x = 5*rand(m.N,6); % 6 explanatory variables

theta_true = [1;4.5;3.3;0;-10;-0.1]; % True coefficients
u = 0.1*randn(m.N,1);

m.y = 1./(1+exp(-m.x*theta_true))+u;
m.z = [m.x, m.x.^2];

csvwrite('data.csv',[m.y,m.x,m.z])

%data: m.x; m.y; m.z
%parameters to estimate: \theta

%% Two-step GMM

[~,K]=size(m.z);
[~,J] = size(m.x);
w_mat = eye(K);

theta0 = zeros(J,1); % initial guess

options = optimoptions('fminunc','Algorithm','quasi-newton','Display','iter',...
    'MaxIterations',500,'MaxFunctionEvaluations',100000,'OptimalityTolerance',1e-8);

[theta_p,fval,exitflag,output] = fminunc(@(theta)func_obj_2sgmm(theta,w_mat,m),theta0,options);

w_mat = func_optw(theta_p,m);

[htheta,fval,exitflag,output] = fminunc(@(theta)func_obj_2sgmm(theta,w_mat,m),theta_p,options);

disp([theta_true,theta_p,htheta])

