clear
clc
close all

%% Predetermined arameters
m.alpha = 1/3;
m.delta = 0.06;
m.beta = 0.9;

m.k =(0.01:0.01:10)'; % state space for capital k
m.J = length(m.k);

%% Simulated shocks
%{
m.S = 10000;
m.T = 300;

m.k0_vec = ones(m.S,1); % initial capital stock
m.z0_index = ones(m.S,1); % initial productivity
for s=1:m.S
    m.k0_vec(s) = m.k(randperm(m.J,1));
    m.z0_index(s) = randperm(2,1);
end

m.u_mat = rand(m.S,m.T-1);

save('temp/sim_shock.mat', 'm');
%}

load 'temp/sim_shock.mat'


%% Test
c = 0.1;
piHH = 0.2;
piLL = 0.6;
zgap = 0.3;
theta = [c;piHH;piLL;zgap];

[k_mat_sim,d_mat_sim] = func_sim(theta,m);











