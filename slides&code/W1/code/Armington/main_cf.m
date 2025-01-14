clear
clc
close all

%% Data and Parameters
m.N = 36;
m.sigma = 3.059;
m.theta = m.sigma-1;

% Load the data
xtr2001 = csvread('data_matlab/xtr2001.csv', 1,1); % dim 1: iso_o; dim 2: iso_d
tar2001 = csvread('data_matlab/tariff2001.csv', 1,1); % dim 1: iso_o; dim 2: iso_d

xtr = xtr2001;
tar = tar2001;

X = sum(xtr, 1)'; % expenditure X_n
lambda = xtr./repmat(sum(xtr,1),[m.N, 1]);
wL = sum(1./(1+tar).*xtr, 2);

m.X = X;
m.lambda = lambda;
m.wL = wL;
m.tar = tar;



%%  Trade balance Adjustment
ttar = m.tar;
htau = ones(m.N, m.N);

tic
[hw, hP, hX, hlambda, hwel] = eqm_iter(ttar, htau, m);
toc

m.lambda = m.lambda.*hlambda; % Balanced trade share
m.X = m.X.*hX; % Balanced trade expenditure
m.wL = m.wL.*hw; % Balanced trade wage

% Test the trade balance
[hw, hP, hX, hlambda, hwel] = eqm_iter(ttar, htau, m);


%% Counterfactual: China eliminates all of its import tariffs

ttar = m.tar;
i = 9; % China
ttar(:,i) = zeros(m.N, 1);
htau = ones(m.N, m.N);

[hw, hP, hX, hlambda, hwel] = eqm_iter(ttar, htau, m);

dwel = (hwel-1)*100;
drw = (hw./hP - 1)*100;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solving the Equilibrium by Iteration
function [hw, hP, hX, hlambda, hwel] = eqm_iter(ttar, htau, m)
% Initial values
hw = ones(m.N, 1);
hX = ones(m.N, 1);

% iteration parameters
tol = 1e-10;
ww = 0.1;
cc = 0;
diff = 1;

while diff > tol && cc < 1000
    htar = (1 + ttar)./(1 + m.tar);
    psi = (htar.*htau).^(-m.theta).*repmat((hw).^(-m.theta), [1,m.N]);
    hlambda = psi./repmat(sum(m.lambda.*psi,1), [m.N, 1]);
    hw_out = sum(hlambda.*m.lambda./(1 + ttar).*repmat(hX'.*m.X',[m.N,1]),2)./m.wL;
    hX_out = (hw.*m.wL + sum(ttar./(1+ttar).*hlambda.*m.lambda.*repmat(hX'.*m.X',[m.N, 1]),1)')./m.X;
    
    r = [hw; hX];
    rr = [hw_out; hX_out];
    diff = max(abs(r-rr));
    disp(diff)
    
    wnum = hw(1);
    hw = hw./wnum;
    hX = hX./wnum;
    
    hw = ww*hw_out + (1 - ww)*hw;
    hX = ww*hX_out + (1 - ww)*hX;
    cc = cc+1;
end

hP = (sum(m.lambda.*psi,1)').^(-1/m.theta);
hP = hP./wnum;


hwel = (hw.*m.wL + sum(ttar./(1+ttar).*hlambda.*m.lambda.*repmat(hX'.*m.X',[m.N, 1]),1)')./(m.wL + sum(m.tar./(1+m.tar).*m.lambda.*repmat(m.X',[m.N, 1]),1)')./hP;

end








