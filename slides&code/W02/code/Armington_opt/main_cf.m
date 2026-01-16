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


%% Counterfactual: China's unilateral optimal import tariffs

i = 9; % China

% Lower and upper bound
hX_lb = 0*ones(m.N, 1);
hP_lb = 0*ones(m.N, 1);
hw_lb = 0*ones(m.N, 1);
ttar_lb = -1*ones(m.N, 1);
ttar_lb(i) = 0;
lb = [hX_lb; hP_lb; hw_lb; ttar_lb];

hX_ub = inf*ones(m.N, 1);
hP_ub = inf*ones(m.N, 1);
hw_ub = inf*ones(m.N, 1);
ttar_ub = inf*ones(m.N, 1);
ttar_ub(i) = 0;
ub = [hX_ub; hP_ub; hw_ub; ttar_ub];

% Initial guesses
x0 = [ones(3*m.N, 1); m.tar(:,i)]; 

tic
[x,fval,exitflag,output,lambda,grad,hess] = ...
   knitromatlab(@(x)obj_tar(x, i, m),x0,[],[],[],[],lb,ub,@(x)con_eqm(x, i ,m),[],[],'nlpoptions.opt');
toc

% Consequences of the unilaterally optimal tariffs

hX = x(1: m.N);
hP = x(m.N + 1: 2*m.N);
hw = x(2*m.N + 1: 3*m.N);
ttar_cn = x(3*m.N + 1: 4*m.N);

ttar = m.tar;
ttar(:,i) = ttar_cn;
htau = ones(m.N, m.N);

[hw, hP, hX, hlambda, hwel] = eqm_iter(ttar, htau, m);
dwel_opt = (hwel-1)*100;

%% Nash tariffs between the U.S. and China

i = 9;  % China
j = 36; % U.S.

nash_tar = [m.tar(:,i), m.tar(:,j)]; % initial guesses for Nash tariffs
htau = ones(m.N, m.N);

% Nash tariffs
diff = 1;
tol = 1e-3;
cc = 0;
ww = 0.2;

while diff>tol && cc <100
    nash_tar_out = zeros(m.N,2);

    % China's optimal tariffs
    hX_lb = 0*ones(m.N, 1);
    hP_lb = 0*ones(m.N, 1);
    hw_lb = 0*ones(m.N, 1);
    ttar_lb = -1*ones(m.N, 1);
    ttar_lb(i) = 0;
    lb = [hX_lb; hP_lb; hw_lb; ttar_lb];
    
    hX_ub = inf*ones(m.N, 1);
    hP_ub = inf*ones(m.N, 1);
    hw_ub = inf*ones(m.N, 1);
    ttar_ub = inf*ones(m.N, 1);
    ttar_ub(i) = 0;
    ub = [hX_ub; hP_ub; hw_ub; ttar_ub];
    x0 = [ones(3*m.N, 1); m.tar(:,i)]; 
    [x,fval,exitflag,output,lambda,grad,hess] = ...
       knitromatlab(@(x)obj_tar(x, i, m),x0,[],[],[],[],lb,ub,@(x)con_eqm(x, i ,m),[],[],'nlpoptions.opt');
    nash_tar_out(:,1)= x(3*m.N + 1: 4*m.N);

    % The US optimal tariffs
    hX_lb = 0*ones(m.N, 1);
    hP_lb = 0*ones(m.N, 1);
    hw_lb = 0*ones(m.N, 1);
    ttar_lb = -1*ones(m.N, 1);
    ttar_lb(j) = 0;
    lb = [hX_lb; hP_lb; hw_lb; ttar_lb];
    
    hX_ub = inf*ones(m.N, 1);
    hP_ub = inf*ones(m.N, 1);
    hw_ub = inf*ones(m.N, 1);
    ttar_ub = inf*ones(m.N, 1);
    ttar_ub(j) = 0;
    ub = [hX_ub; hP_ub; hw_ub; ttar_ub];
    x0 = [ones(3*m.N, 1); m.tar(:,j)]; 
    [x,fval,exitflag,output,lambda,grad,hess] = ...
       knitromatlab(@(x)obj_tar(x, j, m),x0,[],[],[],[],lb,ub,@(x)con_eqm(x, j ,m),[],[],'nlpoptions.opt');
    nash_tar_out(:,2)= x(3*m.N + 1: 4*m.N);

    diff = max(abs(nash_tar_out(:)-nash_tar(:)));
    cc = cc+1;
    disp(diff)
    disp(cc)

    nash_tar = ww*nash_tar_out+(1-ww)*nash_tar;
    m.tar(:,i) = nash_tar(:,1);
    m.tar(:,j) = nash_tar(:,2);
end

disp(nash_tar)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    %disp(diff)
    
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

%% Objective function for the optimal tariff
function [f,g] = obj_tar(x, i, m)

g = [];

% variables
hX = x(1: m.N);
hP = x(m.N + 1: 2*m.N);
hw = x(2*m.N + 1: 3*m.N);
ttar = x(3*m.N + 1: 4*m.N);

% tariff matrix
ttar_mat = m.tar;
ttar_mat(:, i) = ttar;

htar = (1 + ttar_mat)./(1 + m.tar);
psi = htar.^(-m.theta).*repmat((hw).^(-m.theta), [1,m.N]);
hlambda = psi./repmat(sum(m.lambda.*psi,1), [m.N, 1]);

hwel = (hw.*m.wL + sum(ttar_mat./(1 + ttar_mat).*hlambda.*m.lambda.*repmat(hX'.*m.X',[m.N, 1]) ,1 )' )./(m.wL + sum(m.tar./(1 + m.tar).*m.lambda.*repmat(m.X',[m.N, 1]) ,1 )' )./hP; % welfare

f = -hwel(i);
end


%% Equilibrium constraints for the optimal tariff
function [c,ceq] = con_eqm(x, i ,m)

% variables
hX = x(1: m.N);
hP = x(m.N + 1: 2*m.N);
hw = x(2*m.N + 1: 3*m.N);
ttar = x(3*m.N + 1: 4*m.N);

% tariff matrix
ttar_mat = m.tar;
ttar_mat(:, i) = ttar;

htar = (1 + ttar_mat)./(1 + m.tar);
psi = htar.^(-m.theta).*repmat((hw).^(-m.theta), [1,m.N]);
hlambda = psi./repmat(sum(m.lambda.*psi,1), [m.N, 1]);
hw_out = sum(hlambda.*m.lambda./(1 +ttar_mat).*repmat(hX'.*m.X',[m.N,1]),2)./m.wL;
hX_out = (hw.*m.wL + sum(ttar_mat./(1+ttar_mat).*hlambda.*m.lambda.*repmat(hX'.*m.X',[m.N, 1]),1)')./m.X;
hP_out = (sum(m.lambda.*psi,1)').^(-1/m.theta);

ceq = [hw - hw_out; hX - hX_out; hP - hP_out; hw(1) - 1; ttar(i)]; 
c = [];
end






