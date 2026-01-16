clear
close all
clc

%% Parameters
m.sigma = 0.5; % u(c) = c^{1-\sigma}/(1-\sigma)
m.beta = 0.95; % discouting rate

W0 = 1;
m.J = 200; %#grids
w_vec = (1:m.J)'/m.J*W0;

z_vec = [0.8;1.2];
pi_mat = [
    0.8, 0.2
    0.3, 0.7];

%% Value function iteration: deterministic

v_vec = func_utility(w_vec,m); % initial guess

diff = 1;
tol = 1e-6;
cc = 0;
ww = 0.5;

while diff>tol&&cc<10000
    v_vec_out = ones(m.J,1);
    g_vec = ones(m.J,1); % policy function: W_{t+1} = g(W_t)
    for j=1:m.J
        rhs_vec = ones(j,1);
        for k=1:j
            rhs_vec(k) = func_utility(w_vec(j)-w_vec(k),m)+m.beta*v_vec(k);
        end
        v_vec_out(j) = max(rhs_vec);
        g_vec(j) = w_vec(rhs_vec==max(rhs_vec));
    end
    
    diff = max(abs(v_vec-v_vec_out));
    cc = cc+1;
    disp(diff)
    disp(cc)
    
    v_vec = ww*v_vec_out + (1-ww)*v_vec;
end

u_vec = func_utility(w_vec,m); % instantaneous utility

figure(1)
qx = w_vec;
h(1) = line(qx, v_vec, 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
h(2) = line(qx, g_vec, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);

legend(h, {'Value Function','Policy Function'},'Location','northwest','Orientation','vertical', 'FontSize',20,'Interpreter','latex')
xlabel('$W$', 'FontSize',22, 'Interpreter','latex')
set(gca,'fontsize',22)
set(gcf, 'Position', get(0, 'Screensize'))
%set(gca,'Xdir','reverse')
print('result/VFI', '-dpng')

%% Value function iteration: Stochastic

S = length(z_vec);
v_mat = repmat(func_utility(w_vec,m),[1,S]).*repmat(z_vec',[m.J,1]); % initial guess

diff = 1;
tol = 1e-6;
cc = 0;
ww = 0.5;

while diff>tol&&cc<10000
    v_mat_out = ones(m.J,S);
    g_mat = ones(m.J,S); % policy function: W_{t+1} = g(W_t,z_t)
    for j=1:m.J
        rhs_1 = ones(j,1);
        rhs_2 = ones(j,1);
        for k=1:j
            rhs_1(k) = z_vec(1)*func_utility(w_vec(j)-w_vec(k),m)+m.beta*(pi_mat(1,1)*v_mat(k,1)+pi_mat(1,2)*v_mat(k,2));
            rhs_2(k) = z_vec(2)*func_utility(w_vec(j)-w_vec(k),m)+m.beta*(pi_mat(2,1)*v_mat(k,1)+pi_mat(2,2)*v_mat(k,2));
        end
        v_mat_out(j,1) = max(rhs_1);
        v_mat_out(j,2) = max(rhs_2);
        g_mat(j,1) = w_vec(rhs_1==max(rhs_1));
        g_mat(j,2) = w_vec(rhs_2==max(rhs_2));
    end
    
    diff = max(abs(v_mat(:)-v_mat_out(:)));
    cc = cc+1;
    disp(diff)
    disp(cc)
    
    v_mat = ww*v_mat_out + (1-ww)*v_mat;
end

figure(2)
qx = w_vec;
h(1) = line(qx, v_mat(:,1), 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
h(2) = line(qx, v_mat(:,2), 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);

legend(h, {'$V(W,z_1)$','$V(W,z_2)$'},'Location','northwest','Orientation','vertical', 'FontSize',20,'Interpreter','latex')
xlabel('$W$', 'FontSize',22, 'Interpreter','latex')
set(gca,'fontsize',22)
set(gcf, 'Position', get(0, 'Screensize'))
%set(gca,'Xdir','reverse')
print('result/SVFI_VF', '-dpng')

figure(3)
qx = w_vec;
h(1) = line(qx, g_mat(:,1), 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
h(2) = line(qx, g_mat(:,2), 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);

legend(h, {'$g(W,z_1)$','$g(W,z_2)$'},'Location','northwest','Orientation','vertical', 'FontSize',20,'Interpreter','latex')
xlabel('$W$', 'FontSize',22, 'Interpreter','latex')
set(gca,'fontsize',22)
set(gcf, 'Position', get(0, 'Screensize'))
%set(gca,'Xdir','reverse')
print('result/SVFI_PF', '-dpng')


