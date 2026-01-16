clear
clc
close all

%% Parameters

m.N = 2;
m.barL = 1;
m.sigma = 4;
m.beta = 0.5;
m.alpha = 0.1; % alpha < beta: unique equilibirum

m.barB = ones(m.N,1);
m.barA = [1.2; 1];
m.tau = 1.3*(1-eye(m.N))+eye(m.N);

%% Test s_1

s = -0.2;
[w,L,P,welfare] = func_eqm_iter(s,m);

%% Optimal s_1

ss_vec = (-0.2:0.002:0.09)';
S = length(ss_vec);

wel_vec = ones(S,1);
L_mat = ones(S,m.N);

for ss = 1:S
    s = ss_vec(ss);
    [w,L,P,welfare] = func_eqm_iter(s,m);

    wel_vec(ss) = welfare;
    L_mat(ss,:) = L';

    disp(s)
end

figure(1)
qx = ss_vec;
h(1) = line(qx, L_mat(:,1), 'Color', 'r', 'LineStyle', '-.', 'LineWidth', 2);
h(2) = line(qx, L_mat(:,2), 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);

legend(h, {'$L_1$','$L_2$'},'Location','northwest','Orientation','vertical', 'FontSize',20,'Interpreter','latex')
xlabel('$s_1$', 'FontSize',22, 'Interpreter','latex')
ylabel('Labor', 'FontSize',22)
set(gca,'fontsize',22)
set(gcf, 'Position', get(0, 'Screensize'))
%set(gca,'Xdir','reverse')
print('result/L', '-dpng')

figure(2)
qx = ss_vec;
h = line(qx, wel_vec, 'Color', 'r', 'LineStyle', '-.', 'LineWidth', 2);

xlabel('$s_1$', 'FontSize',22, 'Interpreter','latex')
ylabel('Welfare', 'FontSize',22)
set(gca,'fontsize',22)
set(gcf, 'Position', get(0, 'Screensize'))
%set(gca,'Xdir','reverse')
print('result/welfare', '-dpng')

disp('optimal s_1=')
disp(ss_vec(wel_vec==max(wel_vec)))


