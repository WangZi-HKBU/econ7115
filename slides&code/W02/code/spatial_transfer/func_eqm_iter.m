function [w,L,P,welfare] = func_eqm_iter(s,m)
% Solve the equilibrium by iteration
% s = [s_1,s_2,...s_{N-1}]'; s_N is determined by the government budget balance

w = 1/m.N*ones(m.N,1);
L = m.barL/m.N*ones(m.N,1);
P = ones(m.N,1);

diff = 1;
tol = 1e-6;
cc = 0;
ww = 0.1;

while diff>tol && cc < 10000
    [w_new,L_new,P_new,welfare] = func_eqm_update(w,L,P,s,m);
    r = [w;L;P];
    r_new = [w_new;L_new;P_new];

    diff = max(abs(r-r_new));
    cc = cc+1;
    %disp(diff)
    %disp(cc)

    w = ww*w_new + (1-ww)*w;
    L = ww*L_new + (1-ww)*L;
    P = ww*P_new + (1-ww)*P;
end

end