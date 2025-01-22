function [w_new,L_new,P_new,welfare] = func_eqm_update(w,L,P,s,m)
% Update the equilibrium outcomes by equilibrium conditions
% s = [s_1,s_2,...s_{N-1}]'; s_N is determined by the government budget balance

s_last = -sum(s.*w(1:end-1).*L(1:end-1))/(w(end)*L(end));
s_vec = [s;s_last];

K_mat = (repmat(w./(m.barA.*L.^(m.alpha)),[1,m.N]).*m.tau).^(1-m.sigma); % dim1:i; dim2:n
lambda_mat = K_mat./repmat(sum(K_mat,1),[m.N,1]);

P_new = sum(K_mat,1)'.^(1/(1-m.sigma));

w_new = sum(lambda_mat.*repmat((1+s_vec)'.*w'.*L',[m.N,1]),2)./L; 

B_vec = (m.barB.*(1+s_vec).*w./P).^(1/m.beta);
L_new = m.barL*B_vec/sum(B_vec);

welfare = sum(B_vec).^(m.beta);

wnum = sum(w_new);
w_new = w_new/wnum;
P_new = P_new/wnum;


end