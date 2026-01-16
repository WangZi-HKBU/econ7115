function [w_new,X_new,P,welfare,lambda_mat] = func_eqm_update(w,X,tar,m)
% Update the equilibrium outcomes by equilibrium conditions

kappa = m.tau.*(1+tar);
K_mat = (repmat(w./m.A,[1,m.N]).*kappa).^(1-m.sigma); % dim1:i; dim2:n
P = sum(K_mat,1)'.^(1/(1-m.sigma));
welfare = X./m.L./P;

lambda_mat = K_mat./repmat(sum(K_mat,1),[m.N,1]); % dim1:i; dim2:n
Xin_mat = lambda_mat.*repmat(X',[m.N,1]); % dim1:i; dim2:n

w_new = sum(Xin_mat./(1+tar),2)./m.L;
X_new = w.*m.L+sum(tar./(1+tar).*Xin_mat,1)';

wnum = sum(w_new);
w_new = w_new/wnum;
X_new = X_new/wnum;
P = P./wnum;


end

