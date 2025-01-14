function [w,X,P,welfare,lambda_mat] = func_eqm_iter(tar,m)
% Solving the equilibrium outcomes by simple iteration

w = ones(m.N,1)/m.N;
X = w;
diff = 1;
tol = 1e-6;
ww = 0.2;
cc = 0;

while diff>tol && cc<100
    [w_new,X_new,P,welfare,lambda_mat] = func_eqm_update(w,X,tar,m);
    r = [w;X];
    r_new = [w_new;X_new];
    diff = max(abs(r-r_new));
    cc = cc+1;
    %disp(diff)
    %disp(cc)
    
    w = ww*w_new+(1-ww)*w;
    X = ww*X_new+(1-ww)*X;
end

end

