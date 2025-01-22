function y = func_obj_2sgmm(theta,w_mat,m)
% The objective function of the two-step gmm estimator

[hm,~] = func_mm(theta,m);

y = hm'*w_mat*hm;

end