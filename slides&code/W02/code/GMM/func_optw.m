function y = func_optw(theta,m)
% Optimal weighting matrix given the estimate of \theta

u_vec = m.y - 1./(1+exp(-m.x*theta));

[~,K]=size(m.z);
m_mat = m.z.*repmat(u_vec,[1,K]); % dim1: m.N; dim2: K=#instruments

Vm = m_mat'*m_mat/m.N-(sum(m_mat,1)'/m.N)*(sum(m_mat,1)/m.N);

y = inv(Vm);

end