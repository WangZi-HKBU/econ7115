function [y,gy] = func_mm(theta,m)
% Moment function in GMM and its Jacobian matrix

u_vec = m.y - 1./(1+exp(-m.x*theta));

[~,K]=size(m.z);
m_mat = m.z.*repmat(u_vec,[1,K]); % dim1: m.N; dim2: K=#instruments

% Moment
y = sum(m_mat,1)'/m.N; % moment vector: K*1

% Jacobian
[~,J] = size(m.x);

g_3d = ones(K,J,m.N);
for k = 1:K
    for j = 1:J
        for i = 1:m.N
            g_3d(k,j,i) = -exp(-m.x(i,:)*theta)/(1+exp(-m.x(i,:)*theta))^2*m.z(i,k)*m.x(i,j);
        end
    end
end

gy = sum(g_3d,3)/m.N; %dim1:K=#instruments; dim2:J=#parameters

end