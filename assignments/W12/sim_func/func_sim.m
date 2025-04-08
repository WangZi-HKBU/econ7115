function [k_mat_sim,d_mat_sim] = func_sim(theta,m)
% Simulate capital and cash flows

m.c = theta(1);
m.piHH = theta(2);
m.piLL = theta(3);
m.zgap = theta(4);

m.zL = 1-m.zgap;
m.zH = 1+m.zgap;

%% Value function iteration

% initial guess for v_H(k) and v_L(k)
vH = m.zH*m.k.^m.alpha;
vL = m.zL*m.k.^m.alpha;

diff = 1;
tol = 1e-6;
cc = 0;
ww = 1;

while diff>tol && cc<10000
    [vH_out,vL_out,jH,jL] = func_value_update(vH,vL,m);
    r = [vH;vL];
    r_out = [vH_out;vL_out];
    
    diff = max(abs(r-r_out));
    cc = cc+1;
    disp(diff)
    disp(cc)
    
    vH = ww*vH_out + (1-ww)*vH;
    vL = ww*vL_out + (1-ww)*vL;
end

%% Simulations

z_set = [m.zL;m.zH];

z_mat=ones(m.S,m.T); % simulated productivities: dim1:s;dim2:t
for s=1:m.S
    z_mat(s,1)=z_set(m.z0_index(s));
    for t=2:m.T
        u = m.u_mat(s,t-1);
        if z_mat(s,t-1)==m.zL
            if u <=m.piLL
                z_mat(s,t) = m.zL;
            else
                z_mat(s,t) = m.zH;
            end
        else
            if u <=m.piHH
                z_mat(s,t) = m.zH;
            else
                z_mat(s,t) = m.zL;
            end
        end  
    end
end

k_mat_sim = ones(m.S,m.T);
k_mat_sim(:,1) = m.k0_vec;
d_mat_sim = ones(m.S,m.T-1);

for s=1:m.S
    for t=1:m.T-1
        if z_mat(s,t) == m.zL
            k_mat_sim(s,t+1)=m.k(jL(m.k==k_mat_sim(s,t)));
            i = k_mat_sim(s,t+1)-(1-m.delta)*k_mat_sim(s,t);
            d_mat_sim(s,t)= z_mat(s,t)*k_mat_sim(s,t)^m.alpha-i-m.c/2*i^2/k_mat_sim(s,t);
        else
            k_mat_sim(s,t+1)=m.k(jH(m.k==k_mat_sim(s,t)));
            i = k_mat_sim(s,t+1)-(1-m.delta)*k_mat_sim(s,t);
            d_mat_sim(s,t)= z_mat(s,t)*k_mat_sim(s,t)^m.alpha-i-m.c/2*i^2/k_mat_sim(s,t);
        end
    end
end

end

