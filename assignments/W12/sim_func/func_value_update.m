function [vH_out,vL_out,jH,jL] = func_value_update(vH,vL,m)
% Update the value function

vH_out = ones(m.J,1);
vL_out = ones(m.J,1);
jH = ones(m.J,1);
jL = ones(m.J,1);

for j=1:m.J
    kk = m.k(m.k>=(1-m.delta)*m.k(j));
    H = length(kk);
    vH_vec = vH(end-(H-1):end);
    vL_vec = vL(end-(H-1):end);
  
    i_vec  = kk-(1-m.delta)*m.k(j); %investment
    dH_vec = m.zH*m.k(j)^m.alpha-i_vec-m.c/2*i_vec.^2/m.k(j);
    dL_vec = m.zL*m.k(j)^m.alpha-i_vec-m.c/2*i_vec.^2/m.k(j);
    
    [vH_out(j),temp_jH] = max(dH_vec+m.beta*(m.piHH*vH_vec+(1-m.piHH)*vL_vec));
    [vL_out(j),temp_jL] = max(dL_vec+m.beta*(m.piLL*vL_vec+(1-m.piLL)*vH_vec));
    
    jH(j) = m.J-H+temp_jH;
    jL(j) = m.J-H+temp_jL;
end

end

