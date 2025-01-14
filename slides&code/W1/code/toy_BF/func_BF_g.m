function gy = func_BF_g(x)
% Jacobian matrix of Burden and Faire function

gy = ones(length(x),length(x));

gy(1,1) = 3;
gy(2,1) = 2*x(1);
gy(3,1) = -x(2)*exp(-x(1)*x(2));
gy(1,2) = x(3)*sin(x(2)*x(3));
gy(2,2) = -162*(x(2)+0.1);
gy(3,2) = -x(1)*exp(-x(1)*x(3));
gy(1,3) = x(2)*sin(x(2)*x(3));
gy(2,3) = cos(x(3));
gy(3,3) = 20;

end

