clear
close all
clc

%% Test function BF
x0 = [0.5;0;-0.5];
y = func_BF(x0);
disp(y)

%% Simple Iteration
x = x0;
diff = 1;
tol = 1e-10;
ww = 0.1;
cc = 0;

while diff > tol && cc < 10000
    x_new = x - func_BF(x);
    diff = max(abs(x-x_new));
    cc = cc+1;
    disp(diff)
    disp(cc)
    
    x = ww*x_new+(1-ww)*x;
end

x_iter = x;
disp(x_iter)
disp(func_BF(x_iter))

%% fsolve

options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e6,'MaxIter',1000,'TolFun',1e-10,'TolX',1e-10);
[x,fval] = fsolve(@(x)func_BF(x),x0,options);

x_fsolve = x;
disp([x_fsolve,x_iter])
disp([fval,func_BF(x_iter)])

%% Newton's method
x = x0;
diff = 1;
tol = 1e-10;
ww = 0.5;
cc = 0;

while diff > tol && cc < 10000
    x_new = x - func_BF_g(x)\func_BF(x);
    diff = max(abs(x-x_new));
    cc = cc+1;
    disp(diff)
    disp(cc)
    
    x = ww*x_new+(1-ww)*x;
end

x_newton = x;
disp([x_newton,x_fsolve])

%% Quasi-Newton: Broyden's method
diff = 1;
tol = 1e-10;
ww = 0.1;
cc = 0;

x = x0;
B = eye(length(x0));

while diff > tol && cc < 10000
    x_new = x - B*func_BF(x);
    h = x_new - x;
    y = func_BF(x_new) - func_BF(x);
    B_new = B+(h'*B*y).^(-1)*(h-B*y)*h'*B;
    
    diff = max(abs(x-x_new));
    cc = cc+1;
    disp(diff)
    disp(cc)
    
    x = ww*x_new+(1-ww)*x;
    B = ww*B_new+(1-ww)*B;
end

x_broyden = x;

disp([x_iter,x_fsolve,x_newton,x_broyden])

disp(func_BF_g(x_broyden))
disp(B_new)
disp(B_new*y-h)



                   