clear
clc
close all

ObjectiveFunction = @func_obj;
lb = -pi; % Lower bound of x
ub = pi;  % Upper bound of x
x0 = 0;   % Initial guess
rng default % For reproducibility
[x,fval,exitFlag,output] = simulannealbnd(ObjectiveFunction,x0,lb,ub)