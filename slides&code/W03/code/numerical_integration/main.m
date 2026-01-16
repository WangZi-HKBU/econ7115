clear
close all
clc

x = 0.5;
y = -0.5;

func_int = @(x,y) x.^2 + y.^2+2*x.*y;

S1 = quad2d(func_int,0,1,0,2)
S2 = integral2(func_int,0,1,0,2)

