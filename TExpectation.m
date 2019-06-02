clear;clc
%N_T = 144;
t=0;
for N_T = 20:20:160
syms x z
t = t+1;
N = 0:(N_T-1);
n = N;

y = @(x,z) abs( (sum(exp(1j.*n.*pi.*(sin(x)-sin(z))))) );
zz = y(x,z);
yy = @(x,z) eval(char(zz)).^2;

I(t) = 1/N_T/N_T/(4*pi^2)*double(integral2(yy,0,2*pi,0,2*pi));
end
plot(1:t, I);hold on
plot(1:t, ones(length([1:t])))