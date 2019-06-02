clear;
clc
rho = 1;
M=4;
t = 0:0.01:2;
z = sum(normrnd(0,1,[10000,2]).^2,2);
y = sum(normrnd(0,1,[10000,2*M-2]).^2,2);
x = z./(1/rho+y);
%ksdensity(x)
hold on
cdfplot(x)



f = 1/2*exp(-t/2/rho) .* (1+t).^(1-M);
%g = exp(-t/rho).*(1+t).^(-M) .* (1/rho*(1+t)+M-1);
%g_cdf = 1-exp(-t/rho)./((1+t).^(M-1));
g = 1-exp(-t/rho)./(1+t).^(M-1);
g_cdf = (1-exp(-t/2/rho).*(1+t).^(1-M));
 %plot(t,f,'b--')
plot(t,g,'r-')
plot(t,g_cdf,'k-')