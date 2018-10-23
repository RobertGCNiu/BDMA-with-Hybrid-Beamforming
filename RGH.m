function [R] = RGH(G,F,rho)
%RGH 此处显示有关此函数的摘要
%   此处显示详细说明
s=size(F,2);
T=abs(G*F).^2+eye(s)/rho;
s=sum(T,2);
in=s-diag(T)+1/rho;
snr=s./in;
R=sum(log2(snr));
end

