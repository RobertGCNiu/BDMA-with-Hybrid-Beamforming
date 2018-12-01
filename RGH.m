function [R, users_rate, snr] = RGH(G,F,rho)

s=size(F,2);
T=abs(G*F).^2+eye(s)/rho;
s=sum(T,2);
in=s-diag(T)+1/rho;
snr=s./in;
users_rate = log2(snr);
R=sum(log2(snr));
end

