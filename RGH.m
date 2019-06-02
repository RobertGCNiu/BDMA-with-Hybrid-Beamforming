function [R, users_rate, signal] = RGH(G,F,rho,Num_users,lowerbound)
if nargin<=3
s=size(F,2);
sig = abs(G*F).^2;
signal = sum(sig,2);
T=abs(G*F).^2+eye(s)/rho;
s=sum(T,2);
in=s-diag(T)+1/rho;
snr=s./in;
users_rate = log2(snr);
R=sum(log2(snr));
end

if nargin==4
s=size(F,2);
sig = abs(G*F).^2;
signal = sum(sig,2);
T=abs(G*F).^2+eye(s)/rho;
s=sum(T,2);
in=0+1/rho;
snr=s./in;
users_rate = log2(snr);
R=sum(log2(snr));
end

if nargin==5
    R = 0;
for u = 1:Num_users    
sig = abs(G(u,:)*F(:,u)).^2;
in=sig*(Num_users-1)*0.015+1/rho;
snr=sig./in;
R = R+ log2(snr);
end
end

end

