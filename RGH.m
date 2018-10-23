function [R] = RGH(G,F,rho)
%RGH �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
s=size(F,2);
T=abs(G*F).^2+eye(s)/rho;
s=sum(T,2);
in=s-diag(T)+1/rho;
snr=s./in;
R=sum(log2(snr));
end

