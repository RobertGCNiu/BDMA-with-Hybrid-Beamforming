
Rate_HP_mc = zeros(1,length(SNR_dB_range));

system('python maxcut.py test.txt text.txt');
data = load('text.txt');
x  = data(2:end,:);
[~,used_total] = sort(x(:,2));

Frf_mc = a_TX_select(:,used_total);
Wrf_mc =a_RX_select(:,used_total);
H_mc = H(used_total,:,:);
G_mc = effective_H(H_mc,Wrf_mc,Frf_mc);



Fbb_mc = [];
for k = 1:K
    Fbb_k = pinv(G_mc(1+m_k*(k-1):m_k*k, 1+m_k*(k-1):m_k*k));
    Fbb_mc = blkdiag(Fbb_mc, Fbb_k);
end

Fbb_mc = normalize_f(Fbb_mc,Frf_mc);
SNR_index=0;
for SNR_dB=SNR_dB_range
    rho=db2pow(SNR_dB)/Num_users; 
    SNR_index=SNR_index+1;
    Rate_HP_mc(SNR_index)=Rate_HP_mc(SNR_index) + RGH(G_mc,Fbb_mc,rho)/(ITER);
end

hold on
plot(SNR_dB_range,Rate_HP_mc,'--o','linewidth',1.5);
plot(SNR_dB_range,Rate_HP_cl,'--','linewidth',1.5);
legend('mc','greedy')