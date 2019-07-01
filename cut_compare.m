
Rate_HP_mc = zeros(1,length(SNR_dB_range));
x = [1 0
2 1
3 1
4 0
5 1
6 1
7 0
8 1
9 1
10 0
11 0
12 1
13 0
14 0
15 1
16 1
17 0
18 1
19 1
20 1
21 0
22 0
23 0
24 0
25 1
26 0
27 0
28 1
29 1
30 0
31 0
32 1
33 0
34 0
35 1
36 0
37 1
38 1
39 0
40 0
41 0
42 1
43 0
44 0
45 0
46 0
47 1
48 1
49 0
50 1
];


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