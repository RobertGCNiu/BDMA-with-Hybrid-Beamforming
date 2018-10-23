%%%%%%%%%%----- Performance of Adaptive Channel Estimation of MmWave Channels-----%%%%%%%
% Author: Niu Guanchong
% Date: 2018/08/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes
%
%--------------------------------------------------------------------------
clear;clc;close;
% ----------------------------- System Parameters% -------------------------
Num_user_cluster = [16];
Num_users_all = 40;
Rate_SIR=zeros(1,length(Num_user_cluster));
Rate_hb = zeros(1,length(Num_user_cluster));
paths = [1];
Rate_Path = zeros(length(paths),1);
for Num_users_index=1:length(Num_user_cluster) % Number of users
    Num_users = Num_user_cluster(Num_users_index);
    for TX_ant=144 %Number of UPA TX antennas
        TX_ant_w=sqrt(TX_ant); % width
        TX_ant_h=sqrt(TX_ant); % hieght
        ind_TX_w=reshape(repmat([0:1:TX_ant_w-1],TX_ant_h,1),1,TX_ant_w*TX_ant_h);
        ind_TX_h=repmat([0:1:TX_ant_h-1],1,TX_ant_w);
        
        RX_ant=64; %Number of UPA RX antennas
        RX_ant_w=sqrt(RX_ant); % width
        RX_ant_h=sqrt(RX_ant); % hieght
        ind_RX_w=reshape(repmat([0:1:RX_ant_w-1],RX_ant_h,1),1,RX_ant_w*RX_ant_h);
        ind_RX_h=repmat([0:1:RX_ant_h-1],1,RX_ant_w);
        for k_cluster = [8]
            Num_group = Num_users/k_cluster;
            % ----------------------------- Channel Parameters ------------------------
            for Num_paths_index=1:length(paths) %Number of channel paths
                Num_paths = paths(Num_paths_index);
                % ----------------------------- Simulation Parameters ---------------------
                SNR_dB_range=-15:3:30;  % SNR in dB
                Rate_SU=zeros(1,length(SNR_dB_range)); % Will carry the single-user MIMO rate (without interference)
                Rate_BS=zeros(1,length(SNR_dB_range));% Will carry the rate with analog-only beamsteering
                Rate_BS_BDMA = zeros(1,length(SNR_dB_range));
                Rate_HP_cl = zeros(1,length(SNR_dB_range));
                Rate_HP_fzf = zeros(1,length(SNR_dB_range));
                Rate_HP_schedule = zeros(1,length(SNR_dB_range));
                Rate_HP_SLNR = zeros(1,length(SNR_dB_range));
                
                ITER=10; % Number of iterations
                
                % --------------- Simulation starts ---------------------------------------
                for iter=1:ITER
                    T = zeros(Num_users, Num_users);
                    
                    % Generate user channels
                    [H_all,a_TX_all,a_RX_all]=ULAMulPath(Num_users_all,TX_ant_w,RX_ant_w,Num_paths);
                    % H is a 3-dimensional matrix, with Num_users,RX_ant,TX_ant dimensions
                    
                    [a_TX_schedule,a_RX_schedule, ~,~,H_schedule] = Selectusers(Num_users,Num_users_all,a_TX_all,a_RX_all,Num_paths,H_all);%select Num_users from Num_users_all
                    
                    
                    H = H_all(1:Num_users,:,:);
                    a_TX = a_TX_all(:,1:Num_users,:);
                    a_RX = a_RX_all(:,1:Num_users,:);
                    
                    [a_TX_select, a_RX_select, a_TX_select_inf, a_RX_select_inf] = SelectBestBeam(Num_users,a_TX,a_RX,Num_paths,H);
                    
                    
                    Frf_BDMA = a_TX_select;
                    Wrf_BDMA = a_RX_select;
                    
                   G_bdma=effective_H(H,Wrf_BDMA,Frf_BDMA);
                    
                    Frf_fzf=a_TX_select;
                    Wrf_fzf=a_RX_select;
                    
                    % Constructin the effective channels
                    
                    G_fzf=effective_H(H,Wrf_fzf,Frf_fzf);
                    
                    % Baseband zero-forcing precoding
                    
                    Fbb_fzf=pinv(G_fzf);
                    Fbb_fzf = normalize_f(Fbb_fzf,Frf_fzf);
                    
                    
                    %Schedule Selection
                    G_schedule = effective_H(H_schedule,a_RX_schedule,a_TX_schedule);
                    
                    % Baseband zero-forcing precoding
                    %   Fbb_cl = eye(Num_users);
                    Fbb_schedule=pinv(G_schedule);
                    Fbb_schedule = normalize_f(Fbb_schedule,a_TX_schedule);
                    
                    [Wrf_cl, Frf_cl, H_cl]= greedySelection(a_TX_select, a_RX_select,k_cluster,H);
                    % Constructin the effective channels
                    G_cl = effective_H(H_cl,Wrf_cl,Frf_cl);
                    
                    % Baseband zero-forcing precoding
                    %   Fbb_cl = eye(Num_users);
                    Fbb_cl=pinv(G_cl);
                    Fbb_cl=Fbb_cl.*kron(eye(Num_group),ones(k_cluster));
                    Fbb_cl = normalize_f(Fbb_cl,Frf_cl);
                    
 
                    
                    
                    %% Power Allocation
                    for u = 1: Num_users
                        for col = 1: Num_users
                            T(u, col) = abs(G_cl(u,:) * Fbb_cl(:,col)).^2;
                        end
                    end
                    [sir, power] = GMtraixGeneration(T); %The power is orthonomal vector.
                    Rate_SIR(Num_users_index) = Rate_SIR(Num_users_index)+ log2(1+ sir)/(ITER);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Leakage-Based Hybrid Beamforming Design for Downlink Multiuser mmWave MIMO Systems
                    
                    
                    % Spectral efficiency calculations
                    SNR_index=0;
                    for SNR_dB=SNR_dB_range
                        SNR_index=SNR_index+1;
                        rho=10^(.1*SNR_dB)/Num_users; % SNR value
                        for u=1:Num_users
                            Int_set=1:Num_users; % interference index
                            Int_set(u)=[];
                            C = rho*G_fzf(u,:)'*G_fzf(u,:);
                            D = eye(Num_users)+rho*(G_fzf(Int_set,:)'*G_fzf(Int_set,:));
                            [Vector_C, lamda_C]= eigs((D^(-1)*C));
                            F_slnr(:,u) =Vector_C(:,1);
                        end
                        F_slnr=normalize_f(F_slnr,Frf_fzf);
                        
                        
                        interval=1:k_cluster:Num_users+1;
                        for u=1:Num_users
                            Int_set=1:Num_users; % interference index
                            Int_set(u)=[];
                            G_index=sum(u>=interval);
                            G_replace=G_cl(:,interval(G_index):interval(G_index+1)-1);
                            C = rho*G_replace(u,:)'*G_replace(u,:);
                            D = eye(k_cluster)+rho*(G_replace(Int_set,:)'*G_replace(Int_set,:));
                            [Vector_C, lamda_C]= eigs((D^(-1)*C));
                            f_temp=zeros(Num_users,1);
                            f_temp(interval(G_index):interval(G_index+1)-1)=Vector_C(:,1);
                            Fbb_slnr(:,u) =f_temp;
                        end
                        Fbb_slnr=normalize_f(Fbb_slnr,Frf_cl);
                        
                        
                        
                        
                        for u=1:Num_users
                            Int_set=1:Num_users; % interference index
                            Int_set(u)=[];
                            Channel=zeros(RX_ant,TX_ant);
                            Channel(:,:)= H(u,:,:);
                            [U_channel S_channel V_channel]=svd(Channel);
                            Channel_cl(:,:) =H_cl(u,:,:);
                            Rate_SU(SNR_index)=Rate_SU(SNR_index)+log2(1+rho*S_channel(1,1)^2)/(Num_users*ITER);
                        end
                        
                        
                        G_BDMA=effective_H(H,Wrf_BDMA,Frf_BDMA);
                        Rate_BS_BDMA(SNR_index)=Rate_BS_BDMA(SNR_index)+RGH(G_BDMA,eye(size(G_BDMA)),rho)/(Num_users*ITER);
                        Rate_HP_cl(SNR_index)=Rate_HP_cl(SNR_index)+RGH(G_cl,Fbb_cl,rho)/(Num_users*ITER);
                        Rate_HP_fzf(SNR_index) = Rate_HP_fzf(SNR_index)+RGH(G_fzf,Fbb_fzf,rho)/(Num_users*ITER);
                        Rate_HP_schedule(SNR_index) = Rate_HP_schedule(SNR_index) + RGH(G_schedule,Fbb_schedule,rho)/(Num_users*ITER);
                        Rate_HP_SLNR(SNR_index) = Rate_HP_SLNR(SNR_index) + RGH(G_cl,Fbb_slnr,rho)/(Num_users*ITER);
                        
                        
                    end % End of SNR loop
                    
                end % End of ITER loop
                
            end
            
        end
    end
    
end


plot(SNR_dB_range,Rate_SU,'-v');
hold on; plot(SNR_dB_range,Rate_BS_BDMA);
plot(SNR_dB_range,Rate_HP_fzf,'--o');
plot(SNR_dB_range,Rate_HP_cl,'--');
plot(SNR_dB_range,Rate_HP_schedule,'-s');
plot(SNR_dB_range,Rate_HP_SLNR,'-');

legend('signal user','BDMA','full-zf','group','schedule','SLNR')
