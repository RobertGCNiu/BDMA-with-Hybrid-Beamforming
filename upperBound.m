%%%%%%%%%%----- Performance of Adaptive Channel Estimation of MmWave Channels-----%%%%%%%
% Author: Niu Guanchong
% Date: 2018/08/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes
%
%--------------------------------------------------------------------------
clear;clc;
% ----------------------------- System Parameters% -------------------------
Num_user_cluster = [16];
Num_users_all = 40;
Rate_SIR=zeros(1,length(Num_user_cluster));
Rate_hb = zeros(1,length(Num_user_cluster));
paths = 1;
Rate_Path = zeros(length(paths),1);
TX_sets = (8:16).^2;
Rate_HP_ant = zeros(1,length(TX_sets));
SNR_dB_range=-10:3:20;
%SNR_dB_range = 10;
sinr_all =zeros(Num_user_cluster, 50);


for Num_users_index=1:length(Num_user_cluster) % Number of users
    Num_users = Num_user_cluster(Num_users_index);
    RF_sets = 8;
    %Rate_HP_fzf = zeros(1,length(RF_sets));
    for Num_RF_index = 1:length(RF_sets)
        Num_RF = RF_sets(Num_RF_index);
        TX_index = 0;
        for TX_ant=144  %Number of UPA TX antennas
            TX_index = TX_index+1;
            TX_ant_w=sqrt(TX_ant); % width
            TX_ant_h=sqrt(TX_ant); % hieght
            ind_TX_w=reshape(repmat([0:1:TX_ant_w-1],TX_ant_h,1),1,TX_ant_w*TX_ant_h);
            ind_TX_h=repmat([0:1:TX_ant_h-1],1,TX_ant_w);
            
            RX_ant=64; %Number of UPA RX antennas
            RX_ant_w=sqrt(RX_ant); % width
            RX_ant_h=sqrt(RX_ant); % hieght
            ind_RX_w=reshape(repmat([0:1:RX_ant_w-1],RX_ant_h,1),1,RX_ant_w*RX_ant_h);
            ind_RX_h=repmat([0:1:RX_ant_h-1],1,RX_ant_w);
            for K = [1]
                m_k = Num_users/K;
                % ----------------------------- Channel Parameters ------------------------
                    Rate_upper_bound=zeros(1,length(SNR_dB_range)); % Will carry the single-user MIMO rate (without interference)
                    Rate_BS=zeros(1,length(SNR_dB_range));% Will carry the rate with analog-only beamsteering
                    Rate_BS_BDMA = zeros(1,length(SNR_dB_range));
                     Rate_HP_cl = zeros(1,length(SNR_dB_range));
                    Rate_HP_fzf = zeros(1,length(SNR_dB_range));
                    Rate_HP_schedule = zeros(1,length(SNR_dB_range));
                    Rate_HP_SLNR = zeros(1,length(SNR_dB_range));
                    Rate_lower_bound= zeros(1,length(SNR_dB_range));
                for Num_paths_index=1:length(paths) %Number of channel paths
                    Num_paths = paths(Num_paths_index);
                    % ----------------------------- Simulation Parameters ---------------------
                    %SNR_dB_range=-15:3:30;  % SNR in dB
                    
                    
                    ITER=50; % Number of iterations
                    Fbb_all = zeros(Num_users, Num_users, ITER);
                    Fbb_all_fzf = zeros(Num_users, Num_users, ITER);
                    % --------------- Simulation starts ---------------------------------------
                    for iter=1:1:ITER

                        T = zeros(Num_users, Num_users);
                        
                        % Generate user channels
                        [H_all,a_TX_all,a_RX_all]=ULAMulPath(Num_users_all,TX_ant_w,RX_ant_w,Num_paths);
                        % H is a 3-dimensional matrix, with Num_users,RX_ant,TX_ant dimensions
                        
                       % [a_TX_schedule,a_RX_schedule, ~,~,H_schedule] = Selectusers(Num_users,Num_users_all,a_TX_all,a_RX_all,Num_paths,H_all);%select Num_users from Num_users_all
                        [a_TX_schedule,a_RX_schedule, H_schedule]  = Selectusers(Num_users_all, Num_users,a_TX_all,a_RX_all,Num_paths,H_all);
                        
                        H = H_all(1:Num_users,:,:);
                        a_TX = a_TX_all(:,1:Num_users,:);
                        a_RX = a_RX_all(:,1:Num_users,:);
                        
                        [a_TX_select, a_RX_select, a_TX_select_inf, a_RX_select_inf] = SelectBestBeam(Num_users,a_TX,a_RX,Num_paths,H);
                        
                        
                        Frf_BDMA = a_TX_select;
                        Wrf_BDMA = a_RX_select;
                        
                        G_bdma=effective_H(H,Wrf_BDMA,Frf_BDMA);
                        

                        
                        
                        %Schedule Selection
                        G_schedule = effective_H(H_schedule,a_RX_schedule,a_TX_schedule);
                        
                        % Baseband zero-forcing precoding
                        %   Fbb_cl = eye(Num_users);
                        Fbb_schedule=pinv(G_schedule);
                        Fbb_schedule = normalize_f(Fbb_schedule,a_TX_schedule);
                        
                        [Wrf_cl, Frf_cl, H_cl]= greedySelection(a_TX_select, a_RX_select,K,H);
                        % Constructin the effective channels
                        G_cl = effective_H(H_cl,Wrf_cl,Frf_cl);
                        
                        % Baseband zero-forcing precoding
                        %   Fbb_cl = eye(Num_users);
                        Fbb_cl_org=pinv(G_cl);
                        %Fbb_cl=Fbb_cl.*kron(eye(K),ones(m_k));
                        %Fbb_cl = OffDiagonalZero(Num_RF, Num_users,Fbb_cl);
                        Fbb_cl = [];
                        for k = 1:K
                               Fbb_k = pinv(G_cl(1+m_k*(k-1):m_k*k, 1+m_k*(k-1):m_k*k));       
                               Fbb_cl = blkdiag(Fbb_cl, Fbb_k);
                        end
                        
                        Fbb_cl = normalize_f(Fbb_cl,Frf_cl);
                        
                        %%%%%%%%%
                       %% fully zero-forcing
                       Frf_fzf=a_TX_select;
                        Wrf_fzf=a_RX_select;
                        
                        % Constructin the effective channels
                        
                        G_fzf=effective_H(H,Wrf_fzf,Frf_fzf);
                        
                        % Baseband zero-forcing precoding
                        
                        %Fbb_fzf=pinv(G_fzf);
                        %Fbb_fzf = OffDiagonalZero(Num_RF, Num_users,Fbb_fzf);
                        %Fbb_fzf=Fbb_fzf.*kron(eye(K),ones(m_k));
                        Fbb_fzf = [];
                        for k = 1:K
                            Fbb_k = pinv(G_fzf(1+m_k*(k-1):m_k*k, 1+m_k*(k-1):m_k*k));
                            Fbb_fzf = blkdiag(Fbb_fzf, Fbb_k);
                        end
                        Fbb_fzf = normalize_f(Fbb_fzf,Frf_fzf); 

                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Leakage-Based Hybrid Beamforming Design for Downlink Multiuser mmWave MIMO Systems
                        
                        
                        % Spectral efficiency calculations
                        SNR_index=0;
                        for SNR_dB=SNR_dB_range
                            SNR_index=SNR_index+1;
                            rho=db2pow(SNR_dB)/Num_users; % SNR value
                            %rho = db2pow(16)
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       [Wrf_cl_new, Frf_cl_new, H_cl_new]= select_interCluster(a_TX_select, a_RX_select,K,H,rho);
                        % Constructin the effective channels
                        G_cl_new = effective_H(H_cl_new,Wrf_cl_new,Frf_cl_new);
                        
                        % Baseband zero-forcing precoding
                        %   Fbb_cl = eye(Num_users);
                        %Fbb_cl_org=pinv(G_cl_new);
                        %Fbb_cl=Fbb_cl.*kron(eye(K),ones(m_k));
                        %Fbb_cl = OffDiagonalZero(Num_RF, Num_users,Fbb_cl);
                        Fbb_cl_new = [];
                        for k = 1:K
                               Fbb_k_new = pinv(G_cl_new(1+m_k*(k-1):m_k*k, 1+m_k*(k-1):m_k*k));       
                               Fbb_cl_new = blkdiag(Fbb_cl_new, Fbb_k_new);
                        end
                        
                        Fbb_cl_new = normalize_f(Fbb_cl_new,Frf_cl_new);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                            
                            interval=1:m_k:Num_users+1;
                            for u=1:Num_users
                                Int_set=1:Num_users; % interference index
                                Int_set(u)=[];
                                G_index=sum(u>=interval);
                                G_replace=G_cl(:,interval(G_index):interval(G_index+1)-1);
                                C = rho*G_replace(u,:)'*G_replace(u,:);
                                D = eye(m_k)+rho*(G_replace(Int_set,:)'*G_replace(Int_set,:));
                                [Vector_C, lamda_C]= eigs((D^(-1)*C));
                                f_temp=zeros(Num_users,1);
                                f_temp(interval(G_index):interval(G_index+1)-1)=Vector_C(:,1);
                                Fbb_slnr(:,u) =f_temp;
                            end
                            Fbb_slnr=normalize_f(Fbb_slnr,Frf_cl);
                                               
                            
                            for u=1:1:Num_users
                                Int_set=1:Num_users; % interference index
                                Int_set(u)=[];
                                Channel=zeros(RX_ant,TX_ant);
                                Channel(:,:)= H(u,:,:);
                                [U_channel S_channel V_channel]=svd(Channel);
  %                              Channel_cl(:,:) =H_cl(u,:,:);
                                
                                % Single-user rate
                                
                            end
                            
                            G_BDMA=effective_H(H,Wrf_BDMA,Frf_BDMA);
                            [~,~,sinr] = RGH(G_BDMA,eye(size(G_BDMA)),rho);
                             sinr_all(:,iter) = sinr ;
                            Rate_BS_BDMA(SNR_index)=Rate_BS_BDMA(SNR_index) + RGH(G_BDMA,eye(size(G_BDMA)),rho)/(ITER);
                             Rate_HP_fzf(SNR_index)=Rate_HP_fzf(SNR_index) + RGH(G_fzf, Fbb_fzf,rho,1)/(ITER);
                            Rate_upper_bound(SNR_index)=Rate_upper_bound(SNR_index) + RGH(G_fzf, eye(size(Fbb_fzf)),rho,1)/(ITER);
                            %[~,~,sinr] = RGH(G_fzf,Fbb_fzf,rho);
                           
                            Rate_HP_cl(SNR_index)=Rate_HP_cl(SNR_index) + RGH(G_cl,Fbb_cl,rho)/(ITER);
                            Rate_HP_schedule(SNR_index) = Rate_HP_schedule(SNR_index) + RGH(G_cl_new,Fbb_cl_new,rho)/(ITER);
                            Rate_HP_SLNR(SNR_index) = Rate_HP_SLNR(SNR_index) + RGH(G_cl,Fbb_slnr,rho)/(ITER);
                            Rate_lower_bound(SNR_index) = Rate_lower_bound(SNR_index) + RGH(G_cl, blkdiag(eye(size(Fbb_fzf))),rho,Num_users,1)/(ITER);
                         end % End of SNR loop
                        %Rate_HP_fzf(Num_RF_index)=Rate_HP_fzf(Num_RF_index) + RGH(G_fzf,Fbb_fzf,rho)/(ITER);
                        %Rate_HP_ant(TX_index) = Rate_HP_ant(TX_index) + RGH(G_cl,Fbb_cl,rho)/(ITER);
                        Fbb_all(:,:,iter) = abs(Fbb_cl);
                        Fbb_all_fzf(:,:,iter) = abs(Fbb_fzf);
                    end % End of ITER loop
                    
                end
            end
        end
        
    end
end
%cdf_test(sinr_all)

 %plot(RF_sets,Rate_HP_fzf,'-v');
 figure
 %plot(TX_sets,Rate_HP_ant,'-v');
% 
plot(SNR_dB_range,Rate_upper_bound,'-v','linewidth',2);
hold on; plot(SNR_dB_range,Rate_BS_BDMA,'linewidth',2);
plot(SNR_dB_range,Rate_HP_fzf,'--o','linewidth',2);
plot(SNR_dB_range,Rate_HP_cl,'--','linewidth',2);
plot(SNR_dB_range,Rate_HP_SLNR,'-','linewidth',2);
plot(SNR_dB_range,Rate_lower_bound,'-','linewidth',2);
legend('Upper Bound','BDMA','full-zf','group', 'SLNR','lowerbound')
xlabel('SNR')
ylabel('Sum-rate Spectral Efficiency(bps/Hz)')