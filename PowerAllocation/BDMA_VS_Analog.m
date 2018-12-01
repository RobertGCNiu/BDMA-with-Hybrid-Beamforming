%%%%%%%%%%----- Performance of Adaptive Channel Estimation of MmWave Channels-----%%%%%%%
% Author: Niu Guanchong
% Date: 2018/08/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes
%
%--------------------------------------------------------------------------
clear;clc;
% ----------------------------- System Parameters% -------------------------
Num_user_cluster = [32];
Num_users_all = 32;
Rate_SIR=zeros(1,length(Num_user_cluster));
Rate_hb = zeros(1,length(Num_user_cluster));
paths = [1];
Rate_Path = zeros(length(paths),1);
TX_sets = (8:16).^2;
Rate_HP_ant = zeros(1,length(TX_sets));
SNR_dB_range=-15:3:20;
Rate_p = zeros(1,length([3,4,5,6,7,8,9,10]));
all=0;

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
                for Num_paths_index=1:length(paths) %Number of channel paths
                    Num_paths = paths(Num_paths_index);
                    % ----------------------------- Simulation Parameters ---------------------
                    SNR_dB_range=-15:3:20;  % SNR in dB
                    

                    Rate_SU=zeros(1,length(SNR_dB_range)); % Will carry the single-user MIMO rate (without interference)
                    %Rate_p=zeros(1,length(SNR_dB_range));% Will carry the rate with analog-only beamsteering
                    Rate_BS_BDMA = zeros(1,length(SNR_dB_range));
                     Rate_HP_cl = zeros(1,length(SNR_dB_range));
                    Rate_HP_fzf = zeros(1,length(SNR_dB_range));
                    Rate_HP_schedule = zeros(1,length(SNR_dB_range));
                    Rate_HP_SLNR = zeros(1,length(SNR_dB_range));
   
                    
                    ITER=50; % Number of iterations
                    
                    % --------------- Simulation starts ---------------------------------------
                    for iter=1:1:ITER
                        iter
                        T = zeros(Num_users, Num_users);
                        
                        % Generate user channels
                        [H_all,a_TX_all,a_RX_all]=ULAMulPath(Num_users_all,TX_ant_w,RX_ant_w,Num_paths);
                        % H is a 3-dimensional matrix, with Num_users,RX_ant,TX_ant dimensions
                        
                 %       [a_TX_schedule,a_RX_schedule, ~,~,H_schedule] = Selectusers(Num_users,Num_users_all,a_TX_all,a_RX_all,Num_paths,H_all);%select Num_users from Num_users_all
                        
                        
                        H = H_all(1:Num_users,:,:);
                        a_TX = a_TX_all(:,1:Num_users,:);
                        a_RX = a_RX_all(:,1:Num_users,:);
                        
                        [a_TX_select, a_RX_select, a_TX_select_inf, a_RX_select_inf] = SelectBestBeam(Num_users,a_TX,a_RX,Num_paths,H);
                        
                        
                        Frf_BDMA = a_TX_select;
                        Wrf_BDMA = a_RX_select;
                        
                        G_bdma=effective_H(H,Wrf_BDMA,Frf_BDMA);
                        
                       
                        %% Clustering
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
                        
                        Fbb_fzf=pinv(G_fzf);
                        for u =1:Num_users
                            Fbb_fzf(:,u) = Fbb_fzf(:,u)/norm(Fbb_fzf(:,u));
                        end
                        %Fbb_fzf = OffDiagonalZero(Num_RF, Num_users,Fbb_fzf);
                        %Fbb_fzf=Fbb_fzf.*kron(eye(K),ones(m_k));
                        Fbb_fzf = normalize_f(Fbb_fzf,Frf_fzf);  
                      
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %% Power Allocation. We set the priority coefficient gamm to be 1. Thus the denote will be neglected.
                         noise_user = 1/db2pow(20); % the noise is set to be 1/(20db).
                         P_set = [3,4,5,6,7,8,9,10];
                         
                         for p_index = 1:length(P_set)
                        P_max = P_set(p_index); % the power in transmitter is set to be 3 watt.
                        q = 20*rand(1,Num_users);
                        tol = 1e-2; MaxIter = 100;
                        [q_optimal, tor_optimal]= fixPointIter(q,@funcofq,G_fzf,P_max, noise_user,  tol,MaxIter);
                        
                        
                        % Spectral efficiency calculations
%                         f = zeros(Num_users, Num_users);
%                         for u = 1: Num_users
%                             int_set = 1:Num_users;
%                             int_set(u) = [];
%                             inner = 0;
%                             for l = int_set
%                                 inner = inner+q_optimal(l)/Num_users*G_fzf(:,l)*G_fzf(:,l)';
%                             end
%                             f(:, u) =  (inner+noise_user*eye(Num_users))^(-1)*G_fzf(:,u);
%                             T(u,u) = Num_users*norm(f(:,u))^2/abs(G_fzf(:, u)'*f(:,u))^2;
%                         end
%                         F = zeros(Num_users, Num_users);
%                         for k=1:Num_users
%                             for i = 1:Num_users
%                                 if k~=i
%                                 F(k,i) = 1/Num_users*abs(G_fzf(:,k)'*f(:,i))^2/norm(f(:,i))^2;
%                                 end
%                             end
%                         end
%                         p_optimal = tor_optimal * noise_user * (eye(Num_users) - tor_optimal*T*F)^(-1)*T*ones(Num_users,1);

     

                        Rate_p(p_index)=Rate_p(p_index)+log2(1+tor_optimal)/(ITER);
                         end
                            
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %  gamma = ZFBalancedSINR(G_fzf,Fbb_fzf,noise_user);
                         
                        %    all = all+gamma/ITER;

                        %Rate_HP_fzf(Num_RF_index)=Rate_HP_fzf(Num_RF_index) + RGH(G_fzf,Fbb_fzf,rho)/(ITER);
                        %Rate_HP_ant(TX_index) = Rate_HP_ant(TX_index) + RGH(G_cl,Fbb_cl,rho)/(ITER);
                        
                    end % End of ITER loop
                end
            end
        end
        
    end
end
plot(P_set, Rate_p,'*-')

