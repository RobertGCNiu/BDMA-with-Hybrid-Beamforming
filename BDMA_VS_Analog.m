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
Num_users_all = 50;
     Rate_SIR=zeros(1,length(Num_user_cluster));
     Rate_hb = zeros(1,length(Num_user_cluster));
     paths = [4];
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
for k_cluster = [4]
Num_group = Num_users/k_cluster;
% ----------------------------- Channel Parameters ------------------------
for Num_paths_index=1:length(paths) %Number of channel paths
Num_paths = paths(Num_paths_index);
% ----------------------------- Simulation Parameters ---------------------
SNR_dB_range=-15:3:30;  % SNR in dB
Rate_SU=zeros(1,length(SNR_dB_range)); % Will carry the single-user MIMO rate (without interference)
Rate_LB=zeros(1,length(SNR_dB_range));% Will carry the lower bound values
Rate_BS=zeros(1,length(SNR_dB_range));% Will carry the rate with analog-only beamsteering
Rate_HP=zeros(1,length(SNR_dB_range)); % Will carry the rate of the proposed algorithm (with analog 
% and zero-forcing digital precoding)
Rate_BS_BDMA = zeros(1,length(SNR_dB_range));
Rate_HP_cl = zeros(1,length(SNR_dB_range));
Rate_HP_fzf = zeros(1,length(SNR_dB_range));
Rate_HP_schedule = zeros(1,length(SNR_dB_range));

ITER=50; % Number of iterations

% --------------- Simulation starts ---------------------------------------
for iter=1:1:ITER
     T = zeros(Num_users, Num_users);

    % Generate user channels 
    [H_all,a_TX_all,a_RX_all]=ULAMulPath(Num_users_all,TX_ant_w,RX_ant_w,Num_paths); 
    % H is a 3-dimensional matrix, with Num_users,RX_ant,TX_ant dimensions
    
   [a_TX_schedule,a_RX_schedule, H_schedule] = Selectusers(Num_users_all, Num_users,a_TX_all,a_RX_all,Num_paths,H_all);%select Num_users from Num_users_all
    
   
   H = H_all(1:Num_users,:,:);
   a_TX = a_TX_all(:,1:Num_users,:);
   a_RX = a_RX_all(:,1:Num_users,:);
   
    [a_TX_select, a_RX_select, a_TX_select_inf, a_RX_select_inf] = SelectBestBeam(Num_users,a_TX,a_RX,Num_paths,H);
   
    Frf=zeros(TX_ant,Num_users); % BS RF precoders 
    Wrf=zeros(RX_ant,Num_users); % MS RF precoders 
    
    for u=1:1:Num_users
        Frf(:,u)=a_TX(:,u,1);
        Wrf(:,u)=a_RX(:,u,1);
    end   
    
    
       Frf_BDMA = a_TX_select;
       Wrf_BDMA = a_RX_select;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Frf_fzf=zeros(TX_ant,Num_users); % BS RF precoders 
    %Wrf_fzf=zeros(RX_ant,Num_users); % MS RF precoders 
        for u=1:1:Num_users
            Frf_fzf(:,u)=a_TX_select(:,u);
            Wrf_fzf(:,u)=a_RX_select(:,u);
        end      
    %Wrf_fzf= findTheBestResponseVector(a_RX_select);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
     % Constructin the effective channels
     He_fzf = zeros(Num_users, Num_users);
    for u=1:1:Num_users
        Channel=zeros(RX_ant,TX_ant);
        Channel(:,:)= H(u,:,:);
        He_fzf(u,:)=Wrf_fzf(:,u)'*Channel*Frf_fzf ;    % Effective channels
    end
 
    % Baseband zero-forcing precoding
    Fbb_fzf=He_fzf'*(He_fzf*He_fzf')^(-1);   

    for u=1:1:Num_users % Normalization of the hybrid precoders
        Fbb_fzf(:,u)=Fbb_fzf(:,u)/sqrt((Frf_fzf*Fbb_fzf(:,u))'*(Frf_fzf*Fbb_fzf(:,u)));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Schedule Selection
      He_schedule = zeros(Num_users,Num_users);
    for u=1:1:Num_users
        Channel_cl=zeros(RX_ant,TX_ant);
        Channel_cl(:,:)= H_schedule(u,:,:);
        He_schedule(u,:)=a_RX_schedule(:,u)'*Channel_cl*a_TX_schedule ;    % Effective channels
    end
    
    % Baseband zero-forcing precoding
%   Fbb_cl = eye(Num_users);
    Fbb_schedule=He_schedule'*(He_schedule*He_schedule')^(-1);
    
%    for k = 0:k_cluster-1
%         for col = (1+Num_group*(k+1)) : Num_users
%             for row =  (1+Num_group*k) : (Num_group+Num_group*k)
%                     Fbb_schedule(col,row) = 0;
%             end
%         end
%         for row = (1+Num_group*(k+1)) : Num_users
%             for col =  (1+Num_group*k) : (Num_group+Num_group*k)
%                     Fbb_schedule(col,row) = 0;
%             end
%         end
%    end
    
    for u=1:1:Num_users % Normalization of the hybrid precoders
        Fbb_schedule(:,u)=Fbb_schedule(:,u)/sqrt((a_TX_schedule*Fbb_schedule(:,u))'*(a_TX_schedule*Fbb_schedule(:,u)));
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     [Wrf_cl, Frf_cl, H_cl]= greedySelection(a_TX_select, a_RX_select,k_cluster,H);
        % Constructin the effective channels
        He_cl = zeros(Num_users,Num_users);
    for u=1:1:Num_users
        Channel_cl=zeros(RX_ant,TX_ant);
        Channel_cl(:,:)= H_cl(u,:,:);
        He_cl(u,:)=Wrf_cl(:,u)'*Channel_cl*Frf_cl ;    % Effective channels
    end
    
    % Baseband zero-forcing precoding
%   Fbb_cl = eye(Num_users);
    Fbb_cl=He_cl'*(He_cl*He_cl')^(-1);
    
   for k = 0:k_cluster-1
        for col = (1+Num_group*(k+1)) : Num_users
            for row =  (1+Num_group*k) : (Num_group+Num_group*k)
                    Fbb_cl(col,row) = 0;
            end
        end
        for row = (1+Num_group*(k+1)) : Num_users
            for col =  (1+Num_group*k) : (Num_group+Num_group*k)
                    Fbb_cl(col,row) = 0;
            end
        end
   end
    
    for u=1:1:Num_users % Normalization of the hybrid precoders
        Fbb_cl(:,u)=Fbb_cl(:,u)/sqrt((Frf_cl*Fbb_cl(:,u))'*(Frf_cl*Fbb_cl(:,u)));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CVX method

      He_cl_cvx = zeros(Num_users,Num_users);
    for u=1:1:Num_users
        Channel_cl_cvx=zeros(RX_ant,TX_ant);
        Channel_cl_cvx(:,:)= H_cl(u,:,:);
        He_cl_cvx(u,:)=Wrf_cl(:,u)'*Channel_cl_cvx*Frf_cl ;    % Effective channels
    end
     Fbb_cl_cvx=He_cl_cvx'*(He_cl_cvx*He_cl_cvx')^(-1);
    for k = 0:k_cluster-1
        for col = (1+Num_group*(k+1)) : Num_users
            for row =  (1+Num_group*k) : (Num_group+Num_group*k)
                    Fbb_cl_cvx(col,row) = 0;
            end
        end
        for row = (1+Num_group*(k+1)) : Num_users
            for col =  (1+Num_group*k) : (Num_group+Num_group*k)
                    Fbb_cl_cvx(col,row) = 0;
            end
        end
   end
    % Baseband zero-forcing precoding
%   Fbb_cl = eye(Num_users);
%  Fbb_cl_cvx=He_cl_cvx'*(He_cl_cvx*He_cl_cvx')^(-1);
% Fbb_cl_cvx = cvxToBlockOpt(He_cl_cvx, Fbb_fzf_group, k_cluster);

   for u=1:1:Num_users % Normalization of the hybrid precoders
        Fbb_cl_cvx(:,u)=Fbb_cl_cvx(:,u)/sqrt((Frf_cl*Fbb_cl_cvx(:,u))'*(Frf_cl*Fbb_cl_cvx(:,u)));
   end

   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % For the lower bound 
        [Us Ss Vs]=svd(Frf);
        s_min=(min(diag(Ss)))^2;
        s_max=(max(diag(Ss)))^2;
        G_factor=4/(s_max/s_min+s_min/s_max+2);  
      
      
    %% Power Allocation
            for u = 1: Num_users
                for col = 1: Num_users
                    T(u, col) = abs(He_cl(u,:) * Fbb_cl(:,col)).^2;
                 end
            end
            [sir, power] = GMtraixGeneration(T); %The power is orthonomal vector.
            Rate_SIR(Num_users_index) = Rate_SIR(Num_users_index)+ log2(1+ sir)/(ITER);
        
    % Spectral efficiency calculations
    count=0;
    for SNR_dB=SNR_dB_range
        count=count+1;
        SNR=10^(.1*SNR_dB)/Num_users; % SNR value
              
        for u=1:1:Num_users
            Int_set=[]; % interference index
            for i=1:1:Num_users
                if(i~=u)
                    Int_set=[Int_set i]; 
                end
            end
            Channel=zeros(RX_ant,TX_ant);
            Channel(:,:)= H(u,:,:);
            [U_channel S_channel V_channel]=svd(Channel);
            Channel_cl(:,:) =H_cl(u,:,:);
            
            % Single-user rate
            Rate_SU(count)=Rate_SU(count)+log2(1+SNR*S_channel(1,1)^2)/(Num_users*ITER);
            
            % Analog-only beamforming
            SINR_BS=(SNR*(abs(Wrf(:,u)'*Channel*Frf(:,u)).^2))/(SNR*sum((abs(Wrf(:,u)'*Channel*Frf(:,Int_set)).^2))+1);
            Rate_BS(count)=Rate_BS(count)+log2(1+SINR_BS)/(Num_users*ITER);
          
            %%%%%%%%%%%%%%BDMA
            SINR_BS_BDMA=(SNR*(abs(Wrf_BDMA(:,u)'*Channel*Frf_BDMA(:,u)).^2))/(SNR*sum((abs(Wrf_BDMA(:,u)'*Channel*Frf_BDMA(:,Int_set)).^2))+1);
            Rate_BS_BDMA(count)=Rate_BS_BDMA(count)+log2(1+SINR_BS_BDMA)/(Num_users*ITER);
            
            
%             % Derived lower bound
%             Rate_LB(count)=Rate_LB(count)+log2(1+SNR*S_channel(1,1)^2*G_factor)/(Num_users*ITER)*Num_users*4;
             
            % Hybrid Precoding with block diagonal matrix
            SINR_BS_select=(SNR*(abs(He_cl_cvx(u,:)*Fbb_cl_cvx(:,u)).^2))/(SNR*sum((abs(He_cl_cvx(u,:)*Fbb_cl_cvx(:,Int_set)).^2))+1);
            Rate_HP(count)=Rate_HP(count)+log2(1+SINR_BS_select)/(Num_users*ITER);
            
           % Hybrid Precoding with clustering with power allocation
%             SINR_BS_cl=(SNR*power(u)*(abs(He_cl(u,:)*Fbb_cl(:,u)).^2))/(SNR*sum((abs(He_cl(u,:)*Fbb_cl(:,Int_set)).^2)*power(Int_set))+1);
%             Rate_HP_cl(count)=Rate_HP_cl(count)+log2(1+SINR_BS_cl)/(Num_users*ITER);
            %without power allocation
            SINR_BS_cl=(SNR*(abs(He_cl(u,:)*Fbb_cl(:,u)).^2))/(SNR*sum((abs(He_cl(u,:)*Fbb_cl(:,Int_set)).^2))+1);
            Rate_HP_cl(count)=Rate_HP_cl(count)+log2(1+SINR_BS_cl)/(Num_users*ITER);


            SIR_hb=(SNR*(abs(He_cl(u,:)*Fbb_cl(:,u)).^2))/(SNR*sum((abs(He_cl(u,:)*Fbb_cl(:,Int_set)).^2)));
            Rate_hb(Num_users_index)=Rate_hb(Num_users_index)+log2(1+SIR_hb)/(Num_users*ITER+1);
            %%
            
            %Hybrid Precoding with fully zero-forcing
            SINR_BS_fzf=(SNR*(abs(Wrf_fzf(:,u)'*Channel*Frf_fzf*Fbb_fzf(:,u)).^2))/(SNR*sum((abs(Wrf_fzf(:,u)'*Channel*Frf_fzf*Fbb_fzf(:,Int_set)).^2))+1);
            Rate_HP_fzf(count)=Rate_HP_fzf(count)+log2(1+SINR_BS_fzf)/(Num_users*ITER);
            
            %%%%%%%%%%%%%
            %Scheduling is considered.
            SINR_schedule=(SNR*(abs(He_schedule(u,:)*Fbb_schedule(:,u)).^2))/(SNR*sum((abs(He_schedule(u,:)*Fbb_schedule(:,Int_set)).^2))+1);
            Rate_HP_schedule(count)=Rate_HP_schedule(count)+log2(1+SINR_schedule)/(Num_users*ITER);
        end
    
        % Hybrid Precoding
        % Rate_HP(count)=Rate_HP(count)+log2(det(eye(Num_users)+SNR*(He_BDMA*(Fbb_BDMA*Fbb_BDMA')*He_BDMA')))/(Num_users*ITER);
       
    end % End of SNR loop

end % End of ITER loop


%Plotting the spectral efficiencies
  %   plot(SNR_dB_range,Rate_SU,'-v','LineWidth',1.5);
 %       hold on; plot(SNR_dB_range,Rate_HP,'LineWidth',1.5);
 %  hold on; plot(SNR_dB_range,Rate_HP,'-s','LineWidth',1.5);
% if Num_paths==1
%     hold on; plot(SNR_dB_range,Rate_LB,'--k','LineWidth',1.5);
%     hold on; plot(SNR_dB_range,Rate_BS,'-ro','LineWidth',1.5);
%     legend('Single-user (No Interference)','Proposed Hybrid Precoding','Lower Bound (Theorem 1)','Analog-only Beamsteering');
% else
%     hold on; plot(SNR_dB_range,Rate_BS,'-ro','LineWidth',1.5);
 %  hold on;  plot(SNR_dB_range,Rate_BS_BDMA,'LineWidth',1.5);
   
%     legend('Single-user (No Interference)','Proposed Hybrid Precoding','Analog-only Beamsteering','BDMA');    
% end
% hold on; plot(SNR_dB_range,Rate_HP,'LineWidth',1.5);
% hold on; plot(SNR_dB_range,Rate_HP_cl,'LineWidth',1.5);

%Rate_Path(Num_paths_index) = Rate_HP_cl;
end
%hold on; plot(SNR_dB_range,Rate_HP_cl,'--','LineWidth',1.5);
end
end
% hold on; plot(SNR_dB_range,Rate_HP_cl,'-','LineWidth',1.5);
% hold on;  plot(SNR_dB_range,Rate_HP_fzf,'--o','LineWidth',1.5);

% hold on;  plot(SNR_dB_range,Rate_HP,'--','LineWidth',1.5);
end

% plot(Num_user_cluster,Rate_hb,'-v','LineWidth',1.5);
% plot(Num_user_cluster,Rate_SIR,'-s','LineWidth',1.5);
% legend('Uniform Power', 'Balanced Power')
%plot(paths, Rate_Path,'-*','LineWidth',1.5);
plot(SNR_dB_range,Rate_SU,'-v','LineWidth',1.5);
hold on; plot(SNR_dB_range,Rate_BS_BDMA,'LineWidth',1.5);
hold on; plot(SNR_dB_range,Rate_BS,'-ro','LineWidth',1.5);
 hold on; plot(SNR_dB_range,Rate_HP,'LineWidth',1.5);
hold on;  plot(SNR_dB_range,Rate_HP_fzf,'--o','LineWidth',1.5);
hold on; plot(SNR_dB_range,Rate_HP_cl,'--b','LineWidth',1.5);
hold on; plot(SNR_dB_range,Rate_HP_schedule,'-s','LineWidth',1.5);
 legend('signal user','BDMA','analog only', 'block-cvx-ZF','full-zf','group','schedule')
% %legend('2','4','6','8','Single-user','Analog','Hybrid')
% % legend('128BDMA','Single-user','Analog Only','Hybrid','Hybrid cluster')
% % legend('5zf','5group','5','10zf','10group','10','15zf','15group','15','20zf','20group','20')

% %legend('block-cvx-ZF1','group1','block-cvx-ZF2','group2','block-cvx-ZF4','group4')
% % legend('full ZF 5 users','Delay profile 5 users', 'full ZF 10 users','Delay profile 10 users', 'full ZF 15 users','Delay profile 15 users', 'full ZF 20 users','Delay profile 20 users')
% % legend('5','5B','10','10B','15','15B','20','20B')
% xlabel('SNR (dB)','FontSize',12);
% ylabel('Spectral Efficiency (bps/ Hz)','FontSize',12);
% grid;