%%%%%%%%%%----- Performance of Adaptive Channel Estimation of MmWave Channels-----%%%%%%%
% Author: Ahmed Alkhateeb                                              
% Date: March. 1, 2016
% This code calculates the spectral efficincy achieved by the multi-user hybrid precoding
% algorithms in the paper - A. Alkhateeb, G. Leus, and R. W. Heath Jr., "Limited Feedback 
% Hybrid Precoding for Multi-User Millimeter Wave Systems," in IEEE Transactions on 
% Wireless Communications, vol. 14, no. 11, pp. 6481-6494, Nov. 2015.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes
% 
%--------------------------------------------------------------------------
clear;clc;
% ----------------------------- System Parameters -------------------------
Num_users=4; % Number of users
TX_ant=64; %Number of UPA TX antennas
TX_ant_w=sqrt(TX_ant); % width
TX_ant_h=sqrt(TX_ant); % hieght 
ind_TX_w=reshape(repmat([0:1:TX_ant_w-1],TX_ant_h,1),1,TX_ant_w*TX_ant_h);
ind_TX_h=repmat([0:1:TX_ant_h-1],1,TX_ant_w);

RX_ant=16; %Number of UPA RX antennas
RX_ant_w=sqrt(RX_ant); % width 
RX_ant_h=sqrt(RX_ant); % hieght
ind_RX_w=reshape(repmat([0:1:RX_ant_w-1],RX_ant_h,1),1,RX_ant_w*RX_ant_h);
ind_RX_h=repmat([0:1:RX_ant_h-1],1,RX_ant_w);

% ----------------------------- Channel Parameters ------------------------
Num_paths=4; %Number of channel paths

% ----------------------------- Simulation Parameters ---------------------
SNR_dB_range=[-20:3:10 40 100];  % SNR in dB
Rate_SU=zeros(1,length(SNR_dB_range)); % Will carry the single-user MIMO rate (without interference)
Rate_LB=zeros(1,length(SNR_dB_range));% Will carry the lower bound values
Rate_BS=zeros(1,length(SNR_dB_range));% Will carry the rate with analog-only beamsteering
Rate_HP=zeros(1,length(SNR_dB_range)); % Will carry the rate of the proposed algorithm (with analog 
% and zero-forcing digital precoding)

ITER=500; % Number of iterations
    
% --------------- Simulation starts ---------------------------------------
for iter=1:1:ITER
    % Generate user channels 
    [H,a_TX,a_RX]=generate_channels(Num_users,TX_ant_w,TX_ant_h,RX_ant_w,RX_ant_h,Num_paths); 
    % H is a 3-dimensional matrix, with Num_users,RX_ant,TX_ant dimensions

    % Stage 1 of the proposed algorithm (Analog precoding)
    Frf=zeros(TX_ant,Num_users); % BS RF precoders 
    Wrf=zeros(RX_ant,Num_users); % MS RF precoders 
    
    for u=1:1:Num_users
        Frf(:,u)=a_TX(:,u);
        Wrf(:,u)=a_RX(:,u);
    end      
    
    % Constructin the effective channels
    for u=1:1:Num_users
        Channel=zeros(RX_ant,TX_ant);
        Channel(:,:)= H(u,:,:);
        He(u,:)=Wrf(:,u)'*Channel*Frf ;    % Effective channels
    end
    
    % Baseband zero-forcing precoding
    Fbb=He'*(He*He')^(-1);   
    for u=1:1:Num_users % Normalization of the hybrid precoders
        Fbb(:,u)=Fbb(:,u)/sqrt((Frf*Fbb(:,u))'*(Frf*Fbb(:,u)));
    end

     % For the lower bound 
        [Us Ss Vs]=svd(Frf);
        s_min=(min(diag(Ss)))^2;
        s_max=(max(diag(Ss)))^2;
        G_factor=4/(s_max/s_min+s_min/s_max+2);  
      
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
            
            % Single-user rate
            Rate_SU(count)=Rate_SU(count)+log2(1+SNR*S_channel(1,1)^2)/(Num_users*ITER);
            
            % Analog-only beamforming
            SINR_BS=(SNR*(abs(Wrf(:,u)'*Channel*Frf(:,u)).^2))/(SNR*sum((abs(Wrf(:,u)'*Channel*Frf(:,Int_set)).^2))+1);
            Rate_BS(count)=Rate_BS(count)+log2(1+SINR_BS)/(Num_users*ITER);
            
            % Derived lower bound
            Rate_LB(count)=Rate_LB(count)+log2(1+SNR*S_channel(1,1)^2*G_factor)/(Num_users*ITER);
        end
    
        % Hybrid Precoding
        Rate_HP(count)=Rate_HP(count)+log2(det(eye(Num_users)+SNR*(He*(Fbb*Fbb')*He')))/(Num_users*ITER);
       
    end % End of SNR loop
end % End of ITER loop

%Plotting the spectral efficiencies
    plot(SNR_dB_range,Rate_SU,'-v','LineWidth',1.5);
    hold on; plot(SNR_dB_range,Rate_HP,'-s','LineWidth',1.5);
if Num_paths==1
    hold on; plot(SNR_dB_range,Rate_LB,'--k','LineWidth',1.5);
    hold on; plot(SNR_dB_range,Rate_BS,'-ro','LineWidth',1.5);
    legend('Single-user (No Interference)','Proposed Hybrid Precoding','Lower Bound (Theorem 1)','Analog-only Beamsteering');
else
    hold on; plot(SNR_dB_range,Rate_BS,'-ro','LineWidth',1.5);
    legend('Single-user (No Interference)','Proposed Hybrid Precoding','Analog-only Beamsteering');
end
xlabel('SNR (dB)','FontSize',12);
ylabel('Spectral Efficiency (bps/ Hz)','FontSize',12);
grid;