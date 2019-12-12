%%%%%%%%%%----- Performance of Adaptive Channel Estimation of MmWave Channels-----%%%%%%%
% Author: Robert Niu                                        
% Date: Nov. 27, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes
% 
%--------------------------------------------------------------------------
clear;clc;
% ----------------------------- System Parameters -------------------------
Num_users_group=24; % Number of users
TX_ant=16; %Number of UPA TX antennas
TX_ant_w=sqrt(TX_ant); % width
TX_ant_h=sqrt(TX_ant); % hieght 
ind_TX_w=reshape(repmat([0:1:TX_ant_w-1],TX_ant_h,1),1,TX_ant_w*TX_ant_h);
ind_TX_h=repmat([0:1:TX_ant_h-1],1,TX_ant_w);

RX_ant=4; %Number of UPA RX antennas
RX_ant_w=sqrt(RX_ant); % width 
RX_ant_h=sqrt(RX_ant); % hieght
ind_RX_w=reshape(repmat([0:1:RX_ant_w-1],RX_ant_h,1),1,RX_ant_w*RX_ant_h);
ind_RX_h=repmat([0:1:RX_ant_h-1],1,RX_ant_w);

% ----------------------------- Channel Parameters ------------------------
Num_paths=1; %Number of channel paths
Num_groups_all = [2,3,4];
% ----------------------------- Simulation Parameters ---------------------
SNR_dB_range=-10:5:15;  % SNR in dB
Rate_BS_mkc=zeros(1,length(SNR_dB_range)); % Will carry the single-user MIMO rate (without interference)
Rate_BS_org=zeros(1,length(SNR_dB_range));% Will carry the lower bound values
Rate_BS=zeros(1,length(SNR_dB_range));% Will carry the rate with analog-only beamsteering
Rate_HP=zeros(1,length(SNR_dB_range)); % Will carry the rate of the p roposed algorithm (with analog 
% and zero-forcing digital precoding)

ITER=1; % Number of iterations

% --------------- Simulation starts ---------------------------------------
for Num_group_num = 1:length(Num_groups_all)
    Num_users = Num_users_group;
    Num_groups = Num_groups_all(Num_group_num);
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
    

    
    users_inagroup = Num_users/Num_groups;
    flag =1;
    Wrf_all = zeros(RX_ant,users_inagroup, Num_groups); 
    Frf_all = zeros(TX_ant,users_inagroup, Num_groups); 
    H_all = zeros(users_inagroup, RX_ant, TX_ant,Num_groups);
    
    
    
    
    for group = 1:Num_groups
            G_org = effective_H(H((group-1)*users_inagroup+1: group*users_inagroup,:,:),Wrf(:,(group-1)*users_inagroup+1: group*users_inagroup),Frf(:,(group-1)*users_inagroup+1: group*users_inagroup));
           SNR_index=0;
            for SNR_dB=SNR_dB_range
                SNR_index=SNR_index+1;
                rho=db2pow(SNR_dB)/users_inagroup; % SNR value
                Rate_BS_org(SNR_index)=Rate_BS_org(SNR_index) + RGH(G_org,eye(size(G_org)),rho)/(ITER);
            end % End of SNR loop
    end
    
    for group = 1:Num_groups
        Wrf_all(:,:,group) = Wrf(:,(group-1)*users_inagroup+1: group*users_inagroup);
        Frf_all(:,:,group) = Frf(:,(group-1)*users_inagroup+1: group*users_inagroup);
        H_all(:, :,:,group) = H((group-1)*users_inagroup+1: group*users_inagroup,:,:);
        
    end


    
    
    weight_or = 0;
 for group_u = 1:Num_groups
     for u = 1:users_inagroup
         for group_v = 1:Num_groups
             if group_v ~=group_u
                 H_or(:,:) = H_all(u,:,:,group_u);
             for v = 1:users_inagroup
                 weight_or= weight_or + abs(Wrf_all(:,v,group_v)'* H_or*Frf_all(:,v,group_v))^2;
             end
             end
         end
     end
 end
    
 while(flag==1)
     flag = 0;
    for group = 1:Num_groups
        for u = 1:users_inagroup
            w_ui =0;
            H_u(:,:) = H_all(u,:,:,group);
            for i_u =1: users_inagroup
                if i_u~=u
                w_ui = w_ui + abs(Wrf_all(:,i_u,group)' * H_u * Frf_all(:,i_u,group))^2;   %%calculate the WuVi
                end
            end
            for group_vl =1:Num_groups
                if group_vl~=group
                    for v = 1 : users_inagroup
                                w_vl = 0;
                                w_ul = 0;
                                w_vi = 0;
                        H_v(:,:) = H_all(v,:,:,group_vl);
                        for vl = 1: users_inagroup
                            w_ul = w_ul + abs(Wrf_all(:,u,group)' * H_v * Frf_all(:,u,group))^2;   %%calculate the WvVi
                            w_vi = w_vi + abs(Wrf_all(:,vl,group_vl)' * H_u * Frf_all(:,vl,group_vl))^2;
                            if vl~=v
                                w_vl = w_vl +  abs(Wrf_all(:,vl,group_vl)' * H_v * Frf_all(:,vl,group_vl))^2; %% calculate the WvVl
                            end
                        end
                        if w_ui + w_vl > w_ul + w_vi
                            flag = 1;
                            Wrf_all_temp= Wrf_all(:,vl,group_vl);
                            Frf_all_temp = Frf_all(:,vl,group_vl);
                            H_temp = H_v;
                            Wrf_all(:,vl,group_vl) = Wrf_all(:,u,group);
                            Frf_all(:,vl,group_vl) =  Frf_all(:,u,group);
                            H_all(v,:,:,group_vl) = H_all(u,:,:,group);
                            Wrf_all(:,u,group) = Wrf_all_temp;
                             Frf_all(:,u,group) =Frf_all_temp;
                            H_all(u,:,:,group) =H_temp;
                        end
                    end
                end
                
            end
        end
    end
        
 end
 
 weight_mkc = 0;
 
 for group_u = 1:Num_groups
     for u = 1:users_inagroup
         for group_v = 1:Num_groups
             if group_v ~=group_u
                 H_mkc(:,:) = H_all(u,:,:,group_u);
             for v = 1:users_inagroup
                 weight_mkc = weight_mkc + abs(Wrf_all(:,v,group_v)'* H_mkc*Frf_all(:,v,group_v))^2;
             end
             end
         end
     end
 end
    
 
    
end % End of ITER loop
weight_or_all(Num_group_num) = log2(weight_or);
weight_mkc_all(Num_group_num)= log2(weight_mkc);

       for group = 1:Num_groups
            Frf_mkc(:,:) = Frf_all(:,:,group);
            Wrf_mkc(:,:) = Wrf_all(:,:,group);
            H_mkc_tem(:,:,:) = H_all(:,:,:,group);

            G_mkc = effective_H(H_mkc_tem,Wrf_mkc,Frf_mkc);
            SNR_index=0;
            for SNR_dB=SNR_dB_range
                SNR_index=SNR_index+1;
                rho=db2pow(SNR_dB)/users_inagroup; % SNR value

                Rate_BS_mkc(SNR_index)=Rate_BS_mkc(SNR_index) + RGH(G_mkc,eye(size(G_mkc)),rho)/(ITER);
            end % End of SNR loop
        end

end

