%%%It's very difficulty to try the exhastive the possible combination since
%%%the large computation cost. Thus we use greedy method to maximize SINR
function    [a_TX_select,a_RX_select, a_TX_select_inf,a_RX_select_inf]= SelectMaxSinrToBeam(Num_users,a_TX,a_RX,Num_paths,H)
N_TXA = size(a_TX, 1);
N_RXA = size(a_RX,1);
a_TB = zeros(N_TXA,1);
a_RB = zeros(N_RXA,1);
a_TX_select = zeros(N_TXA,Num_users);
a_RX_select = zeros(N_RXA,Num_users);
a_TX_select_inf = zeros(N_TXA, Num_users, Num_paths-1);
a_RX_select_inf= zeros(N_RXA, Num_users, Num_paths-1);

for path_u = 1:Num_paths
    a_TX_select(:,1) = a_TX(:,1,path_u);
    a_RX_select(:,1) = a_RX(:,1,path_u);
        for user = 2 : Num_users
            for path_otheru = 1:Num_paths
            a_TX_try = a_TX(:, user, path_otheru);
            a_RX_try = a_RX(:, user, path_otheru);
            a_TX_select(:,1:user) = [a_TX_select(:,1:user-1) a_TX_try];
            a_RX_select(:,1:user) = [a_RX_select(:,1:user-1) a_RX_try];
                for u = 1:user
                    intf=a_TX_select;
                    intf_all = sum(abs(a_RX_select(:,u)'*H_u*intf).^2);
                    SINR_BS_fzf=(SNR*(abs(Wrf_fzf(:,u)'*Channel*Frf_fzf*Fbb_fzf(:,u)).^2))/(SNR*sum((abs(Wrf_fzf(:,u)'*Channel*Frf_fzf*Fbb_fzf(:,Int_set)).^2)+intf_all)+1);
                    Rate_HP_fzf(count)=Rate_HP_fzf(count)+log2(1+SINR_BS_fzf)/(Num_users*ITER);
                end
            end
        end
end


end

