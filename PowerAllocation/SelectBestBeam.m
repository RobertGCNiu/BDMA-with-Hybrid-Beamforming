%%%If the number of users is less than 12, this function can be used for
%%%exhaustive searching
function [a_TX_select,a_RX_select, a_TX_select_inf,a_RX_select_inf] = SelectBestBeam(Num_users,a_TX,a_RX,Num_paths,H)
N_TXA = size(a_TX, 1);
N_RXA = size(a_RX,1);
a_TB = zeros(N_TXA,1);
a_RB = zeros(N_RXA,1);
a_TX_select = zeros(N_TXA,Num_users);
a_RX_select = zeros(N_RXA,Num_users);
a_TX_select_inf = zeros(N_TXA, Num_users, Num_paths-1);
a_RX_select_inf= zeros(N_RXA, Num_users, Num_paths-1);
%%
for user = 1:Num_users
    H_u(:,:) = H(user,:,:);
    for path = 1:Num_paths
        a_TB =  a_TX(:,user,path);
        a_RB =  a_RX(:,user,path);
    %    if ~(ismemeber(user,USDU)&&ismemeber(path==USDP))
     %       end
      k_inf = 1;
   %%
   %Calculate all the interference path to select the best beam
        for user_inf = 1:Num_users
            for path_inf = 1:Num_paths
%                 if user_inf~=user || path_inf~=path
                if user_inf~=user
                    a_TBInter(:, k_inf) = a_TX(:,user_inf,path_inf);
                    k_inf = k_inf+1;
                end
            end
        end
   %%     
        Signal = abs(a_RB'*H_u*a_TB)^2;
        Interference = sum(abs(a_RB'*H_u*a_TBInter).^2);
        SINR(user,path) = Signal/Interference;
    end
end
%%
[~, max_path] = max(SINR,[],2);
for user = 1:Num_users
        a_TX_select(:,user) = a_TX(:,user,max_path(user));
        a_RX_select(:,user) = a_RX(:,user,max_path(user));
        a_TX_select_inf(:,user,:) = a_TX(:,user,[1:max_path(user)-1 max_path(user)+1:end]);
        a_RX_select_inf(:,user,:) = a_RX(:,user,[1:max_path(user)-1 max_path(user)+1:end]);
end

end

