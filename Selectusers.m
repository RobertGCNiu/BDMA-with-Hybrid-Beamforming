%%%If the number of users is less than 12, this function can be used for
%%%exhaustive searching
% function [a_TX_select,a_RX_select, H_selected] = Selectusers(Num_users_all, Num_users,a_TX,a_RX,Num_paths,H)
% 
% % a_TX_select = a_TX(:,1:Num_users,:);
% % a_RX_select = a_RX(:,1:Num_users,:);
% % H_selected  = H(1:Num_users,:,:);
% 
% %% find the first beam with largest channel gain
% alpha_u = zeros(Num_users_all,Num_paths);
% for u = 1:Num_users_all
%     H_u(:,:) = H(u,:,:);
%     for p = 1:Num_paths
%         alpha_u(u,p) = abs(a_RX(:,u,p)'*H_u*a_TX(:,u,p));
%     end
% end
% [max_u_all, max_u] = max(alpha_u,[],1);
% [~, max_p] = max(max_u_all);
% max_user = max_u(max_p);
% %% The first user with largest channel gain is selected
% a_TX_select(:,1) = a_TX(:,max_user,max_p);
% a_RX_select(:,1) = a_RX(:,max_user,max_p);
% H_selected(1,:,:) = H(max_user,:,:);
% new_add = 1;
% while new_add<Num_users
%     new_add= new_add+1;
%     a_TX(:,max_user,:) = [];
%     H(max_user,:,:) =[]; 
%     a_RX(:,max_user,:)=[];
%     sinr = zeros(Num_users_all,Num_paths);
%     for u_try= 1:(Num_users_all-new_add+1)
%         for path = 1:Num_paths
%             a_TX_select(:,new_add) = a_TX(:,u_try,path);
%             a_RX_select(:,new_add) = a_RX(:,u_try,path);
%             H_selected(new_add,:,:) = H(u_try,:,:);
%             for user = 1:size(a_TX_select,2)
%                 H_u(:,:) = H_selected(user,:,:);
%                 user_int = [];
%                 for int = 1:size(a_TX_select,2)
%                     if int ~=user
%                         user_int = [user_int int];
%                     end
%                 end
%                 sinr(u_try,path) =sinr(u_try,path) + (abs(a_RX_select(:,user)'*H_u*a_TX_select(:,user))^2)/sum(abs(a_RX_select(:,user)'*H_u*a_TX_select(:,user_int)).^2);
%             end
%         end
%     end
% %     [max_u_all,max_u] = max(sinr,[],2) ;
% %     [~,max_p] = max(max_u_all);
% %     max_user = max_u(max_p);
%     [max_user,max_p] = find(sinr == max(max(sinr)));
%     a_TX_select(:,new_add) = a_TX(:,max_user,max_p);
%     a_RX_select(:,new_add) = a_RX(:,max_user,max_p);
%     H_selected(new_add,:,:) = H(max_user,:,:);
%     
% end
% 
% 
% end
function [a_TX_select,a_RX_select, a_TX_select_inf,a_RX_select_inf,H_select] = Selectusers(Num_users_seleted, Num_users,a_TX,a_RX,Num_paths,H)
N_TXA = size(a_TX, 1);
N_RXA = size(a_RX,1);
a_TB = zeros(N_TXA,1);
a_RB = zeros(N_RXA,1);
a_TX_select = zeros(N_TXA,Num_users_seleted);
a_RX_select = zeros(N_RXA,Num_users_seleted);
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
[max_SINR, max_path] = max(SINR,[],2);
[~, max_index] = sort(max_SINR,'descend');
for user = 1:Num_users_seleted
        a_TX_select(:,user) = a_TX(:,max_index(user),max_path(max_index(user)));
        a_RX_select(:,user) = a_RX(:,max_index(user),max_path(max_index(user)));
        H_select(user,:,:) = H(max_index(user),:,:);
        a_TX_select_inf(:,user,:) = a_TX(:,user,[1:max_path(user)-1 max_path(user)+1:end]);
        a_RX_select_inf(:,user,:) = a_RX(:,user,[1:max_path(user)-1 max_path(user)+1:end]);
end

end

