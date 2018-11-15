function [Wrf, Frf, H_u]= greedySelection(a_TX, a_RX,k_cluster,H,SNR)
% if nargin<5
%     SNR = 1;
% end
% 
% Num_users = size(a_TX,2);
% group_size = Num_users/k_cluster;
% user = zeros(k_cluster, group_size);
% 
% for u = 1:Num_users
%     channel(:,:) = H(u,:,:);
%     alpha_u(u) = a_RX(:,u)'*channel*a_TX(:,u);
% end
% [~,max_u1] = max(alpha_u);
% selected_users = max_u1;
% num_selected_users=length(selected_users);
% 
% k=1;
% user_in_cluster = 1;
% user(k,user_in_cluster) = max_u1;
% while(num_selected_users<Num_users)
%     if sum(user(k,:)~=0)==k_cluster
%         k=k+1;
%         user_in_cluster = 1;
%     end
%     R = zeros(Num_users,1);
%     for x = 1:Num_users
%         if ~ismember(x,selected_users)
%             tem_x = [selected_users, x];
%             for t =1:length(tem_x)
%                 if length(tem_x)<=group_size
%                     eta=0;
%                 else 
%                     eta=1;
%                 end
%                 channel(:,:) = H(tem_x(t),:,:);
%                 int_set = tem_x;
%                     if floor(t/group_size)==floor(length(tem_x)/group_size)
%                        int_set(floor(t/group_size)*group_size+1 : end) = [];
%                     else
%                        int_set(floor(t/group_size)*group_size+1 : (floor(t/group_size)+1)*group_size) = [];
%                     end
%                 R(x) =R(x) + log2(1+SNR*abs(a_RX(:,tem_x(t))'*channel*a_TX(:,tem_x(t)))^2/(SNR*sum(abs(a_RX(:,tem_x(t))'*channel*a_TX(:,int_set)).^2)*eta+1));
%             end
%         end
%     end
%     [~,selected_u]=max(R);
%     selected_users = [selected_users,selected_u];
%     num_selected_users = length(selected_users);
%     %user_in_cluster = user_in_cluster+1;
%     user(k,user_in_cluster) = selected_u;
% end
% 
% for u =1:Num_users
%     Frf(:,u) = a_TX(:,selected_users(u));
%     Wrf(:,u) = a_RX(:,selected_users(u));
%     H_u(u,:,:) = H(selected_users(u),:,:);
% end
% 
% 
% 
% end

N_TX = size(a_TX,1);
Num_users = size(a_TX,2);
users_kk = Num_users/k_cluster;
Frf(:,1) = a_TX(:,1);
used= [1];
used_total = [];
used_k = used;
for kk = 1:k_cluster   
    for k_in_cluster = 1:((Num_users/k_cluster)-1)
        residual_max  = zeros(Num_users,1);
        for user = 1:Num_users
            if ~ismember(user, used_k)
            residual_max(user) = norm(a_TX(:,used)'*a_TX(:,user));
            else
                residual_max(user) = 0;
            end
        end
        [~,select_user] = max(residual_max);
        used = [used select_user];
        used_k = [used_k select_user];
    end
    used_total = [used_total used];
   
    residual_min = zeros(Num_users,1);
    for user = 1:Num_users
        if ~ismember(user,used_k)
            residual_min(user) = norm(a_TX(:,used)'*a_TX(:,user));
        else
            residual_min(user) = 1000;
        end
    end
    [~,select_min] = min(residual_min);
    used = select_min;
    used_k = [used_k select_min];
end


Frf = a_TX(:,used_total);
Wrf = a_RX(:,used_total);
H_u = H(used_total,:,:);

end