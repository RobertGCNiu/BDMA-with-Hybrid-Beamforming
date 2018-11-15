function [Wrf, Frf, H_u]= select_interCluster(a_TX, a_RX,k_cluster,H,SNR)
% Num_users = size(a_TX,2);
% group_size = Num_users/k_cluster;
% user = zeros(k_cluster, 1);
% 
% for u = 1:Num_users
%     channel(:,:) = H(u,:,:);
%     alpha_u(u) = a_RX(:,u)'*channel*a_TX(:,u);
% end
% [~,alpha_sort] = sort(alpha_u,'d');
% selected_users = [];
% for k = 1:k_cluster
% user(k,1) = alpha_sort(k);
% selected_users = [selected_users, user(k,1)];
% end
% 
% temp_users = user;
% num_selected_users=length(selected_users);
% cluster_set = 1:k_cluster;
% 
% for x = 1:Num_users
%     if ~ismember(x,selected_users)
%     R = zeros(k_cluster,1);
%     for k = cluster_set
%         tem_k = [user(k,:), x];
%         temp_users(k,1:sum(temp_users(k,:)~=0)+1)  = tem_k;
%         for k_r = cluster_set
%             int_k = cluster_set;
%             int_k(find(int_k==k_r)) = [];     
%             for u_r = 1:sum(temp_users(k_r,:)~=0)
%             interference = 0;
%             if k_r ==k
%                 for k_else = 1:length(int_k)
%                     interference = interference + sum(temp_users(int_k(k_else),:)~=0)*sum(abs(a_RX(:,temp_users(k_r,u_r))'*channel*a_TX(:,temp_users(int_k(k_else),end-1))).^2);
%                 end
%                 R(k) = R(k) + log2(1+SNR*sum(abs(a_RX(:,temp_users(k_r,u_r))'*channel*a_TX(:, temp_users(k_r,:))).^2)/(SNR*interference+1));
%             else
%                 for k_else = 1:length(int_k)
%                     interference = interference + sum(temp_users(int_k(k_else),:)~=0)*sum(abs(a_RX(:,temp_users(k_r,u_r))'*channel*a_TX(:,temp_users(int_k(k_else),end-1))).^2);
%                 end
%                 R(k) = R(k) + log2(1+SNR*sum(abs(a_RX(:,temp_users(k_r,u_r))'*channel*a_TX(:, temp_users(k_r,end-1))).^2)/(SNR*interference+1));
%             end
%             end
%         end
%     end
%     [~,max_k_index]=max(R);
%     user(max_k_index,end+1) = x;
%     if length(user(max_k_index,:))>group_size
%         cluster_set(max_k_index) =[];
%     end
%     end
% end    
% 
% for k = 1:k_cluster
%     Frf(:,k*group_size+1:(k+1)*group_size) = a_TX(:, user(k,:));
%     Wrf(:,k*group_size+1:(k+1)*group_size) = a_RX(:, user(k,:));
%     H_u(k*group_size+1:(k+1)*group_size,:,:) = H(user(k,:),:,:);
% end













if nargin<5
    SNR = 1;
end

Num_users = size(a_TX,2);
group_size = Num_users/k_cluster;
user = zeros(k_cluster, group_size);

for u = 1:Num_users
    channel(:,:) = H(u,:,:);
    alpha_u(u) = a_RX(:,u)'*channel*a_TX(:,u);
end
[~,max_u1] = max(alpha_u);
selected_users = max_u1;
num_selected_users=length(selected_users);

k=1;
user_in_cluster = 1;
user(k,user_in_cluster) = max_u1;
while(num_selected_users<Num_users)
    if sum(user(k,:)~=0)==k_cluster
        k=k+1;
        user_in_cluster = 1;
    end
    R = zeros(Num_users,1);
    for x = 1:Num_users
        if ~ismember(x,selected_users)
            tem_x = [selected_users, x];
            for t =1:length(tem_x)
                if length(tem_x)<=group_size
                    eta=0;
                else 
                    eta=1;
                end
                channel(:,:) = H(tem_x(t),:,:);
                int_set = tem_x;
                    if floor(t/group_size)==floor(length(tem_x)/group_size)
                       int_set(floor(t/group_size)*group_size+1 : end) = [];
                       m_k = length(int_set);
                       v_k = a_TX(:,floor(t/group_size)*group_size+1 : end);
                    else
                       int_set(floor(t/group_size)*group_size+1 : (floor(t/group_size)+1)*group_size) = [];
                       m_k = length(int_set);
                       v_k = a_TX(:, floor(t/group_size)*group_size+1 : (floor(t/group_size)+1)*group_size);
                    end
                R(x) =R(x) + log2(1+SNR*sum(abs(a_RX(:,tem_x(t))'*channel*v_k).^2)/(SNR*m_k*sum(abs(a_RX(:,tem_x(t))'*channel*a_TX(:,int_set)).^2)*eta+1));
            end
        end
    end
    [~,selected_u]=max(R);
    selected_users = [selected_users,selected_u];
    num_selected_users = length(selected_users);
    %user_in_cluster = user_in_cluster+1;
    user(k,user_in_cluster) = selected_u;
end

for u =1:Num_users
    Frf(:,u) = a_TX(:,selected_users(u));
    Wrf(:,u) = a_RX(:,selected_users(u));
    H_u(u,:,:) = H(selected_users(u),:,:);
end



end
