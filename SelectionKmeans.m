function [Wrf_cl, Frf_cl, H_cl] = SelectionKmeans(a_TX_select, a_RX_select,K,H,Num_RF)

Num_users = size(a_TX_select,2);
%% find the first beam with largest channel gain
alpha_u = zeros(Num_users,1);
cluster_index = zeros(Num_users,1);

for u = 1:Num_users
    channel(:,:) = H(u,:,:);
    alpha_u(u) = a_RX_select(:,u)'*channel*a_TX_select(:,u);
end
[~, index_1] = max(alpha_u);
Wrf_cl(:,1) = a_RX_select(:,index_1);
Frf_cl(:,1) = a_TX_select(:,index_1);
H_cl(1,:,:) = H(index_1,:,:);
cluster_index(index_1) = 1;
index_all = index_1;
for k = 2:K
    alpha_u = zeros(Num_users,1);
    for u = 1:Num_users
        if ~ismember(u, index_all)
        channel(:,:) = H(u,:,:);
        alpha_u(u) =sum(abs(a_TX_select(:,u)'*a_TX_select(:,index_all)).^2);
        else
            alpha_u(u) = inf;
        end
    end
    [~, index_k] = min(alpha_u);
    cluster_index(index_k) = k;
    index_all = [index_all, index_k];
end

for u = 1:Num_users
    if ~ismember(u,index_all)
        alpha_k = zeros(K,1);
        for k = 1:K
            if sum(cluster_index==k)<Num_RF
            alpha_k(k) = sum(abs(a_TX_select(:,u)'*a_TX_select(:, find(cluster_index==k))).^2);
            end
        end
        [~,index_k] = max(alpha_k);
        cluster_index(u) = index_k;
        index_all = [index_all, u];
    end
end

    k_all = 0;
    for k = 1:K 
    Frf_cl(:, k_all+1: k_all+sum(cluster_index==k)) = a_TX_select(:, find(cluster_index==k));
    Wrf_cl(:, k_all+1: k_all+sum(cluster_index==k)) = a_RX_select(:, find(cluster_index==k));
    H_cl(k_all+1: k_all+sum(cluster_index==k),:,:) = H(find(cluster_index==k),:, :);
    k_all = k_all+sum(cluster_index==k);
    end



end