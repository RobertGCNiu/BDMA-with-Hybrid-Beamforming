function [Wrf, Frf, H_u]= greedySelection(a_TX, a_RX,k_cluster,H)
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