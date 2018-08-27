function [Wrf, Frf, H_u]= clusterByBeam(a_TX, a_RX,k_cluster,H)
N_TX = size(a_TX,1);

Num_users = size(a_TX,2);
users_kk = Num_users/k_cluster;
all_comb = myperms(1:Num_users);
combnum = size(all_comb,1);
a_TX_kk = zeros(N_TX, users_kk,k_cluster);
sum_residual = zeros(combnum,1);
for comb = 1:combnum
    vecc = all_comb(comb,:);
    for kk = 1:k_cluster
        a_TX_kk(:, :, kk) = a_TX(:, vecc((users_kk*(kk-1)+1):(users_kk*kk)));
    end
    
    
    for kk = 1:k_cluster
        for k_index = 1:k_cluster
            if k_index ~= kk
                sum_residual(comb) = sum_residual(comb) + norm(a_TX_kk(:,:,kk)'*a_TX_kk(:,:,k_index));
            end
        end
    end
    
end
       
[~,min_comb] = min(sum_residual);
vecc_min = all_comb(min_comb,:);

for kk = 1:k_cluster
    Frf(:,(users_kk*(kk-1)+1):users_kk*kk) = a_TX(:, vecc_min((users_kk*(kk-1)+1):users_kk*kk));
    Wrf(:,(users_kk*(kk-1)+1):users_kk*kk) = a_RX(:, vecc_min((users_kk*(kk-1)+1):users_kk*kk));
    H_u = H(vecc_min,:,:);
end

end

