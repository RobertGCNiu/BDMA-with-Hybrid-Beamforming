function [weight_all] = max_cut_selection(a_TX, a_RX,H)

N_RX = size(a_RX,1);

N_TX = size(a_TX,1);
Num_users = size(a_TX,2);
weight_all = [];
for k_1 = 1: Num_users
    H_user = zeros(N_RX, N_TX);
    H_user(:,:) = H(k_1,:,:);
    for k_2 = 1:Num_users
        if k_2 ~= k_1
        weight_t = log2(abs(a_RX(:,k_1)'*H_user*a_TX(:,k_2)));
        weight_all = [weight_all;k_1, k_2, -round(weight_t)];
        end
    end
end


end