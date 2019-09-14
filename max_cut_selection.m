function max_cut_selection(a_TX, a_RX,H)

N_RX = size(a_RX,1);

N_TX = size(a_TX,1);
Num_users = size(a_TX,2);
fid = fopen(['test.txt'], 'w');
total_u = Num_users/2*(Num_users-1);
fprintf(fid, '%d %d\n', [Num_users total_u]);
for k_1 = 1: Num_users
    H_user = zeros(N_RX, N_TX);
    H_user(:,:) = H(k_1,:,:);
    for k_2 = 1:Num_users/2
        if k_2 ~= k_1
        weight_t = 1000*log2(abs(a_RX(:,k_1)'*H_user*a_TX(:,k_2)));
         fprintf(fid, '%d %d %d\n', [k_1  k_2  -round(weight_t)]);
        end
    end
end

fclose(fid);
end