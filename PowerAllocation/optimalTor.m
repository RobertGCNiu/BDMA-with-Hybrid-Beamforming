function tor = optimalTor(q,heff,noise_k, P_max)
N = length(q);
M = size(heff);
below = 0;
for n =1:N
        int_set = 1:N;
        int_set(n) = [];
        below_inner = 0;
        for k = int_set
            below_inner = below_inner+q(k)/N*(heff(:,k)*heff(:,k)'+noise_k * eye(M));
        end
        below = below+(1/N*heff(:,n)'*(below_inner)^(-1)*heff(:,n))^(-1);
end
tor = N*P_max/below;