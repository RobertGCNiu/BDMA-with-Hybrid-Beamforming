function [q_k,tor] = funcofq(q,u,heff,noise_k,P_max)
N = length(q);
tor = optimalTor(q,heff,noise_k, P_max);
    int_set = 1:N;
    int_set(u) = [];
    inner_sum = 0;
    for k = int_set
       inner_sum = inner_sum + q(k)/N*heff(:,k)*heff(:,k)' +noise_k*eye(N);
    end
    below = 1/N*heff(:,u)' * inner_sum^(-1) * heff(:,u);
    q_k = tor/below;