function [maxsir,p] = GMtraixGeneration(T)
N_users = size(T,1);
G = zeros(N_users, N_users);
for u =1:N_users
    G(u,:) = T(u, :)/T(u,u);
end

[power, lamda]=eig(G);
eigvalue = diag(lamda);
sir = 1./(eigvalue-1);
sir(imag(sir)~=0) = -10000;
    [maxsir,maxindex] = max(real(sir));
    p = power(:,maxindex);
end