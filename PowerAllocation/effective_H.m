function [G] = effective_H(H,W,F)
%EFFECTIVE_H Summary of this function goes here
%   Detailed explanation goes here
N=size(H,1);
G = zeros(N, N);
for u=1:N
    C(:,:)=H(u,:,:);
    G(u,:)=W(:,u)'*C*F ;    % Effective channels
end
end

