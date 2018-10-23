function [F] = normalize_f(F,a)
%NORMALIZE_F Summary of this function goes here
%   Detailed explanation goes here
F=F/diag(sqrt(sum(abs(a*F).^2)));
% for u=1:size(F,1)
%     F(:,u)=F(:,u)/sqrt((a*F(:,u))'*(a*F(:,u)));
% end
end

