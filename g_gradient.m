function [g_gradient] = g_gradient(p,H, M)

sigma = 0.0001;
w = ones(M,1);
g_gradient = 0;
e_i = zeros(M,1);


for i = 1 : M
    H_d = 0;
    for j = 1: M
         if j == i
            e_i(j) = 0;
        else
            e_i(j) = w(i)*H(j,i)/log(2);
        end
        if j ~= i
            H_d = H_d + H(j,i)*p(j);
        end
    end
    g_gradient = g_gradient + e_i/(sigma+H_d);
end

end
