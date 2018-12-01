function [eqn] = calculateT(u, vals, T,gama, rho)

    int_set = 1:length(vals);
    int_set(u) = [];
    eqn = [gama == vals(u)*T(u)/(vals(int_set)*T(int_set)+rho)];

end

