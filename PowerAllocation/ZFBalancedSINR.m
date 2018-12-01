function gama = ZFBalancedSINR(G,F,rho)
Num_users = size(G,1);
t = zeros( Num_users,1);

syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 gama
vals = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16];

eqns = [];
for u =1: length(vals)
    for t_else  =1:Num_users
        t(t_else) = abs(G(:,u)'*F(:,t_else))^2;
    end
    eqn = calculateT(u, vals, t,gama, rho);
    eqns = [eqns, eqn];
end

[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16] = solve(eqns, vals);
p_tot = sum([x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16]);
[gama] = solve(p_tot == 3, gama);

end
