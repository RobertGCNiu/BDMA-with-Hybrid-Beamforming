t = 10+rand(10,1);
rho = 1;
p = rand(10,1);
p_old = p+1;

syms x1 x2 x3 gamma
[eqn] = calculateT(x1,t,gamma, rho);
x = [x1,x2,x3];
% gamma = x1*t(1)/(x2*t(2)+rho);
% gamma = x2*t(2)/(x1*t(1)+rho);
eqns = [];
for u =1: length(x)
    eqn = calculateT(x(u),t,gamma, rho);
    eqns = [eqns, eqn];
end
% eqns = [gamma == x1*t(1)/(x3*t(3)+x2*t(2)+rho), gamma == x2*t(2)/(x1*t(1)+x3*t(3)+rho), gamma == x3*t(3)/(x1*t(1)+x2*t(2)+rho)];
% vals = [x1, x2,x3];
 %[x1,x2,x3] = solve(eqns, vals);
% p_tot = x1+x2+x3;
%[gamma] = solve(p_tot == 3, gamma)



% while norm(p_old - p)>0.001
%     p_old = p;
% for u = 1:10
%     t_else = 1:10;
%     t_else(u) = [];
%     p(u) = (p(t_else)'*t(t_else)+rho)/t(u);
% end
% norm(p_old-p)
% end