function p =  dcpower(H, M, p_all,sigma )

w = ones(M,1);
p_max = ones(M,1);
ri = 0;

p_initial = p_max;
t_old = 0;
max_value_dif = 100;
iter = 0;
while(abs(max_value_dif)>0.01)
    
iter = iter + 1;
cvx_begin quiet
variables p(M)
 
 f_p = 0;
 g_p = 0;
 for i = 1 : M
     H_f = 0;
     H_g = 0;
     for j = 1:M
         H_f = H_f +  H(j,i) * p(j);
        if j ~= i
            H_g = H_g + H(j,i) * p_initial(j);
        end
     end
     f_p  = f_p + w(i) * log(sigma + H_f)/log(2);   % f(p)
     g_p  = g_p + w(i) * log(sigma + H_g)/log(2);  % g(p)
 end
    gra_g = g_gradient(p_initial,H, M);   % gradient of g
    t= f_p - g_p - gra_g'*(p-p_initial);
    maximize(t)
    subject to
            for i = 1:M
                He = 0;
                for j = 1:M
                   if j ~= i
                        He = He + H(j,i) * p(j);
                   end
                end
                 H(i,i) * p(i) + (1-2^ri)*(He + sigma) >= 0;
                sum(p)<=p_all;
                p(i)>=0;
            end
cvx_end
    p_initial = p;
    max_value_dif = t-t_old;
    t_old = t;
end