function [p,flag_solve] =  dcpower(H, M, p_all,sigma )

flag_solve = 0;
w = ones(M,1);
p_max = ones(M,1)*p_all/M;
ri = -10;
flag = 0;
p_initial = p_max;
t_old = 0;
max_value_dif = 100;
iter = 0;
while(abs(max_value_dif)>0.01)
    
iter = iter + 1;
cvx_begin quiet
%cvx_precision best
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
    gra_g = g_gradient(p_initial,H, M, sigma);   % gradient of g
    t= f_p - g_p - gra_g'*(p-p_initial);
    maximize(t)
    subject to
          sum(p)<=p_all;
          for i = 1:M
            p(i)>=0.00001;
          end
  %        if flag ==1
%             for i = 1:M
%                 He = 0;
%                 for j = 1:M
%                    if j ~= i
%                         He = He + H(j,i) * p(j);
%                    end
%                 end
%                  H(i,i) * p(i) + (1-2^ri)*(He + sigma) >= 0; 
%             end
  %        end
cvx_end
    p_initial = p;
    max_value_dif = t-t_old;
    t_old = t;
%    if t_old >=30
        flag = 1;
  %  end
end

%if t_old <1000
    flag_solve = 1;
%end