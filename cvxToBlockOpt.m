function x = cvxToBlockOpt(He_fzf, Fbb_fzf, k_cluster)

[m,n] = size(Fbb_fzf);
Num_users = n;
Num_group = Num_users/k_cluster;
cvx_begin quiet
        variable x(m,n) complex
        minimize(norm(He_fzf*x - He_fzf*Fbb_fzf,'fro'))
        subject to 
           for k = 0:k_cluster-1
                for col = (1+Num_group*(k+1)) : Num_users
                    for row =  (1+Num_group*k) : (Num_group+Num_group*k)
                            x(col,row) == 0;
                    end
                end
                for row = (1+Num_group*(k+1)) : Num_users
                    for col =  (1+Num_group*k) : (Num_group+Num_group*k)
                            x(col,row) == 0;
                    end
                end
            end
cvx_end