function [R] = PGH(H,p,sigma)
R = 0;
SNR = p/sigma;
Num_users = length(p);
        for u=1:1:Num_users
            Int_set=[]; % interference index
            for i=1:1:Num_users
                if(i~=u)
                    Int_set=[Int_set i]; 
                end
            end
            
            
            SINR_BS=(SNR(u).*H(u,u))/(sum(SNR(Int_set)'.*H(u,Int_set))+1);
            R= R + log2(1+SINR_BS);
        end

end

