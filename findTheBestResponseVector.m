function wrf_dft = findTheBestResponseVector(Wrf_org)
[m,n] = size(Wrf_org);
wrf_dft = zeros(m,n);
wrf_dft_collection = 1/sqrt(m)*dftmtx(m);
for t = 1:n
     rel = abs(wrf_dft_collection*Wrf_org(:,t));
         [~,t_opt] = min(rel);
    wrf_dft(:,t) = wrf_dft_collection(:,t_opt);
end

end