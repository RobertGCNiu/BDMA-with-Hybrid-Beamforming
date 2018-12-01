
% x = 0:0.1:10;
% 
% r = 1-1./(1+x).^(4-1);%sqrt(0.5*(x.^2));%每个分量的方差为0.5
% F_sc = (1-1./(1+x).^(4-1)).^2;
% plot(x,r,'*-');hold on
% plot(x,F_sc,'r-')
% step = 0.0001;
% range = 0:step:10;
% h = hist(r,range);
% 
% pr_approx_cdf = cumsum(h)/(sum(h));
% plot(range,pr_approx_cdf,'r-','LineWidth',2)
apha = reshape(abs(appha),[1,128]);
hist(abs(apha))