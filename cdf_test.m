
function cdf_test(Fbb_all)
ITER = size(Fbb_all,3);
%a = zeros(size(Fbb_all,1), size(Fbb_all,2));
%range_fbb = zeros(ITER, size(Fbb_all,1));
range_fbb = zeros(ITER*size(Fbb_all,1),1);
count_t = 0;


count_1 = 0;
for x = 1:ITER
for c=1:size(Fbb_all, 1)
    for r = 1:size(Fbb_all,2)
        count_t = count_t+1;
      %  if abs(Fbb_all(r,c,x))<10
   %   if r~=c
    range_fbb(count_t) = abs(Fbb_all(c,r,x));
        count_1 = count_1+1;
 %     end
       % end
    end
end
end

step = 0.0001;
range = 0:step:10;
pap_fbb = range_fbb; 
%for to = 1:size(Fbb_all,1)
h=hist(pap_fbb,range);
fbb_cdf = cumsum(h)/(sum(h));
hold on
 plot(range, fbb_cdf, 'LineWidth',1.5);
 set(gca, 'box', 'on')
end
%end
