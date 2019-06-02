function [right, wrong_i, wrong_index] = count_right_times(D1)
right = 0;
k = 5;
wrong_i = [];
wrong_index = [];
for i = 1: size(D1,3)
    x = D1(:,:,i);
    x_v = x(:);
    d_sort = sort(x_v);
    k_min  = d_sort(k);
    for c_check = 1:size(D1,2)
        x_c = x(:,c_check);
     count(c_check) = length(x_c(x_c<=k_min));
    end
   [~,index_count] = max(count);
if index_count == i 
    right= right+1;
else
    wrong_i = [wrong_i, i];
    wrong_index = [wrong_index, index_count];
end
end
