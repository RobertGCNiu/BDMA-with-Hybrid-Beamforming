clear;clc
forlderlist = dir('C:\Users\robert\Documents\GitHub\BDMA-with-Hybrid-Beamforming\test\*');

forlder_list = 0;
%forlder = 3:112;
for forlder = 3:10
    forlder_list = forlder_list+1;
    name_ford = forlderlist(forlder).name;
    namelist = dir( [string(name_ford)+'\*.txt']);
    addpath(genpath(forlderlist(forlder).name ));
    
    % namelist1 = dir('C:\Users\robert\Documents\GitHub\BDMA-with-Hybrid-Beamforming\test\1e\*.txt');
    % namelist2 = dir('C:\Users\robert\Documents\GitHub\BDMA-with-Hybrid-Beamforming\test\1w\*.txt');
    %%%%%%%%%%%%%%%%%%%%%%
    len = length(namelist);
    test_x1 = load(namelist(len).name);
   
    test_x = sqrt(test_x1(:,1).^2 + test_x1(:,2).^2+test_x1(:,3).^2);
    %%%%%%%%%%%%%%%%%%%%%%
  
    for_different_forlder = 0;
    for forlder_test = 3:10
        for_different_forlder = for_different_forlder+1;
        name_ford = forlderlist(forlder_test).name;
        namelist = dir( [string(name_ford)+'\*.txt']);
        addpath(genpath(forlderlist(forlder_test).name ));
        time_list = 0;
        for i = 4:4:len-1
            time_list = time_list+1;
            file_name1{i}=namelist(i).name;
            x1= load(file_name1{i});
            data1= sqrt(x1(:,1).^2 + x1(:,2).^2+x1(:,3).^2);
            D1(time_list, for_different_forlder, forlder_list) = dtw(test_x, data1);
%           if forlder_test ==3
%            D1 = dtw(test_x, data1)
%           else
%            D2   = dtw(test_x, data1)
%           end
        end
    end
    
end
right = 0;
k = 5;
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
end
end
%sort(D1)

% dist11=dtw(test_x1,data1);
% dist21 = dtw(test_x2,data1);
% dist22 = dtw(test_x2,data2);
% dist12 = dtw(test_x1,data2);