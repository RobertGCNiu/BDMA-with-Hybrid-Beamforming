  clear;clc
forlderlist = dir('C:\Users\robert\Documents\GitHub\BDMA-with-Hybrid-Beamforming\test\*');


right = 0;
k = 5;
%forlder = 3:112;
right_times = 0;
range_folder = 3:112;
tic
for test_index_inter  = 4:4:40
    forlder_list = 0;
for forlder = range_folder
    forlder_list = forlder_list+1;
    name_ford = forlderlist(forlder).name;
    namelist = dir( [string(name_ford)+'\*.txt']);
    addpath(genpath(forlderlist(forlder).name ));
    
    %%%%%%%%%%%%%%%%%%%%%%
   
    len = length(namelist);
    test_times = 0;
    for test_index = test_index_inter
        test_times = test_times+1;
        test_x1 = load(namelist(test_index).name);
        
        test_x = sqrt(test_x1(:,1).^2 + test_x1(:,2).^2+test_x1(:,3).^2);
        %%%%%%%%%%%%%%%%%%%%%%
        
        for_different_forlder = 0;
        for forlder_test =range_folder
            for_different_forlder = for_different_forlder+1;
            name_ford = forlderlist(forlder_test).name;
            namelist = dir( [string(name_ford)+'\*.txt']);
            addpath(genpath(forlderlist(forlder_test).name ));
            time_list = 0;
            for i = 4:4:40
                if test_index ~= i
                time_list = time_list+1;
                file_name1{i}=namelist(i).name;
                x1= load(file_name1{i});
                data1= sqrt(x1(:,1).^2 + x1(:,2).^2+x1(:,3).^2);     
                D1(time_list, for_different_forlder, forlder_list) = dtw(test_x, data1);
                end
            end
        end
    end
end  
[right_count, wrong_count,wrong_index] = count_right_times(D1);
right_times = right_times + right_count;
if ~isempty(wrong_count)
    wrong_all = [test_index_inter wrong_count]
    wrong_index_all =[ test_index_inter wrong_index]
end
end
toc
%sort(D1)

% dist11=dtw(test_x1,data1);
% dist21 = dtw(test_x2,data1);
% dist22 = dtw(test_x2,data2);
% dist12 = dtw(test_x1,data2);