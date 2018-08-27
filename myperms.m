function [ outputlist ] = perm( list_org , chose )

    vecc = nchoosek(list_org, chose);
    k = size(vecc,1);
    for k_v = 1:k
        
         while ~isempty(list_org)
          list_org(find(ismember(vecc(k_v,:), list_org))) = [];
          vecc_2 = nchoosek(list_org,chose);
            for k_v2 = 1:size(vecc_2,1)
                new_vcc = [new_vcc vecc(k_v2)];
            end
         end
          outputlist = [outputlist; new_vcc];
        perm(list_org, chose);
    end


%         for ntimes = 1:length(list/chose)
% 
%         end

    
    
end

