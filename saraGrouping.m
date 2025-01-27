
clear; close all

cd /Users/asf2200/Library/CloudStorage/Dropbox/projectEEGsca
tblData = readtable('saraScores.xlsx');
data{1} = table2array(tblData(:,2:end));
data{2}(1,:) = prctile(data{1},33,1);
data{2}(2,:) = prctile(data{1},66,1);

for i = 1:size(data{1},2)
    count = [1 1 1];
    for m = 1:size(data{1},1)
        if data{1}(m,i) < data{2}(1,i)
            grps{1,i}{count(1)} = tblData{m,1};
            count(1) = count(1) + 1;
        elseif data{1}(m,i) <= data{2}(2,i) & data{1}(m,i) >= data{2}(1,i)
            grps{2,i}{count(2)} = tblData{m,1};
            count(2) = count(2) + 1;
        else
            grps{3,i}{count(3)} = tblData{m,1};
            count(3) = count(3) + 1;
        end

    end
end