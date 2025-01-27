
data = {ctrls,ataxia,et};
for t = 1:3
    for i = 1:size(ctrls,2)
        aucTotal{t}(i) = trapz(data{t}(:,i));
        coiThresh = aucTotal{t}(i)/2;
        cumAuc{t}(i,:) = cumtrapz(data{t}(:,i));
        coi{t}(i) = find(cumAuc{t}(i,:)>coiThresh,1,'first') + 7;
    end
end

aucTotalAvg = cellfun(@(x) mean(x),aucTotal,'UniformOutput',false);
