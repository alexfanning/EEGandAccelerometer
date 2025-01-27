function [] = plotMultiFreqs(sxData,arrData,arrSE,grouped,type,row,titleName)

x = [0.86 1.14];
x1 = [1.86 2.14];
x2 = [2.86 3.14];
x3 = [3.86 4.14];
x4 = 1;
x5 = 2;

if type == 1
    xx = repmat(0.86,size(sxData{1}(:,1)));
    xx(:,2) = repmat(1.14,size(sxData{1}(:,1)));
    xx2 = repmat(1.86,size(sxData{1}(:,1)));
    xx2(:,2) = repmat(2.14,size(sxData{1}(:,1)));
elseif type == 2 || type == 3
    xx = repmat(0.86,size(sxData{1}(:,1)));
    xx(:,2) = repmat(1.14,size(sxData{1}(:,1)));
    xx2 = repmat(1.86,size(sxData{1}(:,1)));
    xx2(:,2) = repmat(2.14,size(sxData{1}(:,1)));
    xx3 = repmat(2.86,size(sxData{1}(:,1)));
    xx3(:,2) = repmat(3.14,size(sxData{1}(:,1)));
    xx4 = repmat(3.86,size(sxData{1}(:,1)));
    xx4(:,2) = repmat(4.14,size(sxData{1}(:,1)));
elseif type == 4
    xx = repmat(0.86,size(sxData(:,1)));
    xx(:,2) = repmat(1.14,size(sxData(:,1)));
    xx2 = repmat(1.86,size(sxData(:,1)));
    xx2(:,2) = repmat(2.14,size(sxData(:,1)));
    xx3 = repmat(1,size(sxData(:,1)));
    xx4 = repmat(2,size(sxData(:,2)));
end

if type == 1
    labels = {'Alpha', 'Beta'};
elseif type == 2 || type == 3
    labels = {'Alpha', 'Beta', 'Delta', 'Theta'};
elseif type == 4
    labels = {'Ctrl','Ataxia'};
end


if type == 4
    figure(); hold on
    h = bar(grouped,'FaceColor','flat');
    h.CData(1,:) = [0.9 0 0.1];
    set(gca,'XTickLabel',labels);
else
    figure(); hold on
    h = bar(grouped,'grouped');
    legend('Alpha','Beta')
    set(gca,'XTickLabel',labels);
end

if type == 1
    errorbar(x,arrData{1},arrSE{1},'.','Color','k','LineWidth',1)
    errorbar(x1,arrData{2},arrSE{2},'.','Color','k','LineWidth',1)
    scatter(xx,sxData{1},'k')
    scatter(xx2,sxData{2},'k')
elseif type == 2 || type == 3
    errorbar(x,arrData{row,1},arrSE{row,1},'.','Color','k','LineWidth',1)
    errorbar(x1,arrData{row,2},arrSE{row,2},'.','Color','k','LineWidth',1)
    scatter(xx,sxData{row,1},'k')
    scatter(xx2,sxData{row,2},'k')

    errorbar(x2,arrData{row,3},arrSE{row,3},'.','Color','k','LineWidth',1)
    errorbar(x3,arrData{row,4},arrSE{row,4},'.','Color','k','LineWidth',1)
    scatter(xx3,sxData{row,3},'k')
    scatter(xx4,sxData{row,4},'k')
elseif type == 4
    errorbar(x4,arrData(1),arrSE(1),'.','Color','k','LineWidth',1)
    errorbar(x5,arrData(2),arrSE(2),'.','Color','k','LineWidth',1)
    scatter(xx3,sxData(:,1),'k')
    scatter(xx4,sxData(:,2),'k')
end

if type == 1
    xticks([1 2])
    xticklabels({'Alpha', 'Beta'})
    legend('Ctrl','Ataxia')
    ylabel('Functional connectivity (Rest - action)')
elseif type == 2 || type == 3
    xticks([1 2 3 4])
    xticklabels({'Alpha', 'Beta','Delta','Theta'})
    legend('Ctrl','Ataxia')
    ylabel('Functional connectivity (Rest)')
elseif type == 4
    xticks([1 2])
    xticklabels({'Ctrl','Ataxia'})
    ylabel('Functional connectivity')
end
title(titleName)
set(gca,'FontSize',16)