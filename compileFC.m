%   Compile functional connectivity data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear

% Get relevant folders
[dirList,numAreas] = getFolders();

% Read in data from excel files
for i = 1:numAreas
    cd(dirList(1).folder)
    cd(dirList(i).name);

    [groupList,numGroups] = getFolders();

    for k = 1:numGroups
        cd(groupList(1).folder)
        cd(groupList(k).name);

        [freqList,numFreqBands] = getFolders();
        numFreqs(i,k) = numFreqBands;

        for m = 1:numFreqBands
            cd(freqList(1).folder)
            cd(freqList(m).name);

            % Read in data from excel sheets
            excelList = dir('*.xlsx');
            if length(excelList) == 0
                excelList = dir('*.csv');
                condExcept = 1;
            end

            numSxs = 0;
            for ii = 1:length(excelList)
                if startsWith(excelList(ii).name,'.') || startsWith(excelList(ii).name,'~')
                    continue
                else
                    numSxs = numSxs + 1;
                end
            end
            numFiles{i,k}(m) = numSxs;
            excelList = excelList(1:numSxs);

            for t = 1:length(excelList)
                if exist('condExcept')
                    data{i,k}{m,t} = csvread(excelList(t).name,1,1);
                else
                    data{i,k}{m,t} = xlsread(excelList(t).name);
                end
            end
        end
    end
end

cd(dirList(1).folder)

%% Extract relevant data

for i = 1:numAreas
    for k = 1:numGroups
        if i == 1
            rowIdxs = 1;
            colIdxs = 2:3;
        elseif i == 2
            rowIdxs = 8;
            rowIdxs2 = 9;
            colIdxs = 6:7;
        elseif i == 3
            rowIdxs = 4;
            colIdxs = 8:9;
        elseif i == 4
            rowIdxs = 5;
            colIdxs = 1:4;
        elseif i == 5 || i == 6 || i == 7
            rowIdxs = 1;
            colIdxs = 1:2;
        end

        for m = 1:numFreqs(i,k)
            for t = 1:numFiles{i,k}(m)
                dataTemp = data{i,k}{m,t}(rowIdxs,colIdxs);
                dataComp{i,k}{1,m}(t,:) = dataTemp;
                if i == 2
                    dataTemp = data{i,k}{m,t}(rowIdxs2,colIdxs);
                    dataComp{i,k}{2,m}(t,:) = dataTemp;
                end
            end

            dataAvgTemp = mean(dataComp{i,k}{1,m},1);
            dataCompAvg{i,k}{1,m} = dataAvgTemp;

            dataSEtemp = (std(dataComp{i,k}{1,m},0,1)) / (sqrt(numFiles{i,k}(m)));
            dataCompSE{i,k}{1,m} = dataSEtemp;

            if i == 2
                dataAvgTemp = mean(dataComp{i,k}{2,m},1);
                dataCompAvg{i,k}{2,m} = dataAvgTemp;
    
                dataSEtemp = (std(dataComp{i,k}{2,m},0,1)) / (sqrt(numFiles{i,k}(m)));
                dataCompSE{i,k}{2,m} = dataSEtemp;
            end
        end

    for m = 1:numFreqs(i,k)
        if i == 1
            figArrTemp = dataCompAvg{i,k}{m}(1,2);
            figArr{m}(k) = figArrTemp;

            figArrSEtemp = dataCompSE{i,k}{m}(1,2);
            figArrSE{m}(k) = figArrSEtemp;

            figArrTemp = dataComp{i,k}{m}(:,2);
            figArrSxs{m}(:,k) = figArrTemp;

            clear figArrTemp figArrSEtemp
        end

        if i == 2
            for j = 1:2 %number of brain regions (Crus1 and Crus 2)
                figArrTemp = dataCompAvg{i,k}{1,m}(1,j);
                figArrCrusPreC{j,m}(k) = figArrTemp;

                figArrTemp = dataCompAvg{i,k}{2,m}(1,j);
                figArrCrusSMA{j,m}(k) = figArrTemp;

                figArrSEtemp = dataCompSE{i,k}{1,m}(1,j);
                figArrSEcrusPreC{j,m}(k) = figArrSEtemp;

                figArrSEtemp = dataCompSE{i,k}{2,m}(1,j);
                figArrSEcrusSMA{j,m}(k) = figArrSEtemp;

                figArrSxsTemp = dataComp{i,k}{1,m}(:,j);
                figArrSxsCrusPreC{j,m}(1:numFiles{i,k}(m),k) = figArrSxsTemp;

                figArrSxsTemp = dataComp{i,k}{2,m}(:,j);
                figArrSxsCrusSMA{j,m}(1:numFiles{i,k}(m),k) = figArrSxsTemp;

                clear figArrTemp figArrSEtemp
            end
        end

        if i == 3
            figArrTemp = dataCompAvg{i,k}{m}(1,1);
            figArrLob8preC{m}(k) = figArrTemp;

            figArrTemp = dataCompAvg{i,k}{m}(1,2);
            figArrLob8sma{m}(k) = figArrTemp;

            figArrSEtemp = dataCompSE{i,k}{m}(1,1);
            figArrSElob8preC{m}(k) = figArrSEtemp;

            figArrSEtemp = dataCompSE{i,k}{m}(1,2);
            figArrSElob8sma{m}(k) = figArrSEtemp;

            figArrLob8Temp = dataComp{i,k}{m}(:,1);
            figArrSxsLob8preC{m}(:,k) = figArrLob8Temp;

            figArrLob8Temp = dataComp{i,k}{m}(:,2);
            figArrSxsLob8sma{m}(:,k) = figArrLob8Temp;

            clear figArrTemp figArrSEtemp
        end

        if i == 4
            figArrTemp = dataCompAvg{i,k}{m}(1,1);
            figArrLob45preC{m}(k) = figArrTemp;

            figArrTemp = dataCompAvg{i,k}{m}(1,2);
            figArrLob6sma{m}(k) = figArrTemp;

            figArrTemp = dataCompAvg{i,k}{m}(1,3);
            figArrCrus1preC{m}(k) = figArrTemp;

            figArrTemp = dataCompAvg{i,k}{m}(1,4);
            figArrCrus2preC{m}(k) = figArrTemp;

            figArrSEtemp = dataCompSE{i,k}{m}(1,1);
            figArrSElob45preC{m}(k) = figArrSEtemp;

            figArrSEtemp = dataCompSE{i,k}{m}(1,2);
            figArrSElob6sma{m}(k) = figArrSEtemp;

            figArrSEtemp = dataCompSE{i,k}{m}(1,3);
            figArrSEcrus1preC{m}(k) = figArrSEtemp;

            figArrSEtemp = dataCompSE{i,k}{m}(1,4);
            figArrSEcrus2preC{m}(k) = figArrSEtemp;

            figArrLob45temp = dataComp{i,k}{m}(:,1);
            figArrSxsLob45preC{m}(:,k) = figArrLob45temp;

            figArrLob6temp = dataComp{i,k}{m}(:,2);
            figArrSxsLob6sma{m}(:,k) = figArrLob6temp;

            figArrCrus1temp = dataComp{i,k}{m}(:,3);
            figArrSxsCrus1preC{m}(:,k) = figArrCrus1temp;

            figArrCrus2temp = dataComp{i,k}{m}(:,4);
            figArrSxsCrus2preC{m}(:,k) = figArrCrus2temp;

            clear figArrTemp figArrSEtemp
        end

        if i == 5
            figArrTemp = dataCompAvg{i,k}{1,m}(1,1);
            figArrCrusPreCnew{1,m}(k) = figArrTemp;

            figArrTemp = dataCompAvg{i,k}{1,m}(1,2);
            figArrCrusSMAnew{1,m}(k) = figArrTemp;

            figArrSEtemp = dataCompSE{i,k}{1,m}(1,1);
            figArrSEcrusPreCnew{1,m}(k) = figArrSEtemp;

            figArrSEtemp = dataCompSE{i,k}{1,m}(1,2);
            figArrSEcrusSMAnew{1,m}(k) = figArrSEtemp;

            figArrSxsTemp = dataComp{i,k}{1,m}(:,1);
            figArrSxsCrusPreCnew{1,m}(1:numFiles{i,k}(m),k) = figArrSxsTemp;

            figArrSxsTemp = dataComp{i,k}{1,m}(:,2);
            figArrSxsCrusSMAnew{1,m}(1:numFiles{i,k}(m),k) = figArrSxsTemp;

            clear figArrTemp figArrSEtemp
        elseif i == 6
            figArrTemp = dataCompAvg{i,k}{1,m}(1,1);
            figArrCrus2preCnew{1,m}(k) = figArrTemp;

            figArrTemp = dataCompAvg{i,k}{1,m}(1,2);
            figArrCrus2smaNew{1,m}(k) = figArrTemp;

            figArrSEtemp = dataCompSE{i,k}{1,m}(1,1);
            figArrSEcrus2preCnew{1,m}(k) = figArrSEtemp;

            figArrSEtemp = dataCompSE{i,k}{1,m}(1,2);
            figArrSEcrus2smaNew{1,m}(k) = figArrSEtemp;

            figArrSxsTemp = dataComp{i,k}{1,m}(:,1);
            figArrSxsCrus2preCnew{1,m}(1:numFiles{i,k}(m),k) = figArrSxsTemp;

            figArrSxsTemp = dataComp{i,k}{1,m}(:,2);
            figArrSxsCrus2smaNew{1,m}(1:numFiles{i,k}(m),k) = figArrSxsTemp;

            clear figArrTemp figArrSEtemp
        elseif i == 7
            figArrTemp = dataCompAvg{i,k}{1,m}(1,1);
            figArrLob8preCnew{1,m}(k) = figArrTemp;

            figArrTemp = dataCompAvg{i,k}{1,m}(1,2);
            figArrLob8smaNew{1,m}(k) = figArrTemp;

            figArrSEtemp = dataCompSE{i,k}{1,m}(1,1);
            figArrSElob8preCnew{1,m}(k) = figArrSEtemp;

            figArrSEtemp = dataCompSE{i,k}{1,m}(1,2);
            figArrSElob8smaNew{1,m}(k) = figArrSEtemp;

            figArrSxsTemp = dataComp{i,k}{1,m}(:,1);
            figArrSxsLob8preCnew{1,m}(1:numFiles{i,k}(m),k) = figArrSxsTemp;

            figArrSxsTemp = dataComp{i,k}{1,m}(:,2);
            figArrSxsLob8smaNew{1,m}(1:numFiles{i,k}(m),k) = figArrSxsTemp;

            clear figArrTemp figArrSEtemp
        end
    end

    end
end

%% Plot difference of function connectivities

new = [figArr{1}; figArr{2}];
plotMultiFreqs(figArrSxs,figArr,figArrSE,new,1,1,'Crus 1 x Precentral')

crus1preC = [figArrCrusPreC{1,1}; figArrCrusPreC{1,2}; figArrCrusPreC{1,3}; figArrCrusPreC{1,4}];
crus2preC = [figArrCrusPreC{2,1}; figArrCrusPreC{2,2}; figArrCrusPreC{2,3}; figArrCrusPreC{2,4}];

plotMultiFreqs(figArrSxsCrusPreC,figArrCrusPreC,figArrSEcrusPreC,crus1preC,2,1,'Crus 1 x Precentral')
plotMultiFreqs(figArrSxsCrusPreC,figArrCrusPreC,figArrSEcrusPreC,crus2preC,2,2,'Crus 2 x Precentral')

crus1sma = [figArrCrusSMA{1,1}; figArrCrusSMA{1,2}; figArrCrusSMA{1,3}; figArrCrusSMA{1,4}];
crus2sma = [figArrCrusSMA{2,1}; figArrCrusSMA{2,2}; figArrCrusSMA{2,3}; figArrCrusSMA{2,4}];

plotMultiFreqs(figArrSxsCrusSMA,figArrCrusSMA,figArrSEcrusSMA,crus1sma,2,1,'Crus 1 x SMA')
plotMultiFreqs(figArrSxsCrusSMA,figArrCrusSMA,figArrSEcrusSMA,crus2sma,2,2,'Crus 2 x SMA')

lob8preC = [figArrLob8preC{1,1}; figArrLob8preC{1,2}; figArrLob8preC{1,3}; figArrLob8preC{1,4}];
plotMultiFreqs(figArrSxsLob8preC,figArrLob8preC,figArrSElob8preC,lob8preC,3,1,'Lob 8 x Precentral')

lob8sma = [figArrLob8sma{1,1}; figArrLob8sma{1,2}; figArrLob8sma{1,3}; figArrLob8sma{1,4}];
plotMultiFreqs(figArrSxsLob8sma,figArrLob8sma,figArrSElob8sma,lob8sma,3,1,'Lob 8 x SMA')

crus1preCnew = [figArrCrusPreCnew{1,1}; figArrCrusPreCnew{1,2}; figArrCrusPreCnew{1,3}];
crus2preCnew = [figArrCrus2preCnew{1,1}; figArrCrus2preCnew{1,2}; figArrCrus2preCnew{1,3}];
lob8preCnew = [figArrLob8preCnew{1,1}; figArrLob8preCnew{1,2}; figArrLob8preCnew{1,3}];

plotMultiFreqs(figArrSxsCrusPreCnew{1},figArrCrusPreCnew{1},figArrSEcrusPreCnew{1},crus1preCnew(1,:),4,1,'Crus 1 x Precentral -- Alpha')
plotMultiFreqs(figArrSxsCrus2preCnew{1},figArrCrus2preCnew{1},figArrSEcrus2preCnew{1},crus2preCnew(1,:),4,1,'Crus 2 x Precentral -- Alpha')
plotMultiFreqs(figArrSxsLob8preCnew{1},figArrLob8preCnew{1},figArrSElob8preCnew{1},lob8preCnew(1,:),4,1,'Lob8 x Precentral -- Alpha')

plotMultiFreqs(figArrSxsCrusPreCnew{2},figArrCrusPreCnew{2},figArrSEcrusPreCnew{2},crus1preCnew(2,:),4,1,'Crus 1 x Precentral -- Beta')
plotMultiFreqs(figArrSxsCrus2preCnew{2},figArrCrus2preCnew{2},figArrSEcrus2preCnew{2},crus2preCnew(2,:),4,1,'Crus 2 x Precentral -- Beta')
plotMultiFreqs(figArrSxsLob8preCnew{2},figArrLob8preCnew{2},figArrSElob8preCnew{2},lob8preCnew(2,:),4,1,'Lob8 x Precentral -- Beta')

plotMultiFreqs(figArrSxsCrusPreCnew{3},figArrCrusPreCnew{3},figArrSEcrusPreCnew{3},crus1preCnew(3,:),4,1,'Crus 1 x Precentral -- Delta')
plotMultiFreqs(figArrSxsCrus2preCnew{3},figArrCrus2preCnew{3},figArrSEcrus2preCnew{3},crus2preCnew(3,:),4,1,'Crus 2 x Precentral -- Delta')
plotMultiFreqs(figArrSxsLob8preCnew{3},figArrLob8preCnew{3},figArrSElob8preCnew{3},lob8preCnew(3,:),4,1,'Lob8 x Precentral -- Delta')

crus1smaNew = [figArrCrusSMAnew{1,1}; figArrCrusSMAnew{1,2}; figArrCrusSMAnew{1,3}];
crus2smaNew = [figArrCrus2smaNew{1,1}; figArrCrus2smaNew{1,2}; figArrCrus2smaNew{1,3}];
lob8smaNew = [figArrLob8smaNew{1,1}; figArrLob8smaNew{1,2}; figArrLob8smaNew{1,3}];

plotMultiFreqs(figArrSxsCrusSMAnew{1},figArrCrusSMAnew{1},figArrSEcrusSMAnew{1},crus1smaNew(1,:),4,1,'Crus 1 x SMA -- Alpha')
plotMultiFreqs(figArrSxsCrus2smaNew{1},figArrCrus2smaNew{1},figArrSEcrus2smaNew{1},crus2smaNew(1,:),4,1,'Crus 2 x SMA -- Alpha')
plotMultiFreqs(figArrSxsLob8smaNew{1},figArrLob8smaNew{1},figArrSElob8smaNew{1},lob8smaNew(1,:),4,1,'Lob8 x SMA -- Alpha')

plotMultiFreqs(figArrSxsCrusSMAnew{2},figArrCrusSMAnew{2},figArrSEcrusSMAnew{2},crus1smaNew(2,:),4,1,'Crus 1 x SMA -- Beta')
plotMultiFreqs(figArrSxsCrus2smaNew{2},figArrCrus2smaNew{2},figArrSEcrus2smaNew{2},crus2smaNew(2,:),4,1,'Crus 2 x SMA -- Beta')
plotMultiFreqs(figArrSxsLob8smaNew{2},figArrLob8smaNew{2},figArrSElob8smaNew{2},lob8smaNew(2,:),4,1,'Lob8 x SMA -- Beta')

plotMultiFreqs(figArrSxsCrusSMAnew{3},figArrCrusSMAnew{3},figArrSEcrusSMAnew{3},crus1smaNew(3,:),4,1,'Crus 1 x SMA -- Delta')
plotMultiFreqs(figArrSxsCrus2smaNew{3},figArrCrus2smaNew{3},figArrSEcrus2smaNew{3},crus2smaNew(3,:),4,1,'Crus 2 x SMA -- Delta')
plotMultiFreqs(figArrSxsLob8smaNew{3},figArrLob8smaNew{3},figArrSElob8smaNew{3},lob8smaNew(3,:),4,1,'Lob8 x SMA -- Delta')

%% Permutation Student's t-test (unequal variances) - Brainstorm version
for i = 1:length(figArr)
end
for i = 1:length(figArrCrusPreC)
    statVecCrus1preCctrl(1,:) = figArrSxsCrusPreC{1,i}(:,1);
    statVecCrus1preCatx(1,:) = figArrSxsCrusPreC{1,i}(:,2);
    
    pValueC1preC(i) = unpairedTtest(statVecCrus1preCctrl,statVecCrus1preCatx);
end

for i = 1:length(figArrCrusPreC)
    statVecCrus2preCctrl(1,:) = figArrSxsCrusPreC{2,i}(:,1);
    statVecCrus2preCatx(1,:) = figArrSxsCrusPreC{2,i}(:,2);
    
    pValueC2preC(i) = unpairedTtest(statVecCrus2preCctrl,statVecCrus2preCatx);
end

for i = 1:length(figArrCrusSMA)
    statVecCrus1smaCtrl(1,:) = figArrSxsCrusSMA{1,i}(:,1);
    statVecCrus1smaAtx(1,:) = figArrSxsCrusSMA{1,i}(:,2);
    
    pValueC1sma(i) = unpairedTtest(statVecCrus1smaCtrl,statVecCrus1smaAtx);
end

for i = 1:length(figArrCrusSMA)
    statVecCrus2smaCtrl(1,:) = figArrSxsCrusSMA{2,i}(:,1);
    statVecCrus2smaAtx(1,:) = figArrSxsCrusSMA{2,i}(:,2);
    
    pValueC2sma(i) = unpairedTtest(statVecCrus2smaCtrl,statVecCrus2smaAtx);
end

for i = 1:length(figArrLob8preC)
    statVecLob8preCctrl(1,:) = figArrSxsLob8preC{1,i}(:,1);
    statVecLob8preCatx(1,:) = figArrSxsLob8preC{1,i}(:,2);
    
    pValueL8preC(i) = unpairedTtest(statVecLob8preCctrl,statVecLob8preCatx);
end

for i = 1:length(figArrLob8sma)
    statVecLob8smaCtrl(1,:) = figArrSxsLob8sma{1,i}(:,1);
    statVecLob8smaAtx(1,:) = figArrSxsLob8sma{1,i}(:,2);
    
    pValueL8sma(i) = unpairedTtest(statVecLob8smaCtrl,statVecLob8smaAtx);
end
% y = [ctrlCrus1alphaPop(1,2) atxCrus1alphaPop(1,2)];
% z = [ctrlCrus1betaPop(1,2) atxCrus1betaPop(1,2)];
% new = [y; z];
% 
% s = [ctrlCrus1alphaPopSE(1,2) atxCrus1alphaPopSE(1,2)];
% s1 = [ctrlCrus1betaPopSE(1,2) atxCrus1betaPopSE(1,2)];
% 
% dd =[ctrlCrus1alpha(:,2) atxCrus1alpha(:,2)];
% ee = [ctrlCrus1beta(:,2) atxCrus1beta(:,2)];
% 
% x = [0.86 1.14];
% x1 = [1.86 2.14];
% xx = repmat(0.86,size(dd(:,1)));
% xx(:,2) = repmat(1.14,size(dd(:,1)));
% x2 = repmat(1.86,size(dd(:,1)));
% x2(:,2) = repmat(2.14,size(dd(:,1)));
% 
% labels = {'Alpha','Beta'};
% figure(); hold on
% h = bar(new,'grouped');
% legend('Ctrl','Ataxia')
% set(gca,'XTickLabel',labels);
% 
% errorbar(x,y,s,'.','Color','k','LineWidth',1)
% errorbar(x1,z,s1,'.','Color','k','LineWidth',1)
% scatter(xx,dd,'k')
% scatter(x2,ee,'k')
% xticks([1 2])
% xticklabels({'Alpha', 'Beta'})
% legend('Ctrl','Ataxia')
% ylabel('Functional connectivity (Rest - action)')
% set(gca,'FontSize',16)
% 
% % h(2).FaceColor = 'r';
% % h(4).FaceColor = 'r';
% 
% %% Plot PR data
% figTitle = {'crus 1','L4/5','L6','Crus 1','Crus 2'};
% for ii = 1:length(ctrlPRalphaPop)
%     y = [ctrlPRalphaPop(1,ii) atxPRalphaPop(1,ii)];
%     z = [ctrlPRbetaPop(1,ii) atxPRbetaPop(1,ii)];
%     new = [y; z];
%     
%     s = [ctrlPRalphaPopSE(1,ii) atxPRalphaPopSE(1,ii)];
%     s1 = [ctrlPRbetaPopSE(1,ii) atxPRbetaPopSE(1,ii)];
%     
%     dd =[ctrlPRalpha(:,ii) atxPRalpha(:,ii)];
%     ee = [ctrlPRbeta(:,ii) atxPRbeta(:,ii)];
%     
%     x = [0.86 1.14];
%     x1 = [1.86 2.14];
%     xx = repmat(0.86,size(dd(:,1)));
%     xx(:,2) = repmat(1.14,size(dd(:,1)));
%     x2 = repmat(1.86,size(dd(:,1)));
%     x2(:,2) = repmat(2.14,size(dd(:,1)));
%     
%     labels = {'Alpha','Beta'};
%     figure(); hold on
%     h = bar(new,'grouped');
%     legend('Ctrl','Ataxia')
%     set(gca,'XTickLabel',labels);
%     
%     errorbar(x,y,s,'.','Color','k','LineWidth',1)
%     errorbar(x1,z,s1,'.','Color','k','LineWidth',1)
%     scatter(xx,dd,'k')
%     scatter(x2,ee,'k')
%     xticks([1 2])
%     xticklabels({'Alpha', 'Beta'})
%     legend('Ctrl','Ataxia')
%     ylabel('Functional connectivity (Rest - action)')
%     title(figTitle{ii+1})
%     set(gca,'FontSize',16)
%     
% %     h(2).FaceColor = 'r';
% %     h(4).FaceColor = 'r';
% end
