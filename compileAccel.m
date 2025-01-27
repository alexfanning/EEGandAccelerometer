%   Compiles data from multiple subjects
%
%   Alex Fanning, 10/19/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grab folder information
clear; close all;

[groupList,~] = getFolders(3);

%% Loop through each subjects file and extract data

params.numSxs = NaN(1);
mPsd = cell(1);
nPsdMean = cell(1);
nPsdMax = cell(1);
nPsdMaxFreq = cell(1);
pkAmp = cell(1);
rngeAmp = cell(1);
nPsdAuc = cell(1);
nPsdAucPct = cell(1);
freq = cell(1);
row = [7,8,9; 13,14,15; 2,3,4; 6,7,8; 10,11,12; 14,15,16; 18,19,20; 23,24,25];
row2 = [2,3,4,5];

for j = 1:length(groupList)
    cd(groupList(j).name)
    [dirList,~] = getFolders(4);
    
    % Load parameter info
    if j == 1
        load(uigetfile('*.mat'),'params')
    end
    
    % Set number of timepoints and epochs
    params.numSxs(j) = length(dirList);

    % Loop through each file and extract data from certain variables
    for k = 1:params.numSxs(j)
        params = rmfield(params,'taskNums');
        fileName = dirList(k).name;
        params.names = sheetnames(fileName);
        params = taskList(params,params,nPsdMax,2);
        for t = 1:numel(params.taskNums)
            if params.taskNums(t) ~= 6
                dataSx{j,t} = readmatrix(fileName,'Sheet',params.names(t));
                nPsdMax{j,params.taskNums(t)}(k,:) = dataSx{j,t}(1:3,4);
                nPsdMaxFreq{j,params.taskNums(t)}(k,:) = dataSx{j,t}(1:3,5);
                nPsdAuc{j,params.taskNums(t)}(k,:) = dataSx{j,t}(1:3,6);
                nPsdAucPct{j,params.taskNums(t)}(k,:) = dataSx{j,t}(1:3,7);
                if params.taskNums(t)~=5
                    freq{j,params.taskNums(t)}(k,:) = dataSx{j,t}(1:3,8);
                    pkAmp{j,params.taskNums(t)}(k,:) = dataSx{j,t}(1:3,9);
                    rngeAmp{j,params.taskNums(t)}(k,:) = dataSx{j,t}(1:3,10);
                    ipi{j,params.taskNums(t)}(k,:) = dataSx{j,t}(1:3,11);
                    cvIpi{j,params.taskNums(t)}(k,:) = dataSx{j,t}(1:3,12);
                    cvPkAmp{j,params.taskNums(t)}(k,:) = dataSx{j,t}(1:3,13);
                    cvSlWndwAmp{j,params.taskNums(t)}(k,:) = dataSx{j,t}(1:3,14);
                    cv2ipi{j,params.taskNums(t)}(k,:) = dataSx{j,t}(1:3,15);
                    cv2pkAmp{j,params.taskNums(t)}(k,:) = dataSx{j,t}(1:3,16);
                    cv2slWndwAmp{j,params.taskNums(t)}(k,:) = dataSx{j,t}(1:3,17);

                    % Extract frequency band data
                    for m = 1:3
                        nPsdSumBands{j,params.taskNums(t)}{m}(k,:) = dataSx{j,t}(row(1,m),2:5);
                        nPsdAucBands{j,params.taskNums(t)}{m}(k,:) = dataSx{j,t}(row(1,m),7:10);
                        nPsdAucPctBands{j,params.taskNums(t)}{m}(k,:) = dataSx{j,t}(row(1,m),12:15);

                        nPsdMaxChnks{j,params.taskNums(t)}{m}(k,:) = dataSx{j,t}(row(2,m),row(3,:));
                        nPsdMaxFreqChnks{j,params.taskNums(t)}{m}(k,:) = dataSx{j,t}(row(2,m),row(4,:));
                        freqChnks{j,params.taskNums(t)}{m}(k,:) = dataSx{j,t}(row(2,m),row(5,:));
                        pkAmpChnks{j,params.taskNums(t)}{m}(k,:) = dataSx{j,t}(row(2,m),row(6,:));
                        slWndwAmpChnks{j,params.taskNums(t)}{m}(k,:) = dataSx{j,t}(row(7,m),row(3,:));
                        ipiChnks{j,params.taskNums(t)}{m}(k,:) = dataSx{j,t}(row(7,m),row(4,:));
                        cvIpiChnks{j,params.taskNums(t)}{m}(k,:) = dataSx{j,t}(row(7,m),row(5,:));
                        cvPkAmpChnks{j,params.taskNums(t)}{m}(k,:) = dataSx{j,t}(row(7,m),row(6,:));
                        cvSlWndwAmpChnks{j,params.taskNums(t)}{m}(k,:) = dataSx{j,t}(row(8,m),row(3,:));
                        cv2ipiChnks{j,params.taskNums(t)}{m}(k,:) = dataSx{j,t}(row(8,m),row(4,:));
                        cv2pkAmpChnks{j,params.taskNums(t)}{m}(k,:) = dataSx{j,t}(row(8,m),row(5,:));
                    end
                end
            elseif params.taskNums(t) == 6
                dataSx{j,t} = readcell(fileName,'Sheet',params.names(t));
                nPsdMax{j,params.taskNums(t)}(k,:) = dataSx{j,t}(2:5,4);
                nPsdMaxFreq{j,params.taskNums(t)}(k,:) = dataSx{j,t}(2:5,5);
                nPsdAuc{j,params.taskNums(t)}(k,:) = dataSx{j,t}(1:5,6);
                nPsdAucPct{j,params.taskNums(t)}(k,:) = dataSx{j,t}(1:5,7);

                for m = 1:3
                    nPsdSumBands{j,params.taskNums(t)}(k,:) = dataSx{j,t}(row2(m),2:5);
                    nPsdAucBands{j,params.taskNums(t)}(k,:) = dataSx{j,t}(row2(m),7:10);
                    nPsdAucPctBands{j,params.taskNums(t)}(k,:) = dataSx{j,t}(row2(m),12:15);
                end
            end
        end

    end
    
    cd(groupList(1).folder)
end

%% Replace zeros with NaNs

baseVars = {nPsdMax,nPsdMaxFreq,nPsdAucPct};
% Average nPsd plots to create population mean
for w = 1:length(baseVars)
    for j = 1:length(groupList)
        for k = 1:length(nPsdMax)-1
            for i = 1:numel(baseVars{w}{j,k})
                if baseVars{w}{j,k}(i) == 0
                    baseVars{w}{j,k}(i) = NaN;
                end
            end
        end
    end
    baseVarsAvg{w} = cellfun(@(x) nanmean(x,1),baseVars{w}(:,1:length(nPsdMax)-1),'UniformOutput',false);
    baseVarsSe{w} = cellfun(@(x) nanstd(x,0,1)/(sqrt(size(x,1))),baseVars{w}(:,1:length(nPsdMax)-1),'UniformOutput',false);
end

for w = 1:length(baseVarsAvg)
    for k = 1:length(nPsdMax)-1
        for j = 1:length(groupList)
            for m = 1:3
                baseVarsBarAvg{m,w}(j,k) = baseVarsAvg{w}{j,k}(m);
                baseVarsBarSe{m,w}(j,k) = baseVarsSe{w}{j,k}(m);
                baseVarsBarSxs{m,w}{k,j}(:) = baseVars{w}{j,k}(:,m);
                xNums{m,w}{k,j} = length(baseVarsBarSxs{m,w}{k,j});
            end
        end
    end
end

%% 

colr = [0, 0, 0; 0.00,0.45,0.74];
c = [0.65,0.65,0.65; 0.46,0.73,0.91];
errCol = [0.8 0.5 0.5];
nums = [0.85 1.15; 1.85 2.15; 2.85 3.15; 3.85 4.15; 4.85 5.15];
titleNames = {'Ch X','Ch Y','Ch Z'};
subTitles = {'Normalized peak PSD','Normalized peak PSD freq. ','AUC (AUs)','AUC (percent of total)'};
for w = 1:length(baseVarsBarAvg)
    figure(); hold on
    sgtitle(subTitles{1,w},'FontSize',20)
    for  m = 1:3
        subplot(1,3,m); hold on
        b = bar(baseVarsBarAvg{m,w}', 'grouped');
        colororder(colr)
        title(titleNames{m})
        for t = 1:length(groupList)
            for k = 1:length(nPsdMax)-1
                xVals{k,t} = repelem(nums(k,t),xNums{m,w}{k,t});
                scatter(xVals{k,t},baseVarsBarSxs{m,w}{k,t},35,c(t,:),'LineWidth',1)
            end
        end
        errorbar(nums,baseVarsBarAvg{m,w}',baseVarsBarSe{m,w}','LineStyle','none','LineWidth',2,'Color',errCol)
        xticklabels(params.sheetNames)
        set(gca,'FontSize',16)
        if m == 1
            ylabel(subTitles{w})
        end
    end
    legend('Ctrl','Patients','FontSize',12)
end

%% Extended variables

expandedVars = {freq,pkAmp,rngeAmp,ipi,cvIpi,cvPkAmp,cvSlWndwAmp,cv2ipi,cv2pkAmp,cv2slWndwAmp};
for w = 1:length(expandedVars)
    for j = 1:length(groupList)
        for k = 1:length(nPsdMax)-2
            for i = 1:numel(expandedVars{w}{j,k})
                if expandedVars{w}{j,k}(i) == 0
                    expandedVars{w}{j,k}(i) = NaN;
                end
            end
        end
    end
    expandedVarsAvg{w} = cellfun(@(x) nanmean(x,1),expandedVars{w}(:,1:length(nPsdMax)-2),'UniformOutput',false);
    expandedVarsSe{w} = cellfun(@(x) nanstd(x,0,1)/(sqrt(size(x,1))),expandedVars{w}(:,1:length(nPsdMax)-2),'UniformOutput',false);
end

for w = 1:length(expandedVarsAvg)
    for k = 1:length(nPsdMax)-2
        for j = 1:length(groupList)
            for m = 1:3
                expandedVarsBarAvg{m,w}(j,k) = expandedVarsAvg{w}{j,k}(m);
                expandedVarsBarSe{m,w}(j,k) = expandedVarsSe{w}{j,k}(m);
                expandedVarsBarSxs{m,w}{k,j}(:) = expandedVars{w}{j,k}(:,m);
                xNums2{m,w}{k,j} = length(expandedVarsBarSxs{m,w}{k,j});
            end
        end
    end
end

%% Plot expanded variables

colr = [0, 0, 0; 0.00,0.45,0.74];
nums = [0.85 1.15; 1.85 2.15; 2.85 3.15; 3.85 4.15];
subTitles = {'Frequency (Hz)','Peak Amp. (g)','SlideWndw Amp (g)','Inter-peak interval (ms)','CV ipi','CV peak amp.'...
    'CV slWndw Amp','CV2 ipi','CV2 peak amp','CV2 slWndw Amp'};
for w = 1:length(expandedVarsBarAvg)
    figure(); hold on
    sgtitle(subTitles{1,w},'FontSize',20)
    for  m = 1:3
        subplot(1,3,m); hold on
        b = bar(expandedVarsBarAvg{m,w}', 'grouped');
        colororder(colr)
        title(titleNames{m})
        for t = 1:length(groupList)
            for k = 1:length(nPsdMax)-2
                xVals{k,t} = repelem(nums(k,t),xNums2{m,w}{k,t});
                scatter(xVals{k,t},expandedVarsBarSxs{m,w}{k,t}',35,c(t,:),'LineWidth',2)
            end
        end
        errorbar(nums,expandedVarsBarAvg{m,w}',expandedVarsBarSe{m,w}','LineStyle','none','LineWidth',2,'Color',errCol)
        xticklabels(params.sheetNames)
        set(gca,'FontSize',16)
        if m == 1
            ylabel(subTitles{w})
        end
    end
    legend('Ctrl','Patients','FontSize',12)
end

%% Summed frequency bands

summedVars = {nPsdSumBands,nPsdAucPctBands};
% Average nPsd plots to create population mean
for w = 1:length(summedVars)
    for j = 1:length(groupList)
        for k = 1:length(nPsdMax)-2
            for m = 1:3
                for i = 1:numel(summedVars{w}{j,k}{m})
                    if summedVars{w}{j,k}{m}(i) == 0
                        summedVars{w}{j,k}{m}(i) = NaN;
                    end
                end
            end
            sumVarsAvg{w}{j,k} = cellfun(@(x) nanmean(x,1),summedVars{w}{j,k}(:),'UniformOutput',false);
            sumVarsSe{w}{j,k} = cellfun(@(x) nanstd(x,0,1)/(sqrt(size(x,1))),summedVars{w}{j,k}(:),'UniformOutput',false);
        end
    end
end

for w = 1:length(sumVarsAvg)
    for k = 1:length(nPsdMax)-2
        for j = 1:length(groupList)
            for m = 1:3
                for t = 1:4
                    sumVarsBarAvg{m,w}{k}(j,t) = sumVarsAvg{w}{j,k}{m}(t);
                    sumVarsBarSe{m,w}{k}(j,t) = sumVarsSe{w}{j,k}{m}(t);
                    sumVarsBarSxs{m,w}{k,j}(t,:) = summedVars{w}{j,k}{m}(:,t);
                    xNums3{m,w}{k,j}(t) = size(sumVarsBarSxs{m,w}{k,j},2);
                end
            end
        end
    end
end

%% Plot summed frequencies

colr = [0, 0, 0; 0.00,0.45,0.74];
nums = [0.85 1.15; 1.85 2.15; 2.85 3.15; 3.85 4.15];
subTitles = {'Summed normalized PSD: ','Summed AUC (pct of total): '};
titleNames = {'Ch X','Ch Y','Ch Z'};
xNames = {'1-3 Hz','4-7 Hz','8-12 Hz','13-30 Hz'};

for w = 1:size(sumVarsBarAvg,2)
    for  m = 1:3
        figure(); hold on
        sgtitle(subTitles{w},'FontSize',20)
        for k = 1:length(nPsdMax)-2
            subplot(1,4,k); hold on
            b = bar(sumVarsBarAvg{m,w}{k}', 'grouped');
            colororder(colr)
            xticks([1 2 3 4])
            xticklabels(xNames)
            set(gca,'FontSize',16)
            for t = 1:length(groupList)
                for j = 1:4
                    xVals{k,t} = repelem(nums(j,t),xNums3{m,w}{k,t}(j));
                    scatter(xVals{k,t},sumVarsBarSxs{m,w}{k,t}(j,:),35,c(t,:),'LineWidth',1.5)
                end
            end
            errorbar(nums,sumVarsBarAvg{m,w}{k}',sumVarsBarSe{m,w}{k}','LineStyle','none','LineWidth',2,'Color',errCol)
            if k == 1
                ylabel(subTitles{w})
            end
            title(params.sheetNames{k})
        end
        sgtitle([subTitles{w} titleNames{m}])
        set(gca,'FontSize',16)
    end
    legend('Ctrl','Patients','FontSize',12)
end

%% Plot epoch data

chnkVars = {nPsdMaxChnks,nPsdMaxFreqChnks,freqChnks,pkAmpChnks,slWndwAmpChnks,ipiChnks,cvIpiChnks,cvPkAmpChnks,cvSlWndwAmpChnks,...
    cv2ipiChnks,cv2pkAmpChnks};

% Average nPsd plots to create population mean
for w = 1:length(chnkVars)
    for j = 1:length(groupList)
        for k = 1:length(nPsdMax)-2
            for m = 1:3
                for i = 1:numel(chnkVars{w}{j,k}{m})
                    if chnkVars{w}{j,k}{m}(i) == 0
                        chnkVars{w}{j,k}{m}(i) = NaN;
                    end
                end
            end
            chnkVarsAvg{w}{j,k} = cellfun(@(x) nanmean(x,1),chnkVars{w}{j,k}(:),'UniformOutput',false);
            chnkVarsSe{w}{j,k} = cellfun(@(x) nanstd(x,0,1)/(sqrt(size(x,1))),chnkVars{w}{j,k}(:),'UniformOutput',false);
        end
    end
end

% for w = 1:length(chnkVarsAvg)
%     for k = 1:length(nPsdMax)-2
%         for j = 1:length(groupList)
%             for m = 1:3
%                 for t = 1:params.numEpochs
%                     chnkVarsBarAvg{m,w}{k}(j,t) = chnkVarsAvg{w}{j,k}{m}(t);
%                     chnkVarsBarSxs{m,w}{k,j}(t,:) = chnkVars{w}{j,k}{m}(:,t);
%                     xNums4{m,w}{k,j}(t) = size(chnkVarsBarSxs{m,w}{k,j},2);
%                 end
%             end
%         end
%     end
% end

colr = [0, 0, 0; 0.00,0.45,0.74];
nums = [0.85 1.15; 1.85 2.15; 2.85 3.15; 3.85 4.15];
subTitles = {'Norm. peak PSD;','Norm. peak PSD freq.','Freq (by peaks)','Peak amplitude','SlWndw Amp','Ipi (ms)','CV ipi',...
    'CV peak amp.','CV slWndw Amp.','CV2 ipi','CV2 peak Amp.'};
titleNames = {'Ch X','Ch Y','Ch Z'};
xNames = {'Epoch 1','Epoch 2','Epoch 3'};

% for w = 1:size(chnkVarsBarAvg,2)
%     for  m = 1:3
%         figure(); hold on
%         sgtitle(subTitles{w},'FontSize',20)
%         for k = 1:length(nPsdMax)-2
%             subplot(1,4,k); hold on
%             b = bar(chnkVarsBarAvg{m,w}{k}', 'grouped');
%             colororder(colr)
%             xticks([1 2 3 4])
%             xticklabels(xNames)
%             set(gca,'FontSize',16)
%             title(params.sheetNames{k})
%             for t = 1:length(groupList)
%                 for j = 1:3
%                     xVals{k,t} = repelem(nums(j,t),xNums4{m,w}{k,t}(j));
%                     scatter(xVals{k,t},chnkVarsBarSxs{m,w}{k,t}(j,:),35,[0.93,0.69,0.13],'LineWidth',2)
%                 end
%             end
%             if k == 1
%                 ylabel(subTitles{w})
%             end
%         end
%         sgtitle([subTitles{w} titleNames{m}])
%         set(gca,'FontSize',16)
%     end
%     legend('Ctrl','Patients','FontSize',12)
% end

for w = 1:size(chnkVars,2)
    for k = 1:length(nPsdMax)-2
        figure(); hold on
        sgtitle([subTitles{w} params.sheetNames{k}],'FontSize',20)
        for t = 1:length(groupList)
            for m = 1:3
                subplot(1,3,m); hold on
                b = plot(chnkVarsAvg{w}{t,k}{m},'Color',colr(t,:),'LineWidth',2,'Marker','o','MarkerSize',10,'MarkerFaceColor',colr(t,:));
                errorbar(chnkVarsAvg{w}{t,k}{m},chnkVarsSe{w}{t,k}{m},'LineStyle','none','LineWidth',2,'Color',errCol)
                colororder(colr)
                xlim([0.5 3.5])
                xticks([1 2 3])
                xticklabels(xNames)
                title(titleNames{m})
                set(gca,'FontSize',16)
            end
            if k == 1
                ylabel(subTitles{w})
            end
        end
        legend('Ctrl','Patients','FontSize',1)
        set(gca,'FontSize',16)
    end
    legend('Ctrl','Patients','FontSize',12)
end

%% Save data

clear i j k t bound1 bound2 count x filColor lnCol
save(['compiled_' datestr(datetime('now'),'mmddyyyy')])
