%   Analyze accelerometer data
%
%   Written by Alex Fanning on 4/29/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;

% Create structures to hold data
params = struct();
data2export = struct();

% User sets parameters
tempInput = inputdlg({'Task number (1 = Action, 2 = Posture, 3 = Spiral, 4 = Tapping): ';'Group (1 = Control, 2 = Ataxia): ';'Subject number: '},'Parameters',1,{'4';'1';'1'});
params.taskNumber = str2double(tempInput{1});
params.group = str2double(tempInput{2});
params.sxNum = str2double(tempInput{3});
params.ch = [64,65,66; 38,41,48; 33,42,51; 44,47,49];
clear tempInput

% Automatic parameters
params.sf = 500;
params.wndwSize = 80;
params.thresh = 20;
params.sheetNames = {'Action','Posture','Spiral','Tapping','Rest'};
params.numEpochs = 3;

%% Accelerometer data extraction and processing

% Read the workspace of the subject
params.fileName = uigetfile('*.mat');
dataSx = load(params.fileName);
params.names = fieldnames(dataSx);

% Separate data based on condition
data = cell(1,5); vel = cell(1,5); pos = cell(1,5); eeg = cell(1,5);
data = taskList(params,data,dataSx,1);
clear dataSx

% Set more parameters
params.taskNums = find(~cellfun(@isempty,data));
params.prompt = {'Threshold (0.5-2.5): ','Diff order (1-3): '};
params.dlgtitle = 'Artifact removal threshold';
params.prompt2 = 'Artifact threshold good? (y/n): ';
params.prompt3 = {'Low: ','High: '};
params.defaults = {'3','30'};
params.base = {'1' '2'};
params.chunkSize = {};
totalInertia = struct();

% Loop through each task
for j = params.taskNumber% params.taskNums

    % User defines the bandpass range for hilbert processing
    if j ~= 5
        params.bpRange = str2double(inputdlg(params.prompt3,'BP range',1,params.defaults));
    end

    % Loop through each channel (X, Y, and Z)
    for m = 1:3

        eeg{1,j}(1,:) = detrend(data{1,j}(params.ch(2,1),:));
        eeg{1,j}(2,:) = detrend(data{1,j}(params.ch(2,2),:));
        eeg{1,j}(3,:) = detrend(data{1,j}(params.ch(2,3),:));
        eeg{2,j}(1,:) = detrend(data{1,j}(params.ch(3,1),:));
        eeg{2,j}(2,:) = detrend(data{1,j}(params.ch(3,2),:));
        eeg{2,j}(3,:) = detrend(data{1,j}(params.ch(3,3),:));
        eeg{3,j}(1,:) = detrend(data{1,j}(params.ch(4,1),:));
        eeg{3,j}(2,:) = detrend(data{1,j}(params.ch(4,2),:));
        eeg{3,j}(3,:) = detrend(data{1,j}(params.ch(4,3),:));

        % Detrend data
        data{2,j}(m,:) = detrend(data{1,j}(params.ch(1,m),:));
        data{3,j}(m,:) = data{2,j}(m,:);

%         vel{1,j}(m,:) = cumtrapz(0.1,data{2,j}(m,:));
%         pos{1,j}(m,:) = cumtrapz(1,vel{1,j}(m,:));

        % Remove artifact from tasks with mechanical noise
        params.continue = 'n';
        if j == 1 || j == 4
%             while params.continue == 'n'
%                 params.artThresh = inputdlg(params.prompt,params.dlgtitle,[1 40],params.base);
%                 params.diffOrder = str2double(params.artThresh{2});
%                 params.artThresh = str2double(params.artThresh{1});
                data = filtArtifact(data,params,[j m]);
%                 params.continue = input(params.prompt2,'s');
%             end
        end

        % Bandpass filter data for best PSD results
        data{4,j}(m,:) = bandpass(data{2,j}(m,:),[params.bpRange(1) params.bpRange(2)],params.sf);
        data{5,j}(m,:) = bandpass(data{2,j}(m,:),[1 30],params.sf);

        % Apply hilbert transform to find instantaneous frequency
        hbert{1,j}(m,:) = hilbert(data{4,j}(m,:));
        hbert{2,j}(m,:) = angle(hbert{1,j}(m,:));
        hbert{3,j}(m,:) = sin(hbert{2,j}(m,:));

        % Find average of each channel (X, Y, and Z)
        chAvg(j,m) = mean(data{2,j}(m,:),2);

        % Calculate inertia for each channel (X, Y, and Z)
        inertia{1,j}(m,:) = data{2,j}(m,:) - chAvg(j,m);

        close all

        %data{2,2}(2,:) = eeg{2,1}(2,:);
        if j ~= 5
             % Find peaks of oscillations
             [data,data2export,params] = findPks(data,params,data2export,[j m]);
        else
            data2export.idxs{j,m}(1) = data2export.idxs{1,m}(1);
        end

        % Calculate size of epochs to break up data into chunks
        temp = length(data{2,j}) - data2export.idxs{j,m}(1);
        params.rem = rem((temp / params.sf),params.numEpochs) * params.sf;
        params.chunkSize{j} = [((temp - params.rem) / params.sf) / params.numEpochs, ((temp - params.rem) / params.sf) / params.numEpochs * 2, (temp - params.rem) / params.sf];

        % Calculate power and peak frequency for action tasks
        [data2export,params] = computePSD(data{5,j}(m,:),params,data2export,1,1,5,3,[j m]);

        if j ~= 5

            % Find frequency of peaks for each channel
            data2export.freq{j}(m) = numel(data2export.peaks{j,m}) / (temp / params.sf);
            clear temp

            % Calculate average of peak amplitudes for each channel
            data2export.pkAmpAvg{j}(m) = mean(data2export.peaks{j,m});

            % Calculate amplitude of signal with sliding window
            [params,data2export] = ampSlWndw(data,params,data2export,[j m]);

            % Find inter-peak intervals for each channel (time in ms)
            for t = 2:length(data2export.idxs{j,m})
                data2export.ipi{j,m}(t-1) = ((data2export.idxs{j,m}(t) - data2export.idxs{j,m}(t-1)) / params.sf) * 1000;
            end

            % Find average inter-peak interval for each channel
            data2export.ipiAvg{j}(m) = mean(data2export.ipi{j,m});

            % Find coefficient of variation for inter-peak interval and amplitude
            data2export.cvIpi{j}(m) = std(data2export.ipi{j,m}) / mean(data2export.ipi{j,m});
            data2export.cvPkAmp{j}(m) = std(data2export.peaks{j,m}) / nanmean(data2export.peaks{j,m});
            data2export.cvSlWndwAmp{j}(m) = std(data2export.wndwRange{j,m}) / mean(data2export.wndwRange{j,m});

            for t = 2:length(data2export.ipi{j,m})
                data2export.allCv2ipi{j,m}(t-1) = (std([data2export.ipi{j,m}(t-1) data2export.ipi{j,m}(t)]) / mean([data2export.ipi{j,m}(t-1) data2export.ipi{j,m}(t)])) * sqrt(2);
                data2export.allCv2pkAmp{j,m}(t-1) = (std([data2export.peaks{j,m}(t-1) data2export.peaks{j,m}(t)]) / mean([data2export.peaks{j,m}(t-1) data2export.peaks{j,m}(t)])) * sqrt(2);
            end
            data2export.cv2ipi{1,j}(m) = mean(data2export.allCv2ipi{j,m});
            data2export.cv2pkAmp{1,j}(m) = mean(data2export.allCv2pkAmp{j,m});

            for t = 2:length(data2export.wndwRange{j,m})
                data2export.allCv2slWndwAmp{j,m}(t-1) = (std([data2export.wndwRange{j,m}(t-1) data2export.wndwRange{j,m}(t)]) / mean([data2export.wndwRange{j,m}(t-1) data2export.wndwRange{j,m}(t)])) * sqrt(2);
            end
            data2export.cv2slWndwAmp{j}(m) = nanmean(data2export.allCv2slWndwAmp{j,m});

            % Sort peaks and indexes by chunk windows
            params.chnkIdx{j}(m,1) = 1;
            params.ipiChnkIdx{j}(m,1) = 1;
            for t = 2:numel(params.chunkSize{j})+1

                    params.chnkIdx{j}(m,t) = interp1(data2export.idxs{j,m}/params.sf,data2export.idxs{j,m}/params.sf,(params.chunkSize{j}(t-1)+data2export.idxs{j,m}(1)/params.sf),'nearest');
                    if ~isnan(params.chnkIdx{j}(m,t))
                        params.chnkIdx{j}(m,t) = find(params.chnkIdx{j}(m,t)==(data2export.idxs{j,m}/params.sf));
                        params.ipiChnkIdx{j}(m,t) = params.chnkIdx{j}(m,t)-1;
                    else
                        params.chnkIdx{j}(m,t) = find(data2export.idxs{j,m}/params.sf,1,'last');
                        params.ipiChnkIdx{j}(m,t) = params.chnkIdx{j}(m,t)-1;
                    end

                    data2export.pkAmpChnkAvg{j}(m,t-1) = mean(data2export.peaks{j,m}(params.chnkIdx{j}(m,t-1):params.chnkIdx{j}(m,t)));
                    data2export.ipiChnkAvg{j}(m,t-1) =  mean(data2export.ipi{j,m}(params.ipiChnkIdx{j}(m,t-1):params.ipiChnkIdx{j}(m,t)));
                    data2export.cvIpiChnkAvg{j}(m,t-1) = std(data2export.ipi{j,m}(params.ipiChnkIdx{j}(m,t-1):params.ipiChnkIdx{j}(m,t))) / mean(data2export.ipi{j,m}(params.ipiChnkIdx{j}(m,t-1):params.ipiChnkIdx{j}(m,t)));
                    data2export.cvPkAmpChnkAvg{j}(m,t-1) = std(data2export.peaks{j,m}(params.chnkIdx{j}(m,t-1):params.chnkIdx{j}(m,t))) / mean(data2export.peaks{j,m}(params.chnkIdx{j}(m,t-1):params.chnkIdx{j}(m,t)));
                    data2export.cvPkAmpChnkAvg{j}(m,t-1) = std(data2export.peaks{j,m}(params.chnkIdx{j}(m,t-1):params.chnkIdx{j}(m,t))) / mean(data2export.peaks{j,m}(params.chnkIdx{j}(m,t-1):params.chnkIdx{j}(m,t)));
                    data2export.freqChnk{j}(m,t-1) = numel(data2export.peaks{j,m}(params.chnkIdx{j}(m,t-1):params.chnkIdx{j}(m,t))) / params.chunkSize{j}(1);
                    for tt = params.ipiChnkIdx{j}(m,t-1)+1:params.ipiChnkIdx{j}(m,t)
                        params.tempCv2chnk{j}{m,t}(tt-1) = (std([data2export.ipi{j,m}(tt-1) data2export.ipi{j,m}(tt)]) / nanmean([data2export.ipi{j,m}(tt-1) data2export.ipi{j,m}(tt)])) * sqrt(2);
                        params.tempCv2pkAmpChnk{j}{m,t}(tt-1) = (std([data2export.peaks{j,m}(tt-1) data2export.peaks{j,m}(tt)]) / nanmean([data2export.peaks{j,m}(tt-1) data2export.peaks{j,m}(tt)])) * sqrt(2);
                    end
                    data2export.cv2ipiChnk{j}(m,t-1) = nanmean(params.tempCv2chnk{j}{m,t});
                    data2export.cv2pkAmpChnk{j}(m,t-1) = nanmean(params.tempCv2pkAmpChnk{j}{m,t});
                    if t-1 == 1
                        data2export.wndwRangeChnk{j}(m,t-1) = nanmean(data2export.wndwRange{j,m}(1:params.slWndwChnk{j}(t-1)));
                        data2export.cvSlWndwAmpChnkAvg{j}(m,t-1) = std(data2export.wndwRange{j,m}(1:params.slWndwChnk{j}(t-1))) / nanmean(data2export.wndwRange{j,m}(1:params.slWndwChnk{j}(t-1)));
                    else
                        data2export.wndwRangeChnk{j}(m,t-1) = nanmean(data2export.wndwRange{j,m}(params.slWndwChnk{j}(t-2):params.slWndwChnk{j}(t-1)));
                        data2export.cvSlWndwAmpChnkAvg{j}(m,t-1) = std(data2export.wndwRange{j,m}(params.slWndwChnk{j}(t-2):params.slWndwChnk{j}(t-1))) / nanmean(data2export.wndwRange{j,m}(params.slWndwChnk{j}(t-2):params.slWndwChnk{j}(t-1)));
                    end

            end

            % reshape data to align each cycle to peak
            data2export.pkAlign{j,m} = NaN(15000,10000);
            for p = 1:3
                data2export.eegAlign{p}{j,m} = NaN(15000,10000);
            end
            params.pkAlignShft = round((1/data2export.freq{j}(m)*params.sf)/2);
            for t = 2:length(data2export.idxs{j,m})
                if data2export.idxs{j,m}(t-1) - params.pkAlignShft > 0
                    data2export.pkAlign{j,m}(1:length(data2export.idxs{j,m}(t-1):data2export.idxs{j,m}(t)),t-1) = data{2,j}(m,data2export.idxs{j,m}(t-1) - params.pkAlignShft:data2export.idxs{j,m}(t) - params.pkAlignShft);
                    for p = 1:3
                        data2export.eegAlign{p}{j,m}(1:length(data2export.idxs{j,m}(t-1):data2export.idxs{j,m}(t)),t-1) = eeg{p,j}(m,data2export.idxs{j,m}(t-1) - params.pkAlignShft:data2export.idxs{j,m}(t) - params.pkAlignShft);
                    end
                end
            end
            data2export.pkAlignAvg{j}(:,m) = nanmean(data2export.pkAlign{j,m},2);
            for p = 1:3
                data2export.eegAlignAvg{p,j}(:,m) = nanmean(data2export.eegAlign{p}{j,m},2);
            end

            figure(); hold on
            plot(data2export.pkAlign{j,m})
            plot(data2export.pkAlignAvg{j}(:,m),'k','LineWidth',3)
            xline(params.pkAlignShft,'--',{'Peak'},'LineWidth',1.5)
            xlim([0 500])
            xlabel('Time (ms)')
            ylabel('Acceleration (g)')

            figure(); hold on
            plot(data2export.eegAlign{p}{j,m})
            plot(data2export.eegAlignAvg{p,j}(:,m),'k','LineWidth',3)
            xlim([0 500])
        end

    end

    % Calculate total inertia and plot
    totalInertia.data{j} = (inertia{j}(1,:) + inertia{j}(2,:) + inertia{j}(3,:))';

    if j ~= 5
        [data2export.amp{1,j},params.bestCh{1,j}] = max(data2export.rangeAvg{j});
        [data2export.amp{2,j},params.bestCh{2,j}] = max(data2export.rangeAvg{j});
    end

    % Calculate power and peak frequency for action tasks for total inertia
    totalInertia.idxs{j,1}(1) = data2export.idxs{j,m}(1);
    [totalInertia,params] = computePSD(totalInertia.data{j},params,totalInertia,1,1,5,3,[j 1]);

end
close all
%% 

eegPsd = struct();
for k = 1:3
    count = 1;
    for m = 1:3
        eegPsd.idxs{j,m}(1) = data2export.idxs{j,m}(1);
        if sum(eeg{k,j}(m,:))==0
            continue
        else
            [eegPsd,params] = computePSD(eeg{k,j}(m,:),params,eegPsd,1,1,5,3,[j m]);
        end
        nPsdBands{k} = eegPsd.nPsdBands{j};
        nPsdBandsAucPct{k} = eegPsd.nPsdBandsAucPct{j};
        nMaxPsd{k} = eegPsd.nMaxPsd{j};
        nMaxPsdFreq{k} = eegPsd.nMaxPsdFreq{j};
        count = count + 1;
    end
    close all
end

nPsdBands = cellfun(@(x) subsasgn(x, substruct('()', {x==0}), NaN), nPsdBands, 'uniform', 0);
nPsdBandsAucPct = cellfun(@(x) subsasgn(x, substruct('()', {x==0}), NaN), nPsdBandsAucPct, 'uniform', 0);
nMaxPsd = cellfun(@(x) subsasgn(x, substruct('()', {x==0}), NaN), nMaxPsd, 'uniform', 0);
nMaxPsdFreq = cellfun(@(x) subsasgn(x, substruct('()', {x==0}), NaN), nMaxPsdFreq, 'uniform', 0);

%% 


params.xlCellId = {['B' num2str(params.sxNum)],['H' num2str(params.sxNum)]; ['B' num2str(params.sxNum+24)],['H' num2str(params.sxNum+24)];...
    ['B' num2str(params.sxNum+49)],['C' num2str(params.sxNum+49)]; ['F' num2str(params.sxNum+49)],['G' num2str(params.sxNum+49)]};
params.xlCellId2 = {['A' num2str(1)];['A' num2str(2)];['A' num2str(25)];['A' num2str(50)];['E' num2str(50)];['G' num2str(1)]};
params.xlTabNames = {'Ch 38','Ch 41','Ch 48','Ch 33','Ch 42','Ch 51'};
params.xlTabNames2 = {'Ch X','Ch Y','Ch Z'};
params.dataNames = {'Ctrl';'nPsdBands';'nPsdBandsAucPct';'nMaxPsd';'nMaxPsdFreq';'Ataxia'};
count2 = 1;
for k = 1:3
    for m = 1:3
        writematrix(nPsdBands{k}(m,:),[params.sheetNames{j} '_eegSummary.xlsx'],'Sheet',params.xlTabNames{count2},'Range',params.xlCellId{1,params.group})
        writematrix(nPsdBandsAucPct{k}(m,:),[params.sheetNames{j} '_eegSummary.xlsx'],'Sheet',params.xlTabNames{count2},'Range',params.xlCellId{2,params.group})
        writematrix(nMaxPsd{k}(m),[params.sheetNames{j} '_eegSummary.xlsx'],'Sheet',params.xlTabNames{count2},'Range',params.xlCellId{3,params.group})
        writematrix(nMaxPsdFreq{k}(m),[params.sheetNames{j} '_eegSummary.xlsx'],'Sheet',params.xlTabNames{count2},'Range',params.xlCellId{4,params.group})
        for d = 1:6
            writecell(cellstr(params.dataNames{d}),[params.sheetNames{j} '_eegSummary.xlsx'],'Sheet',params.xlTabNames{count2},'Range',params.xlCellId2{d})
        end
        count2 = count2 + 1;
    end
end
for k = 1:3
    writematrix(data2export.nPsdBands{j}(k,:),[params.sheetNames{j} '_eegSummary.xlsx'],'Sheet',params.xlTabNames2{k},'Range',params.xlCellId{1,params.group})
    writematrix(data2export.nPsdBandsAucPct{j}(k,:),[params.sheetNames{j} '_eegSummary.xlsx'],'Sheet',params.xlTabNames2{k},'Range',params.xlCellId{2,params.group})
    writematrix(data2export.nMaxPsd{j}(k),[params.sheetNames{j} '_eegSummary.xlsx'],'Sheet',params.xlTabNames2{k},'Range',params.xlCellId{3,params.group})
    writematrix(data2export.nMaxPsdFreq{j}(k),[params.sheetNames{j} '_eegSummary.xlsx'],'Sheet',params.xlTabNames2{k},'Range',params.xlCellId{4,params.group})
    for d = 1:6
        writecell(cellstr(params.dataNames{d}),[params.sheetNames{j} '_eegSummary.xlsx'],'Sheet',params.xlTabNames2{k},'Range',params.xlCellId2{d})
    end
end
% plot(eeg{2,1}(1,:));

%% Coherency

% data{2,2}(2,:) = eeg{2,4}(3,:);
% data{2,3}(2,:) = eeg{2,4}(1,:);
% [data2export] = eegCohere(data,data2export,params,2,1);
% col = copper(3);
% figure(); hold on
% for i = 1:3
%     plot(data2export.mCoh(:,i),'LineWidth',2,'Color',col(i,:))
% end
% ylabel('Coherence')
% xlabel('Frequency (Hz)')
% xlim([0 30])
% set(gca,'FontSize',16)
% legend(params.titleNamesCh)

%% Transfer data to excel

xlVarNames = {'maxPSD','maxPsdFreq','nMaxPsd','nMaxPsdFreq','nPsdAuc','nPsdAucPct','freq (by pks)','pkAmpAvg','rangeAvg','ipi (ms)','cvIpi','cvPkAmp','cvSlWndwAmp',...
    'cv2ipi','cv2pkAmp','cv2slWndwAmp'};
xlVarNames2 = {'nPsd (sum)','nPsdAuc (sum)','nPsdAucPct'};
xlVarNames3 = {'1-3Hz','4-7Hz','8-12Hz','13-30Hz','','1-3Hz','4-7Hz','8-12Hz','13-30Hz','','1-3Hz','4-7Hz','8-12Hz','13-30Hz'};
xlVarNames4 = {'Ch X';'Ch Y';'Ch Z'};
xlVarNames5 = {'Epoch 1','Epoch 2','Epoch 3'};
xlVarNames6 = {'nPsdMax','nPsdMaxFreq','freq','pkAmp','slWndwAmp','ipi (ms)','cvIpi','cvPkAmp','cvSlWndwAmp','cv2ipi','cv2pkAmp','cv2slWndwAmp'};
% xlVarNames7 = {'nPsd','totalInertia'};
data2export.xlVars = {data2export.maxPSD,data2export.maxPsdFreq,data2export.nMaxPsd,data2export.nMaxPsdFreq,data2export.nPsdAuc,data2export.nPsdAucPct,...
    data2export.freq,data2export.pkAmpAvg,data2export.rangeAvg,data2export.ipiAvg,data2export.cvIpi,data2export.cvPkAmp,data2export.cvSlWndwAmp,...
    data2export.cv2ipi,data2export.cv2pkAmp,data2export.cv2slWndwAmp};
data2export.xlVars2 = {data2export.nPsdBands,data2export.nPsdBandAuc,data2export.nPsdBandsAucPct};
data2export.xlVars3 = {data2export.nPsdChnkMax,data2export.nPsdChnkMaxFreq,data2export.freqChnk,data2export.pkAmpChnkAvg,data2export.wndwRangeChnk,data2export.ipiChnkAvg,...
    data2export.cvIpiChnkAvg,data2export.cvPkAmpChnkAvg,data2export.cvSlWndwAmpChnkAvg,data2export.cv2ipiChnk,data2export.cv2pkAmpChnk};
data2export.xlVars4 = {totalInertia.maxPSD,totalInertia.maxPsdFreq,totalInertia.nMaxPsd,data2export.nMaxPsdFreq,totalInertia.nPsdAuc,totalInertia.nPsdAucPct};
data2export.xlVars5 = {totalInertia.nPsdBands,totalInertia.nPsdBandAuc,totalInertia.nPsdBandsAucPct};
% data2export.xlVars6 = {data2export.nMeanPsd,data2export.inertiaTotal'};

params.xlRng = {'B2','C2','D2','E2','F2','G2','H2','I2','J2','K2','L2','M2','N2','O2','P2','Q2'};
params.xlRng2 = {'B8','G8','L8','B6','G6','L6'};
params.xlRng3 = {'B14','F14','J14','N14','B19','F19','J19','N19','B24','F24','J24','B12','F12','J12','N12','B18','F18','J18','N18','B23','F23','J23'};
params.xlRng4 = {'A2','A8','A14','A19','A24'};
params.xlRng5 = {'B13','F13','J13','N13'};
params.xlRng6 = {'B10','G10','L10','B8','G8','L8'};
% params.xlRng7 = {'W2','AA2'};
for j = params.taskNums
    if j ~= 5
        for k = 1:length(data2export.xlVars)
            writematrix(data2export.xlVars{k}{j}',[params.fileName(1:end-4) '_accSummary.xlsx'],'Sheet',params.sheetNames{j},'Range',params.xlRng{k})
        end
    
        for k = 1:length(data2export.xlVars2)
            writematrix(data2export.xlVars2{k}{j},[params.fileName(1:end-4) '_accSummary.xlsx'],'Sheet',params.sheetNames{j},'Range',params.xlRng2{k})
            writecell(cellstr(xlVarNames2{k}),[params.fileName(1:end-4) '_accSummary.xlsx'],'Sheet',params.sheetNames{j},'Range',params.xlRng2{k+3})
        end
        
        for k = 1:length(data2export.xlVars3)
            writematrix(data2export.xlVars3{k}{j},[params.fileName(1:end-4) '_accSummary.xlsx'],'Sheet',params.sheetNames{j},'Range',params.xlRng3{k})
            writecell(cellstr(xlVarNames6(k)),[params.fileName(1:end-4) '_accSummary.xlsx'],'Sheet',params.sheetNames{j},'Range',params.xlRng3{k+length(data2export.xlVars3)})
        end
    
        writecell(cellstr(xlVarNames),[params.fileName(1:end-4) '_accSummary.xlsx'],'Sheet',params.sheetNames{j},'Range','B1')
        writecell(cellstr(xlVarNames3),[params.fileName(1:end-4) '_accSummary.xlsx'],'Sheet',params.sheetNames{j},'Range','B7')

        for k = 1:length(params.xlRng5)
            writecell(cellstr(xlVarNames5),[params.fileName(1:end-4) '_accSummary.xlsx'],'Sheet',params.sheetNames{j},'Range',params.xlRng5{k})
        end
        for k = 1:5
            writecell(cellstr(xlVarNames4),[params.fileName(1:end-4) '_accSummary.xlsx'],'Sheet',params.sheetNames{j},'Range',params.xlRng4{k})
        end
%         for k = 1:length(xlVarNames7)
%             writematrix(data2export.xlVars6{k}{j},[params.fileName(1:end-4) '_accSummary.xlsx'],'Sheet',params.sheetNames{j},'Range',params.xlRng6{k})
%         end

        for k = 1:length(data2export.xlVars4)
            writecell(data2export.xlVars4{k}',[params.fileName(1:end-4) '_accSummary.xlsx'],'Sheet','totalInertia','Range',params.xlRng{k})
        end
        writecell(cellstr(xlVarNames(1:6)),[params.fileName(1:end-4) '_accSummary.xlsx'],'Sheet','totalInertia','Range','B1')
        writecell(cellstr(xlVarNames3),[params.fileName(1:end-4) '_accSummary.xlsx'],'Sheet','totalInertia','Range','B9')
        writecell(cellstr(params.sheetNames'),[params.fileName(1:end-4) '_accSummary.xlsx'],'Sheet','totalInertia','Range','A2')

        for k = 1:length(data2export.xlVars5)
            writecell(data2export.xlVars5{k}',[params.fileName(1:end-4) '_accSummary.xlsx'],'Sheet','totalInertia','Range',params.xlRng6{k})
            writecell(cellstr(xlVarNames2{k}),[params.fileName(1:end-4) '_accSummary.xlsx'],'Sheet','totalInertia','Range',params.xlRng6{k+3})
        end
        writecell(cellstr(params.sheetNames'),[params.fileName(1:end-4) '_accSummary.xlsx'],'Sheet','totalInertia','Range','A10')
    else
        for k = 1:6
            writematrix(data2export.xlVars{k}{j}',[params.fileName(1:end-4) '_accSummary.xlsx'],'Sheet',params.sheetNames{j},'Range',params.xlRng{k})
        end
        writecell(cellstr(xlVarNames(1:6)),[params.fileName(1:end-4) '_accSummary.xlsx'],'Sheet',params.sheetNames{j},'Range','B1')
        writecell(cellstr(xlVarNames4),[params.fileName(1:end-4) '_accSummary.xlsx'],'Sheet',params.sheetNames{j},'Range','A2')
    end
end

%% Save workspace

% clear m t tt ct k p xlVarNames xlVarNames2 xlVarNames xlVarNames3 xlVarNames4 xlVarNames5 xlVarNames6 xlVarNames7
% params.tempName = datetime('now','Format','dd/MM/yyyy');
% params.tempName = datestr(params.tempName,'ddmmmyy');
save([params.fileName(1:end-4) '_' params.sheetNames{j} '.mat'])
