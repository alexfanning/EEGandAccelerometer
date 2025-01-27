
clear; close all

compile = 'y';

if strcmpi(compile,'n')
    filename = uigetfile('*.wav');
    [data,sf] = audioread(filename);
    
    timeframe = 1;
    cohSec = 1;
    timeshift = 0;
    
    n = sf * timeframe;
    freqStep = sf / n; 
    n1s = sf;          
    nCoh = n1s * cohSec; 
    freqCohRes = sf / nCoh;
    shiftedframe = n1s * timeshift;
    
    signal = data;
    endNumber = fix(length(data) / n1s) - timeframe - ceil(timeshift);
    
    %% Compute PSD
    
    freqLabel = NaN(nCoh/2+1,endNumber);
    timeLabel = NaN(nCoh/2+1,endNumber);
    psd = NaN(nCoh/2+1,endNumber);
    
    for s = 1:endNumber
        
        y1 = signal((1 + (s-1) * n1s):(n+(s-1)*n1s));  % Each iteration grabs 20s of data and the starting point is shifted 1 s to the right of the prior starting point
    
        [PxxY, freqPsY] = pwelch(y1,hanning(nCoh),1/2*nCoh,nCoh,sf); % PSD by Wilch method and Hanning window
        psd(:,s) = PxxY;
    
        freqLabel(:,s) = freqPsY;
        timeLabel(1:(nCoh/2+1),s) = s-1;
    
    end
    
    %% Compute mean and normalized PSDs
    
    meanPsd = mean(psd,2);
    sumPsd = sum(meanPsd(1:45));
    normPsd = meanPsd / sumPsd;
    
    %% Extract peak PSD and frequency
    
    % Find maximum value
    [maxPsd,maxPsdFreq] = max(meanPsd(1:35));
    [nMaxPsd,nMaxPsdFreq] = max(normPsd(1:35));
    
    nPsdAll = sum(normPsd(1:35));
    bands = [1 3; 4 7; 8 12; 13 30];
    for t = 1:4
        nPsdBands(t) = sum(normPsd(bands(t,1):bands(t,2)));
        nPsdBandsPct(t) = (nPsdBands(t) / nPsdAll) * 100;
    end
    
%% Bandpass filter and convert to .wav file

%     % Bandpass filter audio data
%     bpData = bandpass(data,[600 1300],sf);
%     % File name for the .wav file
%     outputFileName = [filename(1:end-4) '_bp' '.wav'];
%     
%     % Write the raw data to a .wav file
%     audiowrite(outputFileName, bpData, sf);
    
    %% 

    figure(); hold on
    plot(meanPsd,'LineWidth',1)
    xlabel('Frequency')
    ylabel('Norm. PSD')
    
    figure(); hold on
    plot(data,'LineWidth',1)
    
    figure();
    gca3 = pcolor(timeLabel,freqLabel,psd);
    set(gca3, 'LineStyle','none');
    colorbar('location','eastoutside');
    clim([0 0.0000001]);
    title('PSD');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    axis([-inf,inf,0,50]);
    set(gca,'FontSize',16);
    
    
    %% Export data
    
    dataLabels = {'maxPsd';'maxPsdFreq';'normMaxPsd';'normMaxPsdFreq'};
    data2export = {maxPsd; maxPsdFreq; nMaxPsd; nMaxPsdFreq};
    
    dataLabels2 = {'nPsdBands';'nPsdBandsPct'};
    data2export2 = {nPsdBands;nPsdBandsPct};
    dataLabels4 = {'1-3Hz','4-7Hz','8-12Hz','13-30Hz'};
    
    dataLabels3 = {'meanPsd';'normPsd'};
    data2export3 = {meanPsd(1:2000);normPsd(1:8000)};
    
    writecell(data2export,[filename(1:end-4) '.xlsx'],'Sheet','Key metrics','Range','B2');
    writecell(data2export2,[filename(1:end-4) '.xlsx'],'Sheet','Key metrics','Range','B8');
    writecell(data2export3,[filename(1:end-4) '.xlsx'],'Sheet','Key metrics','Range','B11');
    
    writecell(dataLabels,[filename(1:end-4) '.xlsx'],'Sheet','Key metrics','Range','A2');
    writecell(dataLabels2,[filename(1:end-4) '.xlsx'],'Sheet','Key metrics','Range','A8');
    writecell(dataLabels3,[filename(1:end-4) '.xlsx'],'Sheet','Key metrics','Range','A11');
    writecell(dataLabels4,[filename(1:end-4) '.xlsx'],'Sheet','Key metrics','Range','B7');
    
    save([filename(1:end-4) '.mat']);
else

    [group,numGroups] = getFolders(3);
    for m = 1:numGroups
        cd((group(m).name))
        [sx,numSxs] = getFolders(3);

        for k = 1:numSxs
            cd (sx(k).name)
            filename = dir('*trial_2 -.xlsx');
            data = readcell(filename(1).name);
            data1 = readmatrix(filename(1).name);
            maxPsd(k,m) = data{1,2};
            maxPsdFreq(k,m) = data{2,2};
            normPsd(k,m) = data{3,2};
            normPsdFreq(k,m) = data{4,2};

            for w = 1:4
                nPsdBands{m}(k,w) = data{7,w+1};
                nPsdBandsPct{m}(k,w) = data{8,w+1};
            end
            meanPsdTrace{m}(k,:) = data1(5,2:2001);
            nPsdTrace{m}(k,:) = data1(6,2:2001);

            cd (sx(1).folder)
        end
        cd (group(1).folder)
    end

    meanPsdTracePopAvg = cellfun(@(x) nanmean(x,1),meanPsdTrace,'UniformOutput',false);
    nPsdTracePopAvg = cellfun(@(x) nanmean(x,1),nPsdTrace,'UniformOutput',false);

    figure(); hold on
    for i = 1:numGroups
        plot(meanPsdTracePopAvg{i})
    end

    for i = 1:numGroups
        figure(); hold on
        plot(meanPsdTrace{i}')
    end


end

save('trial2.mat')