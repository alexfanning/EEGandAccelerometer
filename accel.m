%   Analyze accelerometer data
%
%   Written by Alex Fanning on 2/16/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;

params = struct();
data2export = struct();

% Get folder structure
[groupList,params.numFolders] = getFolders(0);

%% Accelerometer data extraction and processing

% Set parameters
params.sf = 500;
params.ch = [64,65,66];
params.wndwSize = 80;
params.thresh = 20;
params.sheetNames = {'Action','Posture','Spiral','Tapping','Rest'};

% Loop through each group
for i = 1:params.numFolders
    
    cd(groupList(1).folder)
    cd(groupList(i).name)
    params.fileList = dir('*.mat');
    params.numSxs(i) = length(params.fileList);

    % Loop through each subject and grab their data
    for ii = 1:params.numSxs(i)
        params.fileName = params.fileList(ii).name;
        dataSx = load(params.fileName);
        params.recName = params.fileName;
        params.names = fieldnames(dataSx);

        % Separate data based on condition
        data{1,i}{ii,1} = dataSx.Action.F;
        data{1,i}{ii,2} = dataSx.Posture.F;
        data{1,i}{ii,3} = dataSx.Spiral.F;
        data{1,i}{ii,4} = dataSx.Tapping.F;
        data{1,i}{ii,5} = dataSx.Rest.F;

        % Loop through each task
        for j = 1:length(params.names)
    
            % Calculate inertia for each channel (X, Y, and Z)
            for m = 1:3
                
                % Detrend data
                data{2,i}{ii,j}(m,:) = detrend(data{1,i}{ii,j}(params.ch(m),:));
                data{3,i}{ii,j}(m,:) = data{2,i}{ii,j}(m,:);

                if j == 1 || j == 4
                    data = filtArtifact(data,params,[i ii j m]);
                end

                data{4,i}{ii,j}(m,:) = hilbert(data{2,i}{ii,j}(m,:));

                 % Find average of each channel (X, Y, and Z)
                data2export.dataAvg{i}{ii,j}(m) = mean(data{2,i}{ii,j}(m,:),2);

                % Calculate inertia for each channel (X, Y, and Z)
                data2export.inertia{i}{ii,j}(m,:) = data{2,i}{ii,j}(m,:) - data2export.dataAvg{i}{ii,j}(m);
    
                % Bandpass filter data for best PSD results
                data2export.bp{1,i}{ii,j}(m,:) = bandpass(data{2,i}{ii,j}(m,:),[1 30],params.sf);

                % Calculate power and peak frequency
                [data2export,params] = computePSD(data2export.bp{i}{ii,j}(m,:),params,data2export,1,1,5,3,[i ii j m]);

                params.paramFix = 'n';
                params.minPeakDist =  10;
                if j == 1
                    params.wndwSize =  250;
                    params.minPeakProm = 0.00005;
                elseif j == 2
                    params.wndwSize =  40;
                    params.minPeakProm = 0.000005;
                elseif j == 3
                    params.wndwSize =  40;
                    params.minPeakProm = 0.000005;
                elseif j == 4
                    params.wndwSize =  15;
                    params.minPeakProm = 0.0001;
                else
                    params.paramFix = 'y';
                end

                while params.paramFix == 'n'

                    % Smooth data
                    data{5,i}{ii,j}(m,:) = smooth(data{2,i}{ii,j}(m,:),params.wndwSize);
                    % data{5,i}{ii,j}(m,:) = fitLine(params,data{2,i}{ii,j}(m,:));

                    [data2export.pks{i}{ii,j}{m},data2export.idxs{i}{ii,j}{m},data2export.w{i}{ii,j}{m},data2export.p{i}{ii,j}{m}] = findpeaks(data{5,i}{ii,j}(m,:),1,'MinPeakProminence',params.minPeakProm,'MinPeakDistance',params.minPeakDist);

                    params.titleNames = {'Action','Posture','Spiral','Tapping','Rest',};
                    params.titleNamesCh = {'Ch X','Ch Y','Ch Z'};
                    figure(); hold on
                    plot(data{2,i}{ii,j}(m,:))
                    plot(data{5,i}{ii,j}(m,:))
                    scatter(data2export.idxs{i}{ii,j}{m},data2export.pks{i}{ii,j}{m})
                    xlabel('Sample number')
                    ylabel('Acceleration (g)')
                    title([params.titleNames{j} params.titleNamesCh{m} ': Raw data'])
                    set(gca,'FontSize',16)

                    prompt2 = 'Analysis good? (y/n): ';
                    params.paramFix = input(prompt2,"s");

                    if params.paramFix == 'n'
                        prompt = {'minPeakProm: ','minPeakInt: ','wndwSize: '};
                        dlgtitle = 'FindPeaks params';
                        if j == 1
                            default = {'0.00005','10','250'};
                        elseif j == 2
                            default = {'0.000005','10','40'};
                        elseif j == 3
                            default = {'0.00005','10','150'};
                        elseif j == 4
                            default = {'0.0001','10','15'};
                        end
                        tempParams = inputdlg(prompt,dlgtitle,1,default);
                        params.minPeakProm = str2double(tempParams{1});
                        params.minPeakDist =  str2double(tempParams{2});
                        params.wndwSize =  str2double(tempParams{3});
                    end
                end

                close all

                if j ~= 5

                    % Calculate size of epochs to break up data into chunks
                    params.temp = length(data{2,i}{ii,j}) - data2export.idxs{i}{ii,j}{m}(1);
                    params.rem = rem((params.temp / params.sf),3) * params.sf;
                    params.chunkSize{i}{ii,j} = [((params.temp - params.rem) / params.sf) / 3, ((params.temp - params.rem) / params.sf) / 3 * 2, (params.temp - params.rem) / params.sf];

                    % Find frequency of peaks for each channel
                    data2export.freq{i}{ii,j}(m) = numel(data2export.pks{i}{ii,j}{m}) / (params.temp / params.sf);

                    % Calculate average of peak amplitudes for each channel
                    data2export.pkAmpAvg{i}{ii,j}(m) = mean(data2export.pks{i}{ii,j}{m});

                    % Calculate amplitude of signal with sliding window
                    params.slWndwLng = floor(params.chunkSize{i}{ii,j}(end) / (1 / data2export.freq{i}{ii,j}(m)));
                    params.wndwTime = floor((1 / data2export.freq{i}{ii,j}(m)) * params.sf);
                    params.slWndwStart = ceil((data2export.idxs{i}{ii,j}{m}(1)/params.sf) / (1 / data2export.freq{i}{ii,j}(m)));
                    for t = params.slWndwStart:params.slWndwLng
                        data2export.wndwRange{i}{ii,j}(m,t-params.slWndwStart+1) = abs(range(data{2,i}{ii,j}(m,params.wndwTime*(t-1)+1:params.wndwTime*(t))));
                    end

                    data2export.rangeAvg{i}{ii,j}(m) = mean(data2export.wndwRange{i}{ii,j}(m,:),2);

                    params.slWndwSpan = params.slWndwLng - params.slWndwStart;
                    params.slWndwBlkLng = floor(params.slWndwSpan / 3);
                    params.slWndwChnk = [params.slWndwBlkLng, params.slWndwBlkLng*2, params.slWndwBlkLng*3];

                    % Find inter-peak intervals for each channel
                    for t = 2:length(data2export.idxs{i}{ii,j}{m})
                        data2export.ipi{i}{ii,j}{m}(t-1) = ((data2export.idxs{i}{ii,j}{m}(t) - data2export.idxs{i}{ii,j}{m}(t-1)) / params.sf) * 1000;
                    end

                    % Find average inter-peak interval for each channel
                    data2export.ipiAvg{i}{ii,j}(m) = mean(data2export.ipi{i}{ii,j}{m});
    
                    % Find coefficient of variation for inter-peak interval
                    data2export.cvIpi{i}{ii,j}(m) = std(data2export.ipi{i}{ii,j}{m}) / mean(data2export.ipi{i}{ii,j}{m});
                    data2export.cvPkAmp{i}{ii,j}(m) = std(data2export.pks{i}{ii,j}{m}) / mean(data2export.pks{i}{ii,j}{m});
    
                    for t = 2:length(data2export.ipi{i}{ii,j}{m})
                        params.tempCv2ipi{i}{ii,j}{m}(t-1) = (std([data2export.ipi{i}{ii,j}{m}(t-1) data2export.ipi{i}{ii,j}{m}(t)]) / mean([data2export.ipi{i}{ii,j}{m}(t-1) data2export.ipi{i}{ii,j}{m}(t)])) * sqrt(2);
                        params.tempCv2pkAmp{i}{ii,j}{m}(t-1) = (std([data2export.pks{i}{ii,j}{m}(t-1) data2export.pks{i}{ii,j}{m}(t)]) / mean([data2export.pks{i}{ii,j}{m}(t-1) data2export.pks{i}{ii,j}{m}(t)])) * sqrt(2);
                    end
                    data2export.cv2ipi{i}{ii,j}(m) = mean(params.tempCv2ipi{i}{ii,j}{m});
                    data2export.cv2pkAmp{i}{ii,j}(m) = mean(params.tempCv2pkAmp{i}{ii,j}{m});

                    % Sort peaks and indexes by chunk windows
                    for t = 1:length(params.chunkSize{i}{ii,j})
                        if t ~= length(params.chunkSize{i}{ii,j})
                            params.chnkIdx{i}{ii,j}{m}(t) = interp1(data2export.idxs{i}{ii,j}{m}/params.sf,data2export.idxs{i}{ii,j}{m}/params.sf,(params.chunkSize{i}{ii,j}(t)+data2export.idxs{i}{ii,j}{m}(1)/params.sf),'nearest');
                            params.chnkIdx{i}{ii,j}{m}(t) = find(params.chnkIdx{i}{ii,j}{m}(t)==(data2export.idxs{i}{ii,j}{m}/params.sf));
                        else
                            params.chnkIdx{i}{ii,j}{m}(t) = find(data2export.idxs{i}{ii,j}{m},1,'last') - 1;
                        end

                        if t == 1
                            data2export.pkAmpChnkAvg{i}{ii,j}{m}(t) = mean(data2export.pks{i}{ii,j}{m}(1:params.chnkIdx{i}{ii,j}{m}(t)));
                            data2export.ipiChnkAvg{i}{ii,j}{m}(t) =  mean(data2export.ipi{i}{ii,j}{m}(1:params.chnkIdx{i}{ii,j}{m}(t)));
                            data2export.cvIpiChnkAvg{i}{ii,j}{m}(t) = std(data2export.ipi{i}{ii,j}{m}(1:params.chnkIdx{i}{ii,j}{m}(t))) / mean(data2export.ipi{i}{ii,j}{m}(1:params.chnkIdx{i}{ii,j}{m}(t)));
                            data2export.cvPkAmpChnkAvg{i}{ii,j}{m}(t) = std(data2export.pks{i}{ii,j}{m}(1:params.chnkIdx{i}{ii,j}{m}(t))) / mean(data2export.pks{i}{ii,j}{m}(1:params.chnkIdx{i}{ii,j}{m}(t)));
                            data2export.freqChnk{i}{ii,j}{m}(t) = numel(data2export.pks{i}{ii,j}{m}(1:params.chnkIdx{i}{ii,j}{m}(t))) / (params.chunkSize{i}{ii,j}(1) - (data2export.idxs{i}{ii,j}{m}(1) / params.sf));
                            for tt = 2:params.chnkIdx{i}{ii,j}{m}(t)
                                params.tempCv2ipiChnk{i}{ii,j}{m,t}(tt-1) = (std([data2export.ipi{i}{ii,j}{m}(tt-1) data2export.ipi{i}{ii,j}{m}(tt)]) / mean([data2export.ipi{i}{ii,j}{m}(tt-1) data2export.ipi{i}{ii,j}{m}(tt)])) * sqrt(2);
                                params.tempCv2pkAmpChnk{i}{ii,j}{m,t}(tt-1) = (std([data2export.pks{i}{ii,j}{m}(tt-1) data2export.pks{i}{ii,j}{m}(tt)]) / mean([data2export.pks{i}{ii,j}{m}(tt-1) data2export.pks{i}{ii,j}{m}(tt)])) * sqrt(2);
                            end
                            data2export.cv2ipiChnk{i}{ii,j}{m}(t) = mean(params.tempCv2ipiChnk{i}{ii,j}{m,t});
                            data2export.cv2pkAmpChnk{i}{ii,j}{m}(t) = mean(params.tempCv2pkAmpChnk{i}{ii,j}{m,t});
                            data2export.inertiaChnk{i}{ii,j}{m}(t,:) = data2export.inertia{i}{ii,j}(m,1:params.chunkSize{i}{ii,j}(t)*params.sf);
                            data2export.wndwRangeChnk{i}{ii,j}{m}(t,:) = data2export.wndwRange{i}{ii,j}(m,1:params.slWndwChnk(t));
                        else
                            data2export.pkAmpChnkAvg{i}{ii,j}{m}(t) = mean(data2export.pks{i}{ii,j}{m}(params.chnkIdx{i}{ii,j}{m}(t-1):params.chnkIdx{i}{ii,j}{m}(t)));
                            data2export.ipiChnkAvg{i}{ii,j}{m}(t) =  mean(data2export.ipi{i}{ii,j}{m}(params.chnkIdx{i}{ii,j}{m}(t-1):params.chnkIdx{i}{ii,j}{m}(t)));
                            data2export.cvIpiChnkAvg{i}{ii,j}{m}(t) = std(data2export.ipi{i}{ii,j}{m}(params.chnkIdx{i}{ii,j}{m}(t-1):params.chnkIdx{i}{ii,j}{m}(t))) / mean(data2export.ipi{i}{ii,j}{m}(params.chnkIdx{i}{ii,j}{m}(t-1):params.chnkIdx{i}{ii,j}{m}(t)));
                            data2export.cvPkAmpChnkAvg{i}{ii,j}{m}(t) = std(data2export.pks{i}{ii,j}{m}(params.chnkIdx{i}{ii,j}{m}(t-1):params.chnkIdx{i}{ii,j}{m}(t))) / mean(data2export.pks{i}{ii,j}{m}(params.chnkIdx{i}{ii,j}{m}(t-1):params.chnkIdx{i}{ii,j}{m}(t)));
                            data2export.freqChnk{i}{ii,j}{m}(t) = numel(data2export.pks{i}{ii,j}{m}(params.chnkIdx{i}{ii,j}{m}(t-1):params.chnkIdx{i}{ii,j}{m}(t))) / params.chunkSize{i}{ii,j}(1);
                            for tt = params.chnkIdx{i}{ii,j}{m}(t-1):params.chnkIdx{i}{ii,j}{m}(t)
                                params.tempCv2chnk{i}{ii,j}{m,t}(tt-1) = (std([data2export.ipi{i}{ii,j}{m}(tt-1) data2export.ipi{i}{ii,j}{m}(tt)]) / mean([data2export.ipi{i}{ii,j}{m}(tt-1) data2export.ipi{i}{ii,j}{m}(tt)])) * sqrt(2);
                                params.tempCv2pkAmpChnk{i}{ii,j}{m,t}(tt-1) = (std([data2export.pks{i}{ii,j}{m}(tt-1) data2export.pks{i}{ii,j}{m}(tt)]) / mean([data2export.pks{i}{ii,j}{m}(tt-1) data2export.pks{i}{ii,j}{m}(tt)])) * sqrt(2);
                            end
                            data2export.cv2ipiChnk{i}{ii,j}{m}(t) = mean(params.tempCv2chnk{i}{ii,j}{m,t});
                            data2export.cv2pkAmpChnk{i}{ii,j}{m}(t) = mean(params.tempCv2pkAmpChnk{i}{ii,j}{m,t});
                            data2export.inertiaChnk{i}{ii,j}{m}(t,params.chunkSize{i}{ii,j}(t-1)*params.sf+1:params.chunkSize{i}{ii,j}(t)*params.sf) = data2export.inertia{i}{ii,j}(m,params.chunkSize{i}{ii,j}(t-1)*params.sf+1:params.chunkSize{i}{ii,j}(t)*params.sf);
                            data2export.wndwRangeChnk{i}{ii,j}{m}(t) = mean(data2export.wndwRange{i}{ii,j}(m,params.slWndwChnk(t-1)+1:params.slWndwChnk(t)));
                        end

                    end

                    % reshape data to align each cycle to peak
                    data2export.pkAlign{i}{ii,j}{m} = NaN(15000,10000);
                    params.pkAlignShft = round((1/data2export.freq{i}{ii,j}(m)*params.sf)/4);
                    for t = 2:length(data2export.idxs{i}{ii,j}{m})
                        if data2export.idxs{i}{ii,j}{m}(t-1) - params.pkAlignShft > 0
                            data2export.pkAlign{i}{ii,j}{m}(1:length(data2export.idxs{i}{ii,j}{m}(t-1):data2export.idxs{i}{ii,j}{m}(t)),t-1) = data{2,i}{ii,j}(m,data2export.idxs{i}{ii,j}{m}(t-1) - params.pkAlignShft:data2export.idxs{i}{ii,j}{m}(t) - params.pkAlignShft);
                        end
                    end
                    data2export.pkAlignAvg{i}{ii,j}(:,m) = nanmean(data2export.pkAlign{i}{ii,j}{m},2);

                    params.counter = [70, 71, 72];
                    figure(params.counter(m)); hold on
                    plot(data2export.pkAlign{i}{ii,j}{m})
                    plot(data2export.pkAlignAvg{i}{ii,j}(:,m),'k','LineWidth',3)
                    xline(params.pkAlignShft,'--',{'Peak'},'LineWidth',1.5)
                    xlabel('Time (ms)')
                    ylabel('Acceleration (g)')
                end
                
            end

            if j ~= 5
                [data2export.amp{i}{ii,j},params.bestCh{i}{ii,j}] = max(data2export.rangeAvg{i}{ii,j});
                data2export.inertiaTotal{i}{ii,j} = data2export.inertia{i}{ii,j}(1,:) + data2export.inertia{i}{ii,j}(2,:) + data2export.inertia{i}{ii,j}(3,:);
            end

        end
        
    end
    
    % Transfer data to excel
    %params.cellFormat{1} = ['E', num2str(1)];
    %writematrix(params.sheetNames{i}, params.fileName,'Sheet','toDoList','Range',params.cellFormat{1});
    %writecell(cellstr(params.endOfMonitor),params.fileName,'Sheet',params.sheetNames{i},'Range',params.cellFormat{1})

    % Plot summary statistics

    % Fit line and add peaks to total inertia

end

cd(groupList(1).folder)

%% Calculate population averages

params.maxNumSxs = max(params.numSxs);
for i = 1:params.numFolders

    for j = 1:length(params.names)

        %
        data2export.freqPop{i,j} = NaN(params.maxNumSxs,3);
        data2export.cvIpiPop{i,j} = NaN(params.maxNumSxs,3);
        data2export.cv2ipiPop{i,j} = NaN(params.maxNumSxs,3);
        data2export.peakAmpPop{i,j} = NaN(params.maxNumSxs,3);
        data2export.rangeAmpPop{i,j} = NaN(params.maxNumSxs,3);
        data2export.cvPkAmpPop{i,j} = NaN(params.maxNumSxs,3);
        data2export.cv2pkAmpPop{i,j} = NaN(params.maxNumSxs,3);
        %data2export.inertiaTotalAmpPop{i,j} = NaN(params.maxNumSxs,3);
        
        for ii = 1:params.numSxs(i)

            for m = 1:3

                %data2export.inertiaTotalRange{i}{ii,j}(t-1) = range(data2export.inertiaTotal{i}{ii,j}(params.wndwTime*(t-1)+1:params.wndwTime*(t)));


                if j ~= 5

                    %data2export.inertiaTotalAmp{i}{ii,j} = mean(data2export.inertiaTotalRange{i}{ii,j},2);

                    data2export.freqPop{i,j}(ii,m) =  data2export.freq{i}{ii,j}(m);
                    data2export.cvIpiPop{i,j}(ii,m) =  data2export.cvIpi{i}{ii,j}(m);
                    data2export.cv2ipiPop{i,j}(ii,m) =  data2export.cv2ipi{i}{ii,j}(m);
                    data2export.pkAmpPop{i,j}(ii,m) =  data2export.pkAmpAvg{i}{ii,j}(m);
                    data2export.rangeAmpPop{i,j}(ii,m) =  data2export.rangeAvg{i}{ii,j}(m);
                    data2export.cvPkAmpPop{i,j}(ii,m) =  data2export.cvPkAmp{i}{ii,j}(m);
                    data2export.cv2pkAmpPop{i,j}(ii,m) =  data2export.cv2pkAmp{i}{ii,j}(m);
                    %data2export.inertiaTotalAmpPop{i,j}(ii) = data2export.inertiaTotalAmp{i}{ii,j};
                end

            end

        end

    end

end

%% 
for i = 1:params.numFolders

    for j = 1:length(params.names)

        for m = 1:3

            data2export.freqPopChnk{i,j}{m} =  NaN(params.maxNumSxs,length(params.chunkSize{1}{1}));
            data2export.cvIpiPopChnk{i,j}{m} =  NaN(params.maxNumSxs,length(params.chunkSize{1}{1}));
            data2export.cv2ipiPopChnk{i,j}{m} =  NaN(params.maxNumSxs,length(params.chunkSize{1}{1}));
            data2export.cvPkAmpPopChnk{i,j}{m} =  NaN(params.maxNumSxs,length(params.chunkSize{1}{1}));
            data2export.cv2PkAmpPopChnk{i,j}{m} =  NaN(params.maxNumSxs,length(params.chunkSize{1}{1}));
            data2export.pkAmpPopChnk{i,j}{m} =  NaN(params.maxNumSxs,length(params.chunkSize{1}{1}));
            data2export.rangeAmpPopChnk{i,j}{m} =  NaN(params.maxNumSxs,length(params.chunkSize{1}{1}));
            %data2export.inertiaTotalAmpPopChnk{i,j} =  NaN(params.maxNumSxs,length(params.chunkSize{1}{1}));

            for ii = 1:params.numSxs(i)

                if j ~= 5

                    for t = 1:length(params.chunkSize{1}{1})

                        data2export.freqPopChnk{i,j}{m}(ii,t) =  data2export.freqChnk{i}{ii,j}{m}(t);
                        data2export.cvIpiPopChnk{i,j}{m}(ii,t) =  data2export.cvIpiChnkAvg{i}{ii,j}{m}(t);
                        data2export.cv2ipiPopChnk{i,j}{m}(ii,t) =  data2export.cv2ipiChnk{i}{ii,j}{m}(t);
                        data2export.pkAmpPopChnk{i,j}{m}(ii,t) =  data2export.pkAmpChnkAvg{i}{ii,j}{m}(t);
                        data2export.cvPkAmpPopChnk{i,j}{m}(ii,t) =  data2export.cvPkAmpChnkAvg{i}{ii,j}{m}(t);
                        data2export.cv2pkAmpPopChnk{i,j}{m}(ii,t) =  data2export.cv2pkAmpChnk{i}{ii,j}{m}(t);
                        data2export.rangeAmpPopChnk{i,j}{m}(ii,t) =  data2export.wndwRangeChnk{i}{ii,j}{m}(t);
                        %data2export.inertiaTotalAmpPopChnk{i,j}(ii,t) =  data2export.inertiaTotalAmpChnk{i}{ii,j}(t);
                    
                        if t == 1
                            %data2export.inertiaTotalAmpChnk{i}{ii,j}(t) = data2export.inertiaTotalRange{i}{ii,j}(1:params.slWndwChnk(t));
                        else
                            %data2export.inertiaTotalAmpChnk{i}{ii,j}(t) = data2export.inertiaTotalRange(params.slWndwChnk(t-1)+1:params.slWndwChnk(t));
                        end
                    end

                end

            end

        end

    end

end

%% 

% data2export.freqPop = cellfun(@(M) subsasgn(M, substruct('()', {M==0}), NaN), data2export.freqPop, 'uniform', 0);
% data2export.cvIpiPop = cellfun(@(M) subsasgn(M, substruct('()', {M==0}), NaN), data2export.cvIpiPop, 'uniform', 0);
% data2export.cv2ipiPop = cellfun(@(M) subsasgn(M, substruct('()', {M==0}), NaN), data2export.cv2ipiPop, 'uniform', 0);
% data2export.pkAmpPop = cellfun(@(M) subsasgn(M, substruct('()', {M==0}), NaN), data2export.pkAmpPop, 'uniform', 0);
% data2export.rangeAmpPop = cellfun(@(M) subsasgn(M, substruct('()', {M==0}), NaN), data2export.rangeAmpPop, 'uniform', 0);

data2export.freqPopAvg = cellfun(@(x) nanmean(x,1),data2export.freqPop,'UniformOutput',false);
data2export.cvIpiPopAvg = cellfun(@(x) nanmean(x,1),data2export.cvIpiPop,'UniformOutput',false);
data2export.cv2ipiPopAvg = cellfun(@(x) nanmean(x,1),data2export.cv2ipiPop,'UniformOutput',false);
data2export.pkAmpPopAvg = cellfun(@(x) nanmean(x,1),data2export.pkAmpPop,'UniformOutput',false);
data2export.cvPkAmpPopAvg = cellfun(@(x) nanmean(x,1),data2export.cvPkAmpPop,'UniformOutput',false);
data2export.cv2pkAmpPopAvg = cellfun(@(x) nanmean(x,1),data2export.cv2pkAmpPop,'UniformOutput',false);
data2export.rangeAmpPopAvg = cellfun(@(x) nanmean(x,1),data2export.rangeAmpPop,'UniformOutput',false);
%data2export.inertiaTotalAmpPopAvg = cellfun(@(x) nanmean(x,1),data2export.inertiaTotalAmpPop,'UniformOutput',false);

for i = 1:params.numFolders
    for j = 1:length(data2export.freqPopChnk)-1
%             data2export.freqPopChnk{i,j} = cellfun(@(M) subsasgn(M, substruct('()', {M==0}), NaN), data2export.freqPopChnk{i,j}, 'uniform', 0);
%             data2export.cvIpiPopChnk{i,j} = cellfun(@(M) subsasgn(M, substruct('()', {M==0}), NaN), data2export.cvIpiPopChnk{i,j}, 'uniform', 0);
%             data2export.cv2ipiPopChnk{i,j} = cellfun(@(M) subsasgn(M, substruct('()', {M==0}), NaN), data2export.cv2ipiPopChnk{i,j}, 'uniform', 0);
%             data2export.pkAmpPopChnk{i,j} = cellfun(@(M) subsasgn(M, substruct('()', {M==0}), NaN), data2export.pkAmpPopChnk{i,j}, 'uniform', 0);
%             data2export.rangeAmpPop{i,j} = cellfun(@(M) subsasgn(M, substruct('()', {M==0}), NaN), data2export.rangeAmpPopChnk{i,j}, 'uniform', 0);
        
        
            data2export.freqPopChnkAvg{i,j} = cellfun(@(x) nanmean(x,1),data2export.freqPopChnk{i,j},'UniformOutput',false);
            data2export.cvIpiPopChnkAvg{i,j} = cellfun(@(x) nanmean(x,1),data2export.cvIpiPopChnk{i,j},'UniformOutput',false);
            data2export.cv2ipiPopChnkAvg{i,j} = cellfun(@(x) nanmean(x,1),data2export.cv2ipiPopChnk{i,j},'UniformOutput',false);
            data2export.pkAmpPopChnkAvg{i,j} = cellfun(@(x) nanmean(x,1),data2export.pkAmpPopChnk{i,j},'UniformOutput',false);
            data2export.rangeAmpPopChnkAvg{i,j} = cellfun(@(x) nanmean(x,1),data2export.rangeAmpPopChnk{i,j},'UniformOutput',false);
            data2export.cvPkAmpPopChnkAvg{i,j} = cellfun(@(x) nanmean(x,1),data2export.cvPkAmpPopChnk{i,j},'UniformOutput',false);
            data2export.cv2pkAmpPopChnkAvg{i,j} = cellfun(@(x) nanmean(x,1),data2export.cv2pkAmpPopChnk{i,j},'UniformOutput',false);
            %data2export.inertialTotalAmpPopChnkAvg{i,j} = cellfun(@(x) nanmean(x,1),data2export.inertiaTotalAmpPopChnk{i,j},'UniformOutput',false);

    end
end

params.measures = {'Freq.','Freq CV','Freq CV2','Amp.','Sliding Amp','Amp CV','Amp CV2','Inertia Amp'};
plotChs(data2export.freqPopAvg,data2export.freqPop,params,groupList,1,1);
plotChs(data2export.cvIpiPopAvg,data2export.cvIpiPop,params,groupList,2,1);
plotChs(data2export.cv2ipiPopAvg,data2export.cv2ipiPop,params,groupList,3,1);
plotChs(data2export.pkAmpPopAvg,data2export.pkAmpPop,params,groupList,4,1);
plotChs(data2export.rangeAmpPopAvg,data2export.rangeAmpPop,params,groupList,5,1);
plotChs(data2export.cvPkAmpPopAvg,data2export.cvPkAmpPop,params,groupList,2,1);
plotChs(data2export.cv2pkAmpPopAvg,data2export.cv2pkAmpPop,params,groupList,3,1);

plotChs(data2export.freqPopChnkAvg,data2export.freqPopChnk,params,groupList,1,4);
plotChs(data2export.cvIpiPopChnkAvg,data2export.cvIpiPopChnk,params,groupList,2,4);
plotChs(data2export.cv2ipiPopChnkAvg,data2export.cv2ipiPopChnk,params,groupList,3,4);
plotChs(data2export.pkAmpPopChnkAvg,data2export.pkAmpPopChnk,params,groupList,4,4);
plotChs(data2export.rangeAmpPopChnkAvg,data2export.rangeAmpPopChnk,params,groupList,5,4);

params.tempName = datetime('now','Format','dd/MM/yyyy');
params.tempName = datestr(params.tempName,'ddmmmyy');
save(['accel' params.tempName '.mat'])
