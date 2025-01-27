%
%   Find peaks
%
%   Written by Alex Fanning on 4/30/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dataVar,dataOut,parameters] = findPks(dataVar,parameters,dataOut,count)

parameters.continue = 'n';
minPeakDist =  10;
if count(1) == 1
    wndwSize =  250;
    minPeakProm = 0.00002;
    minPeakDist =  200;
elseif count(1) == 2
    wndwSize =  10;
    minPeakProm = 0.000001;
    minPeakDist =  25;
elseif count(1) == 3
    wndwSize =  500;
    minPeakProm = 0.000002;
    minPeakDist =  800;
elseif count(1) == 4
    wndwSize =  70;
    minPeakProm = 0.00002;
    minPeakDist =  80;
else
    parameters.continue = 'y';
end

% Smooth data
dataVar{6,count(1)}(count(2),:) = smooth(dataVar{2,count(1)}(count(2),:),wndwSize);
% dataVar{6,count(1)}(count(2),:) = fitLine(parameters,dataVar{2,count(1)}(count(2),:));

while parameters.continue == 'n'

    [dataOut.pks{count(1),count(2)},dataOut.idxs{count(1),count(2)},dataOut.w{count(1),count(2)},dataOut.peaks{count(1),count(2)}] = findpeaks(dataVar{6,count(1)}(count(2),:),1,'MinPeakProminence',minPeakProm,'MinPeakDistance',minPeakDist);

    parameters.titleNamesCh = {'Ch X','Ch Y','Ch Z'};
    figure('Renderer', 'painters', 'Position', [10 500 2000 600]); hold on
    plot(dataVar{2,count(1)}(count(2),:))
    plot(dataVar{6,count(1)}(count(2),:),'LineWidth',1.5)
    scatter(dataOut.idxs{count(1),count(2)},dataOut.pks{count(1),count(2)},200,'k','LineWidth',2)
    xlabel('Sample number')
    ylabel('Acceleration (g)')
    title([parameters.sheetNames{count(1)} parameters.titleNamesCh{count(2)} ': Raw data'])
    set(gca,'FontSize',16)

    prompt = 'Peaks and smoothing accurate? (y/n): ';
    parameters.continue = input(prompt,"s");

    if parameters.continue == 'n'
        prompt2 = {'minPeakProm: ','minInterPkInterval: ','wndwSize: '};
        dlgtitle = 'Find peaks parameters';
        if count(1) == 1
            default = {'0.00002','200','250'};
        elseif count(1) == 2
            default = {'0.000001','25','10'};
        elseif count(1) == 3
            default = {'0.000002','800','500'};
        elseif count(1) == 4
            default = {'0.00002','80','70'};
        end

        tempParams = inputdlg(prompt2,dlgtitle,1,default);
        minPeakProm = str2double(tempParams{1});
        minPeakDist =  str2double(tempParams{2});
        wndwSize =  str2double(tempParams{3});

        % Smooth data
        dataVar{6,count(1)}(count(2),:) = smooth(dataVar{2,count(1)}(count(2),:),wndwSize);
        % dataVar{6,count(1)}(count(2),:) = fitLine(parameters,dataVar{2,count(1)}(count(2),:));
    
        % Get user's selection
        zoom out
        tempDataRemove = input('Remove data? (y/n): ','s');
        if tempDataRemove == 'y'
            rect = getrect(gca);
            rect = round(rect);
            dataVar{6,count(1)}(count(2),1:rect(1)+rect(3)) = dataVar{6,count(1)}(count(2),rect(3));
        end
        close all
    end
end
