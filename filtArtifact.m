%   
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dataOut = filtArtifact(dataOut,parameters,counters)

moveOn = 'n';
artThresh = 1;
diffOrder = 2;
j = counters(1);
m = counters(2);

while moveOn == 'n'
    dataIn = dataOut{3,j}(m,:);
    diffData = diff(dataIn,diffOrder);
    diffAvg = mean(diffData);
    avgSubtract = diffData - diffAvg;
    diffStd = std(avgSubtract)*artThresh;
    
    threshPts{1} = find(avgSubtract > diffAvg + diffStd);
    threshPts{2} = find(avgSubtract < diffAvg - diffStd);
    threshPts{3} = sort(cat(2,threshPts{1},threshPts{2}));
    
    consecPts{1} = diff(threshPts{3});
    %histogram(consecPts{1},1000)
    
    for t = 1:length(consecPts{1})
        if consecPts{1}(t) < parameters.thresh
            consecPts{1}(t) = 1;
        end
    end
    
    consecPts{2} = find(consecPts{1}~=1);
    
    for t = 1:length(consecPts{2})
        if t == 1
            x = threshPts{3}(t):threshPts{3}(consecPts{2}(t));
            if threshPts{3}(t)-1 > 0
                y(1) = dataIn(threshPts{3}(t)-1);
            else
                y(1) = dataIn(threshPts{3}(t));
            end
    
            if threshPts{3}(consecPts{2}(t)+1) <= length(dataOut{2,j})
                y(2) = dataIn(threshPts{3}(consecPts{2}(t)+1));
            else
                 y(2) = length(dataOut{2,j});
            end
    
            y2{t} = linspace(y(1),y(2),length(x));
    
            dataIn(x) = y2{t};
        else
            x = threshPts{3}(consecPts{2}(t-1)+1):threshPts{3}(consecPts{2}(t));
            y(1) = dataIn(threshPts{3}(consecPts{2}(t-1)+1)-1);
            y(2) = dataIn(threshPts{3}(consecPts{2}(t))+1);
            y2{t} = linspace(y(1),y(2),length(x));
            dataIn(x) = y2{t};
        end
    end
    dataOut{2,j}(m,:) = dataIn;
    
    figure(); hold on
    plot(dataOut{3,j}(m,:),'LineWidth',1)
    plot(dataIn(:),'LineWidth',1)
    plot(diffData,'LineWidth',1)
    yline(diffAvg,'LineWidth',1)
    yline(diffAvg+diffStd,'k','LineWidth',1)
    yline(diffAvg-diffStd','r','LineWidth',1)
    legend('Raw','Filtered','Diff','Rec. avg','Avg + std', 'Avg - std')
    xlabel('Sample number')
    ylabel('Acceleration (g)')
    set(gca,'FontSize',16)

    moveOn = input(parameters.prompt2,'s');
    if moveOn == 'n'
        params = inputdlg(parameters.prompt,parameters.dlgtitle,[1 40],parameters.base);
        diffOrder = str2double(params{2});
        artThresh = str2double(params{1});
    else
        close all
    end
end
