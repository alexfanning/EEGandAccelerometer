

close all

data2plot = {ctrl, ataxia};
for i = 1:length(data2plot)

    for m = 1:length(data2plot{i})

        popAvg{i}(m,:) = mean(data2plot{i}{m},1,'omitnan');
        popSe{i}(m,:) = std(data2plot{i}{m},0,1,'omitnan') / sqrt(size(data2plot{i}{m},1));

    end
end

colr = {[0 0 0];[0.93,0.69,0.13]; [0.00,0.45,0.74]};
chNames = {'PO10';'PO12';'P10';'P12';'TP12'};
legNames = {'Ctrl','Ataxia'};
for m = 1:length(data2plot{i})
    figure(m); hold on
    for i = 1:length(data2plot)
        bound1 = popAvg{i}(m,:) + popSe{i}(m,:);
        bound2 = popAvg{i}(m,:) - popSe{i}(m,:);
        x = (1:length(popAvg{i}(m,:)));
        plot(popAvg{i}(m,:),'LineWidth',1.5,'Color',colr{i})
        shade(x,bound1,x,bound2,'FillType',[1 2],'FillColor',colr{i},'FillAlpha',0.3)
    end
    ax = gca;
    ax.Box = 'off';
    ax.LineWidth = 1.5;
    ax.FontSize = 10;
    ax.FontName = 'Times New Roman';
    ax.TickDir = 'out';
    ax.TickLength = [0.02, 0.02];
    ax.XAxis.MinorTick = 'off';
    ax.YAxis.MinorTick = 'off';
    ylabel('Norm PSD','FontSize',20)
    legend(legNames{1},'','','',legNames{2},'FontSize',16,'Location','best')
    xlabel('Frequency (Hz)','FontSize',20)
    set(gca,'FontSize',18)
    title(['Ch ' chNames{m}],'FontSize',20)
    ylim([0 0.1])
end

%% 

grpNames = {'Ctrl';'Ataxia'};
for m = 1:length(data2plot{i})
    for i = 1:length(data2plot)
        figure(); hold on
        colr = copper(size(data2plot{i}{m},1));
        x = 1:size(data2plot{i}{m},2);
        for w = 1:size(data2plot{i}{m},1)
            plot(data2plot{i}{m}(w,:)','LineWidth',1.5,'Color',colr(w,:))
%             text(10,data2plot{i}{m}(w,10),num2str(w),'FontSize',14,'Color',colr(w,:))
        end
        title(['Ch ' chNames{m} ' ' grpNames{i}],'FontSize',20)
        ax = gca;
        ax.Box = 'off';
        ax.LineWidth = 1.5;
        ax.FontSize = 10;
        ax.FontName = 'Times New Roman';
        ax.TickDir = 'out';
        ax.TickLength = [0.02, 0.02];
        ax.XAxis.MinorTick = 'off';
        ax.YAxis.MinorTick = 'off';
        ylabel('Norm PSD','FontSize',20)
        xlabel('Frequency (Hz)','FontSize',20)
        set(gca,'FontSize',18)
        ylim([0 0.1])
    end
end

%% 

for i = 1:length(data2plot)
    for m = 1:length(data2plot{i})
        figure(); hold on
        bx2 = boxplot(data2plot{i}{m},'Notch','on','whisker',7,'Colors','k');
        for k = 1:size(data2plot{i}{m},2)
            x = repmat(k,size(data2plot{i}{m},1),1);
            scatter(k+(rand(size(data2plot{i}{m},1),1)-0.5)/10,data2plot{i}{m}(:,k),50,'b','LineWidth',1)
        end
        set(bx2,'LineWidth',1.5)

        ax = gca;
        ax.Box = 'off';
        ax.LineWidth = 1.5;
        ax.FontSize = 10;
        ax.FontName = 'Times New Roman';
        ax.TickDir = 'out';
        ax.TickLength = [0.02, 0.02];
        ax.XAxis.MinorTick = 'off';
        ax.YAxis.MinorTick = 'off';
        title(['Ch ' chNames{m} ' ' grpNames{i}],'FontSize',20)
        ylim([0 0.3])
        ylabel('Norm. PSD','FontSize',16)
        xlabel('Frequency (Hz)','FontSize',20)
        set(gca,'FontSize',14)
    end
end

save('tappingCompiled.mat')
