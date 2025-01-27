
function parameters = plotChs(dataInAvg,dataInSx,parameters,numTasks,mesure,type)

if mesure == 1
    parameters.ylim = [0 7];
elseif mesure == 2
    parameters.ylim = [0 2];
elseif mesure == 3
    parameters.ylim = [0 1];
elseif mesure == 4
    parameters.ylim = [0 0.0003];
elseif mesure == 5
    parameters.ylim = [0 0.0003];
end

if type == 1
    for j = 1:length(dataInAvg)
        if j ~= 5
            figure()
            for m = 1:length(dataInAvg{1})
                subplot(1,3,m); hold on
                h = bar([dataInAvg{1,j}(m) dataInAvg{2,j}(m)],'FaceColor','flat');
                h.CData(1,:) = [0.9290 0.6940 0.1250];
                h.CData(2,:) = [0.25, 0.50, 0.75];
                xlim([0.25 2.75])
                ylim(parameters.ylim)
                xticks(1:length(numTasks))
                xticklabels({numTasks(1).name; numTasks(2).name})
                xtickangle(45)
                for i = 1:size(dataInAvg,1)
                    parameters.x = repelem(i,size(dataInSx{i,j},1));
                    scatter(parameters.x,dataInSx{i,j}(:,m),50,'k','LineWidth',3)
                end
                if m == 1
                    ylabel(parameters.measures{mesure},'FontSize',14)
                end
                title(parameters.titleNamesCh{m},'FontSize',14)
                sgtitle(parameters.titleNames{j},'FontSize',16)
                set(gca,'FontSize',14)
            end
        end
    end
elseif type == 2
    c = {[0.9290 0.6940 0.1250],[0.25, 0.50, 0.75]};
    for j = 1:length(dataInAvg)
        if j ~= 5
            figure()
            for m = 1:length(dataInAvg{1})
                subplot(1,length(dataInAvg{1}),m); hold on
                plot(dataInAvg{1,j}{m},'Color',c{1},'LineWidth',3);
                plot(dataInAvg{2,j}{m},'Color',c{2},'LineWidth',3);
                for t = 1:size(dataInSx{1,j}{1},1)
                    plot(dataInSx{1,j}{m}(t,:),'Color',c{1},'LineWidth',0.5)
                    plot(dataInSx{2,j}{m}(t,:),'Color',c{2},'LineWidth',0.5)
                end
                xlim([0.5 3.5])
                ylim([parameters.ylim])
                if m == 1
                    ylabel(parameters.measures{mesure},'FontSize',14)
                elseif m == 2
                    xlabel('Block number')
                elseif m == 3
                    legend('Ataxia','Ctrl')
                end
                title(parameters.titleNamesCh{m},'FontSize',14)
                sgtitle(parameters.titleNames{j},'FontSize',16)
                set(gca,'FontSize',14)
            end
        end
    end
elseif type == 3
    for j = 2%1:length(dataInAvg)
        if j ~= 5
            figure()
            for m = 1:length(dataInAvg{1,2})
                subplot(1,3,m); hold on
                h = bar([dataInAvg{1,j}(m) dataInAvg{2,j}(m) dataInAvg{3,j}(m) dataInAvg{4,j}(m) dataInAvg{5,j}(m)],'FaceColor','flat');
                h.CData(1,:) = [0.9290 0.6940 0.1250];
                h.CData(2,:) = [0.25, 0.50, 0.75];
                xlim([0.25 5.75])
                ylim(parameters.ylim)
                xticks(1:length(numTasks))
                xticklabels({numTasks(1).name; numTasks(2).name; numTasks(3).name; numTasks(4).name; numTasks(5).name;})
                xtickangle(45)
                for i = 1:size(dataInAvg,1)
                    parameters.x = repelem(i,size(dataInSx{i,j},1));
                    scatter(parameters.x,dataInSx{i,j}(:,m),50,'k','LineWidth',3)
                end
                if m == 1
                    ylabel(parameters.measures{mesure},'FontSize',14)
                end
                title(parameters.titleNamesCh{m},'FontSize',14)
                sgtitle(parameters.titleNames{j},'FontSize',16)
                set(gca,'FontSize',14)
            end
        end
    end
elseif type == 4
    c = cool(5);
    for j = 2%1:length(dataInAvg)
        if j ~= 5
            figure()
            for m = 1:length(dataInAvg{1,j})
                subplot(1,length(dataInAvg{1,j}),m); hold on
                for i = 1:length(dataInAvg)
                    plot(dataInAvg{i,j}{m},'Color',c(i,:),'LineWidth',3);
%                     for t = 1:size(dataInSx{1,j}{1},1)
%                         plot(dataInSx{i,j}{m}(t,:),'Color','k','LineWidth',0.5)
%                     end
                end
              
                xlim([0.5 3.5])
                ylim([parameters.ylim])
                if m == 1
                    ylabel(parameters.measures{mesure},'FontSize',14)
                elseif m == 2
                    xlabel('Block number')
                elseif m == 3
                    legend('Baseline','Sham','10hz','15hz','30hz')
                end
                title(parameters.titleNamesCh{m},'FontSize',14)
                sgtitle(parameters.titleNames{j},'FontSize',16)
                set(gca,'FontSize',14)
            end
        end
    end
end
