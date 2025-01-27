%
%
%
%
clear; close all

data = [{ctrl} {atx} {et}];
for i = 1:length(data)
    dataAvg{i} = mean(data{i},1);
    dataStdEr{i} = std(data{i},0,1)/sqrt(size(data{i},1));
    dataAvg{i}(1) = NaN;
end

xLength = 30;
yLength = 0.06;
init = 0;

filColor = {'r','g','k'};
figure(); hold on
lnCol = {[1 0.65 0 1],[.05 .73 .41],[0.28,0.04,0.33 1]};
for ii = 1:length(data)
    bound1 = dataAvg{ii} + dataStdEr{ii};
    bound2 = dataAvg{ii} - dataStdEr{ii};
    x = (1:length(dataAvg{i}));
    plot(dataAvg{ii}, 'LineWidth',3,'Color',lnCol{ii})
    shade(x(2:end),bound1(2:end),x(2:end),bound2(2:end),'FillType',[1 2],'FillColor',filColor{ii},'FillAlpha',0.1)
end
xlim([init xLength])
ylim([init yLength])
xlabel('Frequency')
ylabel('PSD (norm)')
yticks([0 yLength])
set(gca,'TickDir','out');
set(gca,"FontSize",14)
legend('','','','Ctrl','','','','Atx','','','','PD')
