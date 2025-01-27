%   Sliding window function
%
%   Written by Alex Fanning on 2/22/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fitOfData = fitLine(parameters,dataIn)

N = length(dataIn);

for i = 1:N
    dom = i-parameters(1).wndwSize/2:i+parameters(1).wndwSize/2;
    dom = dom(find(dom>=1 & dom<=N));
    pc = dataIn(dom);
    mi = prctile(pc,75);
    ma = prctile(pc,85);
    id = find(pc>mi & pc<ma);
    fitOfData(i) = median(pc(id));
end