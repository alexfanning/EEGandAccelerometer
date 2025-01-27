%
%   Calculate amplitude of oscillations in a sliding window
%
%   Written by Alex Fanning on 5/1/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parameters,dataOut] = ampSlWndw(dataIn,parameters,dataOut,count)

slWndwLng = floor(parameters.chunkSize{count(1)}(end) / (1 / dataOut.freq{count(1)}(count(2))));
wndwTime = floor((1 / dataOut.freq{1,count(1)}(count(2))) * parameters.sf);
slWndwStart = ceil((dataOut.idxs{count(1),count(2)}(1)/parameters.sf) / (1 / dataOut.freq{count(1)}(count(2))));
for t = slWndwStart:slWndwLng-1
    dataOut.wndwRange{count(1),count(2)}(t) = abs(range(dataIn{2,count(1)}(count(2),wndwTime*(t):wndwTime*(t+1))));
    %dataOut.wndwRange{count(1),count(2)}(t-slWndwStart) = abs(range(dataIn{2,count(1)}(count(2),wndwTime*(t-1)+1:wndwTime*(t))));
end

dataOut.rangeAvg{count(1)}(count(2)) = mean(dataOut.wndwRange{count(1),count(2)}(:));

slWndwSpan = (slWndwLng - slWndwStart);
slWndwBlkLng = floor(slWndwSpan / 3);
parameters.slWndwChnk{count(1)} = [slWndwBlkLng, slWndwBlkLng*2, slWndwBlkLng*3];