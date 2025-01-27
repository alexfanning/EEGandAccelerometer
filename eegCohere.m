function dataOut = eegCohere(dataIn,dataOut,parameters,timeFrame,sfNum)

n = parameters(sfNum).sf * timeFrame; 
n1s = parameters(sfNum).sf;          
nCoh = n1s * parameters(1).cohSec; 
shiftedframe = n1s * parameters(1).timeshift;
endNum = fix(length(dataIn{2,4}) / n1s) - timeFrame - ceil(parameters(1).timeshift);

for m = 1:3
    vec = {dataIn{2,4}(m,:)' dataIn{2,2}(2,:)'};
    %vec = {dataIn{2,3}(2,:)' dataIn{2,2}(2,:)'};
    for s = 1:endNum
        y1 = vec{1}((1+(s-1)*n1s):(n+(s-1)*n1s));  % Each iteration grabs 20s of data and the starting point is shifted 1 s to the right of the prior starting point
        y2 = vec{2}((1+(s-1)*n1s+shiftedframe):(n+(s-1)*n1s+shiftedframe));
        
        [Cy1y2,FreqCoh] = mscohere(y1,y2,hanning(nCoh),1/2*nCoh,nCoh,parameters(sfNum).sf); % (y1,y2 coherece, with hamming window, take NCoh points, half points overlap, number of fft:NCoh, sampling frequency:fs
        coh(:,s) = Cy1y2;
    end

    dataOut.mCoh(:,m) = mean(coh,2);
    
    [dataOut.maxCoh(1,m),dataOut.maxCohIdx(1,m)] = max(dataOut.mCoh(parameters(1).extractWndw(1):parameters(1).extractWndw(2),m));
    dataOut.maxCohIdx(1,m) = dataOut.maxCohIdx(1,m) + parameters.extractWndw(1) - 1;
end