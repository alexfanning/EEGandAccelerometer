function [pValue] = unpairedTtest(sample1,sample2)



% Set the number of randomizations (iterations)
numIterations = 1000;

% Calculate the observed test statistic
[h,p,ci,observedTStat] = ttest2(sample1, sample2);
tStat = abs(observedTStat.tstat);

% Combine the samples
combinedSamples = [sample1, sample2];
numSubjects = length(combinedSamples);

% Initialize an array to store the randomization test statistics
randTestStats = zeros(numIterations, 1);

% Perform Monte Carlo randomizations
for i = 1:numIterations
    % Generate a random permutation index
    permIndices = randperm(numSubjects);
    
    % Permute the combined samples
    permutedSamples = combinedSamples(permIndices);
    
    % Split the permuted samples into two groups
    permSample1 = permutedSamples(1:length(sample1));
    permSample2 = permutedSamples(length(sample1)+1:end);
    
    % Calculate the test statistic for the permuted samples
    [h1,p1,ci1,permTStat] = ttest2(permSample1, permSample2);
    permTStat = abs(permTStat.tstat);
    
    % Store the test statistic
    randTestStats(i) = permTStat;
end

% Calculate the p-value
pValue = (sum(randTestStats >= tStat) + sum(randTestStats <= -tStat) + 1) / (numIterations + 1);

% Display the results
% disp(['Observed t-statistic: ', num2str(observedTStat)]);
% disp(['P-value: ', num2str(pValue)]);