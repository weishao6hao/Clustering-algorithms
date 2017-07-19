% BuildWeightedKNNMap: build weighted KNN map in which the K nearest
% neighbors for each training data are found.
% Input: trainingData: a N times D matrix storing training data.
%        trainingLabel: a N times 1 vector storing corresponding  training class label.
%        NumClass: the number of class.
%        K: the considered number of nearest neighbors.
% Return: 
%        weightedKNNDistance: a N times max(KNumList) matrix storing the
%                   weights (distance) to max(KNumList) nearest neighbors for each training data.
%        weightedKNNLabel: a N times max(KNumList) matrix storing the class
%                   label of their nearest neighbors.
%        TSOri: a vector storing the class-wise satistics.
function [weightedKNNDistance, weightedKNNLabel, TSOri] = buildWeightedKNNMap(trainingData, trainingLabel, NumClass, K)
numData = size(trainingData, 1);

%% compute the mutual distance for training Data
Inf = 1e+10;
TSOri = zeros(NumClass,length(K));
maxK = max(K);
weightedKNNDistance = zeros(numData, maxK);  
weightedKNNLabel = zeros(numData, maxK);
for i = 1 : numData
    trainingData1 = trainingData(i,:);
    disNorm2 = sqrt(sum((ones(numData, 1) * trainingData1 -  trainingData).^2,2));
    disNorm2(i,1) = max(disNorm2) + 1;  %% exclude itself
    [normSort sortIX] = sort(disNorm2);
    weightedKNNDistance(i, :) = normSort(1:maxK);
    weightedKNNLabel(i, :) = trainingLabel(sortIX(1:maxK));
end

for i = 1 : NumClass
    count = zeros(length(K),1);
    [c,v] = find(trainingLabel == i);
    numTraining = length(c);
    labelTmp = weightedKNNLabel(c,:);
    for m = 1:length(K)
        KNum = K(m);
        %% compute the original TS 
        tmpKNNClass = labelTmp(:, 1:KNum);
        count(m,1) = count(m,1) + length(find(tmpKNNClass == i)) / (numTraining*KNum);
    end
    TSOri(i,:) = count.';    
end

end