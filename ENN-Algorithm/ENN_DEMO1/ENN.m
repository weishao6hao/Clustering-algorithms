% ENNTest: The implementation of ENN method.
% Input: TrainingData, Nr-by-D matrix for training data
%        TrainingLabel, Nr-by-1 vector for training class label
%        TestingData,Nt-by-D matrix for testing data
%        K, the number of nearest neighbors are considered
% Output: PredictionLabel, the predicted class label using ENN method
% Date: 11/08/2014
% By Bo Tang (btang@ele.uri.edu) and Haibo He (he@ele.uri.edu)
% For any questions and/or comments for this code/paper, please feel free
% to contact Prof. Haibo He, Electrical Engineering, University of Rhode Island,
% Email: he@ele.uri.edu
% Web:  http://www.ele.uri.edu/faculty/he/
function [PredictionLabel] = ENN(TrainingData, TrainingLabel, TestingData, K)

NumClass = length(unique(TrainingLabel));

%% Preprocessing stage: Build weighted KNN map in training data set
[weightedKNNDistance, weightedKNNLabel, TSOri] = buildWeightedKNNMap(TrainingData, TrainingLabel, NumClass, K);

NumTrainingEachClass = zeros(NumClass, 1);
for i = 1 : NumClass
    NumTrainingEachClass(i) = length(find(TrainingLabel == i));
end

numTrainingData = size(TrainingData, 1);

PredictionLabel = zeros(size(TestingData, 1), 1);
for i = 1:size(TestingData, 1)
    testingData1 = TestingData(i, :);
    
    disNorm2 = sqrt(sum((ones(numTrainingData, 1) * testingData1 -  TrainingData).^2,2));
    [normSort, sortIX] = sort(disNorm2);
    
    classNNTest = TrainingLabel(sortIX(1:K),1);
    
    %% search for the k nearest points from the new test
    hitNumKNN = zeros(NumClass,1);
    for m = 1:NumClass
        hitNumKNN(m,1) = length(find(classNNTest == m));
    end
    
    %% search for the k nearest points from the new test
    TSENN = zeros(NumClass,1);
    numTrainingNN = zeros(NumClass,1);
    numSameTrainingNN = zeros(NumClass,1);
    for ii = 1 : NumClass %% count the changes for class ii (i in the paper).
        [c,v] = find(TrainingLabel == ii);
        testingMuDis = disNorm2(c,1);
        trainingMuDis = weightedKNNDistance(c, K);  %% Only examine the KNum-th nearest neighbor
        trainingMuClass = weightedKNNLabel(c, K);
        difDis = testingMuDis - trainingMuDis;
        
        %                 testingMuDis = testingMuClassSort(i,k).mu(ii,1:NumTrainingEachClass(ii));   %% get the mutual distance between this testing data and all the training data
        %                 trainingMuDis = reshape(trainingMuDisSort(ii,1:NumTrainingEachClass(ii),KNum),NumTrainingEachClass(ii),1);  %% get the distances between the training data and their k-th nearest neighbor.
        %                 trainingMuClass =  reshape(trainingMuClassSort(ii, 1:NumTrainingEachClass(ii), KNum),NumTrainingEachClass(ii),1);  %% get the class label of the k-th nearest neighbor for all the training data
        %                 difDis = reshape(testingMuDis,NumTrainingEachClass(ii),1) - trainingMuDis;  %% check whether the testing data will be in training data's K nearest neighbors.
        [C,V] = find(difDis <=0);
        numTrainingNN(ii,1) = length(C);
        if (numTrainingNN(ii,1) > 0)
            numSameTrainingNN(ii,1) = length(find(trainingMuClass(C) == ii));
        end
    end
    
    for jj = 1 : NumClass   %% assume that the testing data to be class jj (j in the paper).
        deltaNumSame = numTrainingNN(jj,1) - numSameTrainingNN(jj,1);   %% count the training data if the k-th nearest neighbor doesn't belong to the same class because the number of nearest neighbors from the same class will increase by assigning class jj to the testing data
        difTmp = numSameTrainingNN ./ (NumTrainingEachClass*K); %% count the training data if k-th nearest neighbor belongs to the same class because the number of nearest neighbors from the same class will decrease by assigning class jj to the testing data
        deltaNumDif = sum(difTmp,1) - numSameTrainingNN(jj)/(NumTrainingEachClass(jj)*K);
        TSENN(jj,1) = (deltaNumSame + hitNumKNN(jj,1) - TSOri(jj,1) * K) / ((NumTrainingEachClass(jj)+1)*K) - deltaNumDif;
    end
    [maxv, maxix] = max(TSENN);
    if (length(maxix) > 1)
        PredictionLabel(i, 1) = maxix(ceil(rand(1) * length(maxix)));
    else
        PredictionLabel(i, 1) = maxix;
    end
    
end

end
