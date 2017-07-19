% ENNTest: The implementation of ENN method using N-Folds cross validation approach.
% Input: Data, a matrix for training data
%        Label, a vector for training class label
%        K, the number of nearest neighbors are considered
%        NFold, number of folds for cross validation
% Output: ACC_KNN, return the performance of KNN method
%         ACC_ENN, return the performance of ENN method
% Date: 11/08/2014
% By Bo Tang (btang@ele.uri.edu) and Haibo He (he@ele.uri.edu)
% For any questions and/or comments for this code/paper, please feel free
% to contact Prof. Haibo He, Electrical Engineering, University of Rhode Island,
% Email: he@ele.uri.edu
% Web:  http://www.ele.uri.edu/faculty/he/
function [ACC_KNN, ACC_ENN] = ENNTest(Data, Label, K, NFold)
if (size(Data, 1) ~= length(Label))
    error('error: please check the size of data and corresponding class label.');
end
if (NFold <= 0)
    error('error: invalid parameter of NFolds.'); 
end

NumClass = length(unique(Label));
NumSample = size(Data, 1);
if (NFold >= NumSample)
    error('error: invalid parameter of NFolds which is larger than the size of data set.');
end

if (K >= ceil((NFold-1) * NumSample/NFold))
    error('error: K is too large.');
end

%% N-Folds cross validation approach
ACC_KNN = zeros(NFold, 1);
ACC_ENN = zeros(NFold, 1);
disp('N-Fold cross validation: .............')
for nf = 1 : NFold
    disp(['   ' num2str(nf) '-th fold .....' ])
    startIX = ceil((nf-1) * NumSample/NFold) + 1; endIX = ceil(nf * NumSample/NFold);
    if (endIX > NumSample)
        endIX = NumSample;
    end
    TestingData = Data(startIX:endIX ,:); TestingLabel = Label(startIX:endIX ,:);
    TrainingData = Data;  TrainingData(startIX:endIX,:) = [];
    TrainingLabel = Label; TrainingLabel(startIX:endIX,:) = [];
    
    %% Preprocessing stage: Build weighted KNN map in training data set
    [weightedKNNDistance, weightedKNNLabel, TSOri] = buildWeightedKNNMap(TrainingData, TrainingLabel, NumClass, K);
    
    NumTrainingEachClass = zeros(NumClass, 1);
    NumTestingEachClass = zeros(NumClass, 1);
    for i = 1 : NumClass
        NumTrainingEachClass(i) = length(find(TrainingLabel == i));
        NumTestingEachClass(i) = length(find(TestingLabel == i));
    end
    
    numTrainingData = size(TrainingData, 1);
    
    countKNN = zeros(NumClass,1);
    countENN = zeros(NumClass,1);
    for i = 1:NumClass   %% for the i-th testing class
        [classIdx,v] = find(TestingLabel == i);
        numTestingClass = length(classIdx);
        testingDataInd = TestingData(classIdx,:);  %% for the i-th class data
        for k = 1: numTestingClass
            testingData1 = testingDataInd(k, :);
            %% calculate the distance to all training data
            disNorm2 = sqrt(sum((ones(numTrainingData, 1) * testingData1 -  TrainingData).^2,2));
            [normSort, sortIX] = sort(disNorm2);
            
            classNNTest = TrainingLabel(sortIX(1:K),1);
            
            %% search for the k nearest points from the new test
            hitNumKNN = zeros(NumClass,1);
            for m = 1:NumClass
                hitNumKNN(m,1) = length(find(classNNTest == m));
            end
            [C0,V] = find(max(hitNumKNN) == hitNumKNN);
            if (C0 == i)
                countKNN(i,1) = countKNN(i,1) + 1;
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
            [C1,V] = find(max(TSENN) == TSENN);
            if (C1 == i)
                countENN(i,1) = countENN(i,1) + 1;
            end    
        end
    end
    ACC_KNN(nf, 1) = sum(countKNN) / sum(NumTestingEachClass);
    ACC_ENN(nf, 1) = sum(countENN) / sum(NumTestingEachClass);
end
ACC_KNN = mean(ACC_KNN); ACC_ENN = mean(ACC_ENN);
end
