clear all
close all
clc

%% The Implementation of Extended Nearest Neighbor (ENN)
% Date: 11/08/2014
% By Bo Tang (btang@ele.uri.edu) and Haibo He (he@ele.uri.edu)
% For any questions and/or comments for this code/paper, please feel free
% to contact Prof. Haibo He, Electrical Engineering, University of Rhode Island,
% Email: he@ele.uri.edu
% Web:  http://www.ele.uri.edu/faculty/he/

%% choose one data set
datasetName = 'wine';
% datasetName = 'knowledge';
% datasetName = 'Vertebral';
% datasetName = 'vowel';
% datasetName = 'breasttissue';
% datasetName = 'breast-cancer';

%% load the data set
dataPath = 'dataori/';
matFile = [dataPath datasetName '.mat'];
load(matFile)

ENN = {};
ENN.data = Data;
ENN.label = Class;
ENN.numFeatures = size(Data,2);
ENN.numClass = length(unique(Class));

%%
% KAll = 5;
ENN.K = 3;   %% The number of nearest neighbors considered in ENN and KNN methods.
ENN.NFold = 10; %% The number of experiments runing in Monte-Carlo mode.

%% initialize the parameters
initPara

if (ENN.normalizeFlag)
    disp('normalize the data ......')
    ENN.data = normalizeData(ENN.data);
    disp('normalize the data ...... Done!')
end

if (ENN.frFlag)
    disp('feature reduction using PCA ......')
    [ENN.data, ENN.frPC, ENN.frVar] = minPCA(ENN.data);
    disp('feature reduction using PCA ...... Done!')
    ENN.numFeatures = size(ENN.data, 2);
end

%% Evaluate classification performance of ENN method using N-Folds Cross Validation
[ACC_KNN, ACC_ENN] = ENNTest(ENN.data, ENN.label, ENN.K, ENN.NFold);

disp(['Averaged accuracy using KNN:' num2str(ACC_KNN)])
disp(['Averaged accuracy using ENN:' num2str(ACC_ENN)])
