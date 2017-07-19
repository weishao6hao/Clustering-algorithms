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

%% This implementation is based on ENN.V1 as describled in our following paper.
% B. Tang and H. He, "ENN: Extended Nearest Neighbor Method for Multivariate Pattern Classification," 
% IEEE Computational Intelligence Magazine, 2014(under review).

load spam_training
load spam_testing

%% Using ENN method to predict the group membership of test data
K = 3;
[PredictionLabel] = ENN(TrainingData, TrainingLabel, TestingData, K);

disp(['predicted results using ENN are saved in vector "PredictionLabel".'])

if (exist('TestingLabel'))
    if (length(PredictionLabel) == length(TestingData))
        accuracy = length(find(PredictionLabel == TestingLabel)) / length(TestingLabel);
        disp(['Averaged accuracy using ENN:' num2str(accuracy)])
    end
end
