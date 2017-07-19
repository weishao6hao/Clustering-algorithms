function [labels,centerIndexes, rhoArray, deltaArray] = STClusteringAlgorithm(data)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
N = size(data, 1);
% K = N/7;
% K = ceil(log(N)/log(2)) + 1;
%K = ceil(log(N)/log(2));
% K = max([ceil(log(N)/log(2)) + 1,20]);
% N = size(data,1);
% K = ceil(log(N)/log(2)) + 1;
K = round(sqrt(N));
%     K = 10;
[rhoArray, distMatrix] = RhoCalculation(data,K);
% [rhoArray, distMatrix] = RhoCalculationAA(data, 2);
[deltaArray, nneighArray] = DeltaCalculation(distMatrix, rhoArray);
gama = rhoArray.*deltaArray;
% plot(rhoArray, deltaArray,'.')
% plot(sort(gama),'r.')
% hist(gama)
[centerIndexes] = OutwardStatTestCenterDetection(gama, 0.05);

[labels] = CluteringBasedOnCenter(rhoArray, nneighArray, centerIndexes);
end

