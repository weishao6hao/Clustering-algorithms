%Normalize Data: normalize data to the range of 0 to 1. 
function [inputData] = normalizeData(inputData)
    [C,V] = size(inputData);
    numFeatures = V;
    for k = 1:numFeatures
        tmp = inputData(:,k);
        if (max(tmp) - min(tmp) <= 1e-5)
            inputData(:,k) = zeros(C,1);
        else
            inputData(:,k) = (tmp - min(tmp)) / (max(tmp) - min(tmp));
        end
    end
end