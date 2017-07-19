%minPCA   PCA the data.  
%   [frData, frPC, frVar] = minPCA(data)
%   Input data, N-by-D matrix
function [frData, frPC, frVar] = minPCA(data)
[numData,numFeatures] = size(data); 
% subtract off the mean for each dimension 
avg = mean(data,1); 
data = data - ones(numData,1) * avg; 
% calculate the covariance matrix 
covariance = 1 / (numData-1) * data.' * data; 
% find the eigenvectors and eigenvalues 
[PC, V] = eig(covariance);
% extract diagonal of matrix as vector 
V = diag(V); 
% sort the variances in decreasing order 
[junk, rindices] = sort(V,'descend');
frVar = V(rindices); frPC = PC(:,rindices); 
[c,v] = find(junk > 1e-5);
if (isempty(c))
    frNum = 1;
else
    frNum = c(end);
end

% project the original data set 
% frData = frPC.' * data.'; 
frData = data * frPC(:,1:frNum); 
end