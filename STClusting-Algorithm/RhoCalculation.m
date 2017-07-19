function [rhoArray, distMatrix] = RhoCalculation(data,K)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% distArray = pdist(data,'cosine');%for picture
distArray = pdist(data);
distMatrix = squareform(distArray);
N = size(distMatrix, 1);
rhoArray = zeros(N,1);
if K >= N/2
    K = ceil(N/2);
end
for i = 1:N
    sortValues = sort(distMatrix(i,:));
    if mean(sortValues([1:(K+1)])) == 0
       rhoArray(i) = 10*K;
    else
       rhoArray(i) = K/sum(sortValues([1:(K+1)]));
    end
end
% 
% for i = 1:N
% %     sortValues = sort(distMatrix(i,:));
%     rhoArray(i) = mean(exp(-(distMatrix(i,:)).^2));
% end

% distArray = pdist(data);
% N = length(distArray);
% position=round(N*K/100);
% sda=sort(distArray);
% dc=sda(position);
% distMatrix = squareform(distArray);
% 
% ND = size(data,1);
% RhoArray = zeros(ND,1);
% %
% % Gaussian kernel
% for i=1:ND-1
%   for j=i+1:ND
%      RhoArray(i)=RhoArray(i)+exp(-(distMatrix(i,j)/dc)*(distMatrix(i,j)/dc));
%      RhoArray(j)=RhoArray(j)+exp(-(distMatrix(i,j)/dc)*(distMatrix(i,j)/dc));
%   end
% end
%"Cut off" kernel
% for i=1:ND-1
%  for j=i+1:ND
%    if (distMatrix(i,j)<dc)
%       RhoArray(i)=RhoArray(i)+1.;
%       RhoArray(j)=RhoArray(j)+1.;
%    end
%  end
% end

end

