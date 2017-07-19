function [CNN]=CommonNearNeighbor(X,number_sample,rad_density)
%计算样本间的距离
Distance_sample=zeros(number_sample,number_sample);
for i=1:(number_sample-1)
    for j=(i+1):number_sample
        Distance_sample(i,j)=sqrt(DistanceSquare(X(:,i),X(:,j)));%计算样本间的欧式距离
        Distance_sample(j,i)=Distance_sample(i,j);
    end
    Distance_sample(i,i)=inf;%为了方便找样本的近邻
end
%找到每个样本密度半径内的近邻样本
Neighbor=zeros(number_sample,number_sample);
[index_neighbor]=find(Distance_sample<=rad_density);%找到样本的近邻
Neighbor(index_neighbor)=1;%把近邻矩阵中的是近邻的样本赋值为1
CNN=zeros(number_sample,number_sample);
CNN=Neighbor*Neighbor;%计算得到两个样本的近邻的交集个数