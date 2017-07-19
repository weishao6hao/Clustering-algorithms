function [rad_density]=GetRad(X,number_sample)
%计算样本间的距离
 Distance_sample=zeros(number_sample,number_sample);
for i=1:(number_sample-1)
    for j=(i+1):number_sample
        Distance_sample(i,j)=sqrt(DistanceSquare(X(:,i),X(:,j)));%计算样本间的欧式距离
        Distance_sample(j,i)=Distance_sample(i,j);
    end
end
mean_distance=sum(sum(Distance_sample))/((number_sample-1)*number_sample);%计算平均值
max_distance=max(Distance_sample(:));%计算距离矩阵的最大值
for i=1:number_sample      
    Distance_sample(i,i)=inf;% 为了计算样本间距离的最小值，则为距离矩阵的对角线赋值为INF
end
min_distance=min(Distance_sample(:));%计算距离矩阵的最小值
distance_nearest=min(Distance_sample);%计算样本和其最近邻间的距离
mean_nearest=mean(distance_nearest);%计算样本和其最近邻间的距离的均值
max_nearest=max(distance_nearest);%计算样本和其最近邻间的距离的最大值
rad_density=20*mean_distance+54*min_distance+13*max_nearest-6*max_distance-65*mean_nearest;%计算密度半径参数