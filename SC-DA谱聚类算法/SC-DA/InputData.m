function [X,number_sample,k,class_original]=InputData(route)
X_original=importdata(route);%每一列是一个数据对象
X_original=Standard(X_original);
[m,n]=size(X_original);%m维特征，n个数据对象
number_sample=n;
class_original=X_original(m,:);
k=max(class_original);
X=X_original(1:(m-1),:);
