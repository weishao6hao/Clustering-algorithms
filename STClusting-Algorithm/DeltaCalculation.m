function [deltaArray,nneighArray] = DeltaCalculation(distMatrix, rhoArray)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

N = size(distMatrix,1);
maxd=max(max(distMatrix));%两两间距离的最大值，最大距离
[~,ordrho]=sort(rhoArray,'descend');%将对象根据密度从大到小进行排序
deltaArray = zeros(N,1);%用于记录每个对象对应的距离
nneighArray = zeros(N,1);%用于记录每个对象的最近邻
%以下代码为计算delta的代码
deltaArray(ordrho(1))=-1.;%密度最大对象对应的初始化
nneighArray(ordrho(1))=0; %不存在近邻

for ii=2:N
   deltaArray(ordrho(ii))=maxd;
   for jj=1:ii-1%扫描所有比该对象对应密度大的对象，找出最小距离以及对应的最近邻
     if(distMatrix(ordrho(ii),ordrho(jj))<deltaArray(ordrho(ii)))
        deltaArray(ordrho(ii))=distMatrix(ordrho(ii),ordrho(jj));
        nneighArray(ordrho(ii))=ordrho(jj);
     end
   end
end
deltaArray(ordrho(1))=max(deltaArray(:));%设置为最大的值
end