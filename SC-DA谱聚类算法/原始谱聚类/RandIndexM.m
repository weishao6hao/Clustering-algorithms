function [RI]=RandIndexM(A,IDX,patternNum,k)
index_original=A(end,:);%把原始标号存入矩阵
index_new=IDX;%把聚类结果的标号存入该矩阵
k_original=max(index_original);%原始类别数
k_new=max(index_new);%聚类得到的类别数
number_original=zeros(k_original,1);%记录原始数据集每类的样本数
number_new=zeros(k_new,1);%记录聚类结果中每类的样本数
number_mutual=zeros(k_original,k_new);%记录原始类别和聚类后的类别之间的交集，行代表原始类别数，列代表聚类的类别数
for i=1:patternNum
    x_temp=index_original(i);%得到该样本的原始类别号
    y_temp=index_new(i);%得到该样本的聚类后的类别号
    number_original(x_temp)=number_original(x_temp)+1;
    number_new(y_temp)=number_new(y_temp)+1;
    number_mutual(x_temp,y_temp)=number_mutual(x_temp,y_temp)+1;
end
%计算总的点对数
total_pairs=(patternNum*(patternNum-1))/2;
%计算原始同一类中的对象被划分到同一类的点对的数目TP
%计算同一类对象被划分到不同类的点对数目FN
TP=0;
FN=0;
for i=1:k_new
    for j=1:k
        temp_mutual=number_mutual(j,i);
        TP=TP+(temp_mutual*(temp_mutual-1))/2;
        FN=FN+(temp_mutual*(number_original(j)-temp_mutual))/2;
    end
end
%计算不同类之间的点对总数
all_diff=0;
for i=1:k_new
    all_diff=all_diff+number_new(i)*(patternNum-number_new(i))/2;
end
%计算原始不同类的对象被划分到不同类的点对的数目
TN=all_diff-FN;
%计算评价指标RandIndex的值
RI=(TP+TN)/total_pairs;