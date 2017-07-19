function [nmi]=NormalizedMutualInformation(A,IDX,patternNum,k)
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
sum_up=0;
sum_down_left=0;
sum_down_right=0;
for i=1:k_new
    if number_new(i)>0
        if number_new(i)==patternNum
            sum_down_right=1;
        else
            sum_down_right=sum_down_right+number_new(i)*log(number_new(i)/patternNum);%计算NMI公式分子的右半部分
        end
    end
end
for i=1:k_original
    sum_down_left=sum_down_left+number_original(i)*log(number_original(i)/patternNum);%计算NMI公式分子的左半部分    
    for j=1:k_new
        if number_mutual(i,j)==0
            temp=0;
        else
            temp=number_mutual(i,j)*log((patternNum*number_mutual(i,j))/(number_original(i)*number_new(j)));
        end
        sum_up=sum_up+temp;%计算NMI公式的分母
    end
end
nmi=sum_up/sqrt(sum_down_left*sum_down_right);
    
