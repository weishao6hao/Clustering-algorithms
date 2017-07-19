%评价指标Adjusted Rand Index
function [ARI]=AdjustedRandIndexM(A,IDX,patternNum,k)
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
number_mutual_trans=number_mutual';%把矩阵转换成：行代表划分后的类别，列代表原始类别
%ARI公式中的四项
ARI_up_left=0;
ARI_up_right=0;
ARI_down_left=0;
ARI_down_right=0;
%计算分子的左侧项
for i=1:k_new
    for j=1:k_original
        ARI_up_left=ARI_up_left+number_mutual_trans(i,j)*(number_mutual_trans(i,j)-1)/2;
    end
end
%计算分子的右侧项
sum_number_new_2_pick=0;%计算划分得到的结果中在同一类的点对和
for i=1:k_new
    sum_number_new_2_pick=sum_number_new_2_pick+number_new(i)*(number_new(i)-1)/2;
end
sum_number_original_2_pick=0;%计算原始数据中在同一类的点对和
for j=1:k_original
    sum_number_original_2_pick=sum_number_original_2_pick+number_original(j)*(number_original(j)-1)/2;
end
total_pairs=(patternNum*(patternNum-1))/2;%计算总的点对数
ARI_up_right=sum_number_new_2_pick*sum_number_original_2_pick/total_pairs;
%计算分母的左侧项
ARI_down_left=(sum_number_new_2_pick+sum_number_original_2_pick)/2;
%计算分母的右侧项
ARI_down_right=sum_number_new_2_pick*sum_number_original_2_pick/total_pairs;
ARI=(ARI_up_left-ARI_up_right)/(ARI_down_left-ARI_down_right);