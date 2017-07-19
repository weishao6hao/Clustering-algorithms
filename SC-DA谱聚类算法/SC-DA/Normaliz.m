function [Y]=Normaliz(M,number_sample,k)
M_square=M.^2;%计算M中每个元素的平方
sum_M_columns=(sum(M_square,2)).^(1/2);%计算矩阵M_square的每列的和的开方值
for i=1:number_sample
    for j=1:k
        Y(i,j)=M(i,j)/sum_M_columns(i);
    end
end