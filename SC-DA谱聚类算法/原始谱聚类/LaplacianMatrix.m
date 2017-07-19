function [D,L]=LaplacianMatrix(A,number_sample)
D=zeros(number_sample,number_sample);
L=zeros(number_sample,number_sample);
value_diag_D=zeros(number_sample,1);
value_diag_D=sum(A,2);%对相似矩阵的每行进行加和
D=diag(value_diag_D);%得到度矩阵
L=(D^(-1/2))*A*(D^(-1/2));%得到L矩阵