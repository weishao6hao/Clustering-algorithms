function [A]=AffinityMatrix(X,number_sample,k,Standard_Deviation,CNN)
A=zeros(number_sample,number_sample);
for i=1:(number_sample-1)
    for j=(i+1):number_sample
        distance_square=DistanceSquare(X(:,i),X(:,j));
        A(i,j)=exp((-distance_square)/(2*(Standard_Deviation^2)*(CNN(i,j)+1)));
        A(j,i)=A(i,j);
    end
    A(i,i)=0;
end
