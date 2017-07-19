function [MSN]=SchmidtNormaliz(V)
number_vector=size(V,2);%获得特征向量的个数
for i=1:number_vector
    U(:,i)=V(:,i);
    if i==1
        MSN(:,i)=U(:,i)/sqrt(dot(U(:,i),U(:,i))); 
    else
        temp=U(:,i);
        for j=1:(i-1)
            temp=temp-dot(V(:,i),U(:,j))/dot(U(:,j),U(:,j))*U(:,j);
        end
        U(:,i)=temp;
        if dot(U(:,i),U(:,i))==0
            MSN(:,i)=zeros(size(U(:,i)));
        else
            MSN(:,i)=U(:,i)/sqrt(dot(U(:,i),U(:,i))); 
        end
    end
end