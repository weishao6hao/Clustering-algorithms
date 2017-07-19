function [J]=CAA(data,beita,T)
[m n]=size(data);
for i=1:m
     temp=pdist2(data(i,:),data).^2;
     
     J(i)=sum((exp(-temp./beita)).^T);
    
end
J=J';