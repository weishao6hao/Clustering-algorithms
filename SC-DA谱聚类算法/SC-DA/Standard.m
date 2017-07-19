function [Newdata]=Standard(A)
[m,n]=size(A);
Newdata=zeros(m,n);
for i=1:m-1
    minvalue(i)=min(A(i,:));
    maxvalue(i)=max(A(i,:));
end
for i=1:m-1
    if (maxvalue(i)-minvalue(i))==0
        Newdata(i,:)=1;
    else
        for j=1:n
            Newdata(i,j)=(A(i,j)-minvalue(i))/(maxvalue(i)-minvalue(i));
        end
    end
end

Newdata(m,:)=A(m,:);