function  [z]=SCA(data,beita,gamma)
[m n]=size(data);
sigma=0.5;
T=2;
z_temp(1,1)={zeros(m,n)};
for i=1:m
    for j=1:m
         temp=dist(data(i,:),data(j,:)).^2;
         s(i,j)=exp(-temp./beita);
    end
end
while T
   for i=1:m
       for   k=1:n
       X(:,k)=(s(i,:).^gamma)'.*data(:,k);
       end
       z(i,:)=sum(X)/sum(s(i,:).^gamma);
   end
   z_temp(T,1)={z};
   if max(dist_multidimension(z_temp(T,1),z_temp(T-1,1)))<sigma
       break;
   end
   max(dist_multidimension(z_temp(T,1),z_temp(T-1,1)))
   T=T+1;
   for i=1:m
       for j=1:m
          temp=dist(z(i,:),data(j,:)).^2;
          s(i,j)=exp(-temp./beita);
       end
   end
   T

end
