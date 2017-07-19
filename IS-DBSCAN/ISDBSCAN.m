clear all;
k=15;
data=load('D31.txt');
[m n]=size(data);
D=zeros(m,m);temp=zeros(1,m);list=zeros(m,m);
R_num=zeros(m,15);den=zeros(m,1);
% Calculate density
for i=1:m
    temp=dist(data(i,:),data);
    [c r]=sort(temp);
    list(i,:)=r;
    D(i,:)=c;
    den(i,1)=1/D(i,k);

end
% Calculate influence space
for i=1:m
    [rc rr]=find(i==list(:,1:k));
    temp1=length(rc); 
    if temp1~=0
		R_num(i,1:temp1)=rc;
    end
end
[R L]=find(R_num);
ma=max(L);
IS=zeros(m,k);class=zeros(m,1);
R_list=zeros(m,ma-1);
for i=1:m
	R_list(i,1:ma-1)=R_num(i,2:ma);
	temp2=intersect(R_list(i,1:ma-1),list(i,2:k));
	len=length(temp2);
	IS(i,1:len)=temp2;
end
% make cluster
[class] =makecluster(IS,m,k,class);
gscatter(data(:,1),data(:,2),class,'rgbyckm','xd*+pso<>');
title('IS-DBSCAN(k=15)')
xlabel('X');
ylabel('Y');

grid on;
toc
