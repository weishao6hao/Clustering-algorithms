clear all;
tic;
%  load data_test
%  load data_3*************NAN algorithm
data=load('t4.8k.dat');
data=normalize(data);
X=data;
[m n]=size(X);
D=pdist2(X,X);
[B I]=sort(D,2);

r=1;flag=0;
before=m;after=0;c2=0;
nb=zeros(m,1);nnr=cell(m,1);nan=cell(m,1);

% Search neighborhood
while flag==0||~isempty(find(nb==0))
    if(before-after==0)
        flag=1;
        c2=1;
    end
    
    for i=1:m
        nb(I(i,r+1),1)=nb(I(i,r+1),1)+1;
        nnr{i,1}=[nnr{i,1} I(i,r+1)];
         nan{I(i,r+1),1}=[nan{I(i,r+1),1} i];
    end
	
    if(r>=2)
        before=after;
    end
	
    r=r+1;
    after=length(find(nb));
    if(before-after==0)&&c2==1
        break;
    else
        flag==0;
    end  
end
nb_temp=nb;
nanmat=zeros(r,1);
for i=1:m
    nanmat(i,1:length(cell2mat(nan(i))))=cell2mat(nan(i));
end
toc
nnr=nanmat;
touch=zeros(m,1);
class=zeros(m,1);
lay=1;

% Create graph and clustering
while(1)
    x=find(nb==max(nb));
    x=x(1,1);
    nb(x)=-1;
    touch(x)=1;
    ind=nnr(x,1:length(find(x~=0)));
    
    while ~isempty(ind)
        class(x)=lay;
        if touch(ind(1))==0
            class(ind(1))=lay;
            touch(ind(1))=1;
            nb(ind(1))=-1;
            len1=length(find(nnr(ind(1),:)~=0));
            
            
            for i=1:len1
                if(touch(nnr(ind(1),i))==0)
                    class(nnr(ind(1),i))=lay;
                    ind=[ind nnr(ind(1),i)];
                end
            end
        end
        ind(1)=[];
    end
    lay=lay+1;
    if(isempty(find(nb~=-1)))
        break;
    end
end

figure(2)
gscatter(X(:,1),X(:,2),class,'rgbmykc','*+pso.>')
title('ENaN')

toc