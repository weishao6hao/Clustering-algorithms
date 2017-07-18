function [class nn]=Nature_neighbor(data)
% load data1_3N
X=data;
[m n]=size(X);
D=pdist2(X,X);
[B I]=sort(D,2);
toc;
count=1;
c2=0;
r=1;flag=0;
before=m;after=0;
nb=zeros(m,1);nnr=cell(m,1);
while flag==0||~isempty(find(nb==0))

    if(before-after==0)
        flag=1;
        c2=1;
    end
    
    for i=1:m
        nb(I(i,r+1),1)=nb(I(i,r+1),1)+1;
        nnr{i,1}=[nnr{i,1} I(i,r+1)];
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

toc
nnr=cell2mat(nnr);
touch=zeros(m,1);
class=zeros(m,1);
len=size(nnr,2);
nnnr(:,1:ceil(len*2/3))=nnr(:,1:ceil(len*2/3));
lay=1;
% Clustering from the highest point of density
while(1)
    x=find(nb==max(nb));
    x=x(1,1);
    nb(x)=-1;
    temp=find(class(nnr(x,1:len))~=0);
    if ~isempty(temp)
        class(x)=class(nnr(x,temp(1,1)));
    else
        class(x)=lay;
    end
    ind=nnnr(x,:);
    %     if(L>k*2/3)
    
    while ~isempty(ind)
        if nb(ind(1))~=-1
            temp=find(class(nnr(ind(1),1:len))~=0);
            if ~isempty(temp)
                class(ind(1))=class(nnr(ind(1),temp(1,1)));
            else
                class(ind(1))=lay;
            end
			
            nb(ind(1))=-1;

            
            len1=length(find(nnnr(ind(1),:)));

            for i=1:len1
                if(nb(nnnr(ind(1),i))~=-1)
                    ind=[ind nnnr(ind(1),i)];
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
nn=size(nnnr,2);
toc