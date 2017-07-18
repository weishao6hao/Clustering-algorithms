function [class,type]=dbscan(x,k,Eps)

[m,n]=size(x);

if nargin<3 | isempty(Eps)
   [Eps]=epsilon(x,k);
end


x=[[1:m]' x];
[m,n]=size(x);
type=zeros(1,m);
no=1;
touched=zeros(m,1);


for i=1:m
    if touched(i)==0;
       ob=x(i,:);
       D=dist(ob(2:n),x(:,2:n));
       ind=find(D<=Eps);
       if length(ind)<=k
          type(i)=-1;             
          class(i)=-1;             
          touched(i)=1;
       end

      if length(ind)>(k)
          while ~isempty(ind)
                ob=x(ind(1),:);
                touched(ind(1))=1;
                ind(1)=[];
                D=dist(ob(2:n),x(:,2:n));
                i1=find(D<=Eps);
				
                if length(i1)>k
                   class(i1)=no;
                   if length(i1)>k;
                      type(ob(1))=1;
                     
                   else
                      type(ob(1))=0;
                     
                   end

                   for i=1:length(i1)
                       if touched(i1(i))==0
                          touched(i1(i))=1;
                          ind=[ind i1(i)];   
                          class(i1(i))=no;
                       end                    
                   end
				   
                end
          end
          no=no+1;
      end
       end
   end
end
