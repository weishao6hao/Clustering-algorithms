function  [type] =makecluster(IS,m,k,type)
lay=1;
for p=1:m
   if(type(p)==0)
     L=length(find(IS(p,:)));

     ind=IS(p,(find(IS(p,:))>0));
    if(L>k*2/3)
        indd=[];
       while ~isempty(ind)
                if(length(find(IS(ind(1),:)))>k*2/3)
                    len1=length(find(IS(ind(1),:)));
                    for i=1:len1
                         if(type(IS(ind(1),i))==0)
                          type(IS(ind(1),i))=lay;
                          ind=[ind IS(ind(1),i)];
                         end
                    end
                end
                indd=[indd ind(1)];
              ind(1)=[];
       end
       if length(indd)>k
           type(p)=lay;
       else
           type(indd',1)=-1;
       end
       lay=lay+1;
    end
     type(p,1)=-1;
  end
end
type(find(type==0),1)=-1;

