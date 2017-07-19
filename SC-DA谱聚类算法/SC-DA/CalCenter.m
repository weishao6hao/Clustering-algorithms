function [new_center,number_sample_class]=CalCenter(index_center,A,number_sample,index_class)
a=0;
temp_feature=zeros(size(A(:,1)));
for i=1:number_sample
    if (index_class(i)==index_center)
        a=a+1;
        temp_feature=temp_feature+A(:,i);
    end
end
number_sample_class=a;
if (a~=0)
    new_center=temp_feature/a;
else
    new_center=temp_feature;
end
