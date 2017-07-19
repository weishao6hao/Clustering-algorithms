function  [Kappaindex]=KappaIndex(A,IDX,patternNum,k)
%记录原始属于同一类的样本的位置
OL=cell(k,1);
for i=1:k
    OL{i,1}=find(A(end,:)==i);
end
%统计原始属于同一类的样本的现在被标的类别
SLabel=cell(k,1);
for i=1:k
    number_class(i)=max(size(OL{i,1}));
    for j=1:number_class(i)
        temp=OL{i,1}(j);
        SLabel{i,1}(j)=IDX(temp);
    end
end
CL=zeros(k,k);%统计原始类别在每类中的分布情况矩阵，其中每行的行数是原始类别号
for i=1:k
    for j=1:k
        [m,n]=size(find(SLabel{i,1}==j));
        CL(i,j)=n;
    end
end
temp=0;
ON=zeros(1,k);%记录原始坐标和新的坐标的对应关系
number_correct=zeros(1,k);%记录原始坐标中的正确分类的样本数
while (temp<k)%确定样本聚类后的类别号及与原始类别号的对杨关系，以及分类正确的样本个数。
    temp=temp+1;
    maxval=max(CL(:));
    if maxval==0%假如某类中没有样本,则从剩下的还未找到对应项的标号中选择
        class_unlabel=find(ON(:)==0);%找到未标记的样本
        number_class_un=max(size(class_unlabel));%统计未标记样本的个数
        index_label_match=find(ON(:)~=0);        
        label_unmatch=1:1:k;%找到为进行匹配的标号
        number_valid=0;
        for i=1:(k-number_class_un)
            value_label_match=ON(index_label_match(i));
            for j=1:k
                if label_unmatch(j)==value_label_match
                   label_unmatch(j)=0;
                else
                number_valid=number_valid+1;%当该标号还未匹配，则定义为有效标号，则存储在label_unmatch_valid中，有效标号的个数增加1
                label_unmatch_valid(number_valid)=label_unmatch(j);
                end
            end
        end
        for i=1:number_class_un
            ON(class_unlabel(i))=label_unmatch_valid(i);%对未标号样本对应为本身的样本标号
        end
        temp=k;
    else
        [sm,sn]=find(CL==maxval);
        if max(size(sm))>1
            fsm=sm(1);
            fsn=sn(1);
        else
            fsm=sm;
            fsn=sn;
        end
        ON(1,fsm)=fsn;
        number_correct(1,fsm)=maxval;
        CL(fsm,:)=0;
        CL(:,fsn)=0;
    end
end
%记录原始数据和聚类后每类中的类别数
index_original=A(end,:);%把原始标号存入矩阵
index_new=IDX;%把聚类结果的标号存入该矩阵
number_original=zeros(k,1);%记录原始数据集每类的样本数
number_new=zeros(k,1);%记录聚类结果中每类的样本数
for i=1:k
    number_original(i)=length(find(index_original==i));
    number_new(i)=length(find(index_new==i));
end
sum_product_match=0;
for i=1:k
    sum_product_match=sum_product_match+number_original(i)*number_new(ON(i));
end
Ka_up_left=sum(number_correct)*patternNum;
Ka_up_right=sum_product_match;
Ka_down_left=patternNum*patternNum;
Ka_down_right=sum_product_match;
Kappaindex=(Ka_up_left-Ka_up_right)/(Ka_down_left-Ka_down_right);
