function [result]=ShowFigure(X,Y,class_original,index_class,number_sample,k)
color_choose=cell(4,1);%先给出颜色备选集
color_choose{1,1}='c.';
color_choose{2,1}='y.';
color_choose{3,1}='g.';
color_choose{4,1}='m.';
figure(1)
for i=1:number_sample
    plot(X(1,i),X(2,i),color_choose{class_original(i),1});
    hold on
end
figure(2)
Y=Y';%使Y的每一列表示一个样本
for i=1:number_sample
    plot(Y(1,i),Y(2,i),color_choose{class_original(i),1});
    hold on
end
figure(3)
for i=1:number_sample
    plot(X(1,i),X(2,i),color_choose{index_class(i),1});
    hold on
end
result=0;
