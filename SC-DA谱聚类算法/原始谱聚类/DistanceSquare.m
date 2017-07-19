function [distance]=DistanceSquare(x,y)
distance=0;
distance=sum((x-y).^2);