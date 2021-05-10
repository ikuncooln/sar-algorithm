function [ zmed,zmin,zmax] = getdata( imgdata)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
[m,n]=size(imgdata);
zmed=median(reshape(imgdata,1,m*n));
zmin=min(min(imgdata));
zmax=max(max(imgdata));
       

end

