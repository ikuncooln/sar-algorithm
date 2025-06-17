function [ y ] = A( zmed,zmin,zmax )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
A1=zmed-zmin;
A2=zmed-zmax;
if A1>0 && A2<0
    y=1;
else
    y=-1;


end

