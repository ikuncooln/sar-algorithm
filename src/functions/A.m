function [ y ] = A( zmed,zmin,zmax )
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
A1=zmed-zmin;
A2=zmed-zmax;
if A1>0 && A2<0
    y=1;
else
    y=-1;


end

