function [ zmed,zmin,zmax] = getdata( imgdata)
%UNTITLED3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[m,n]=size(imgdata);
zmed=median(reshape(imgdata,1,m*n));
zmin=min(min(imgdata));
zmax=max(max(imgdata));
       

end

