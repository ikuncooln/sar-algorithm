function [ NEWimg ] = selfAdaption_Medianfilter( imgdata )
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[rows,columns]=size(imgdata);
n=3;%�趨���ڳߴ磬Ĭ����n*n�Ĵ�
nmax=15;%���ߴ�
flag=0;
%���µ�start1��start2�趨���Լ�min(j+(n-1)/2,columns)��Ϊ���ܹ��˳�ͼ��߿��ϵ�������
%���庯��AΪA���̡������������������д��������B����ֱ�ӷ���������elseif�С����zmed�Ĺ��̷���else��
for i=1:rows
    for j=1:columns
         zij=imgdata(i,j);
         if i-(n-1)/2<=0;
             start1=1;
             
         else
             start1=i-(n-1)/2;
         end
         if j-(n-1)/2<=0;
             start2=1;
             
         else
             start2=j-(n-1)/2;
         end
         [zmed,zmin,zmax]=getdata(imgdata(start1:min(i+(n-1)/2,rows),start2:min(j+(n-1)/2,columns)));
         y=A(zmed,zmin,zmax);
         if y==-1
            while (flag==4)
             n=n+2;
             [zmed,zmin,zmax]=getdata(imgdata(start1:min(i+(n-1)/2,rows),start2:min(j+(n-1)/2,columns)));
             y=A(zmed,zmin,zmax);
             s=[y==-1, n<=nmax, n<2*i-1,n<2*j-1];
             flag=sum(s);
            end
             if y==1
                 B1=zij-zmin;
                 B2=zij-zmax;
                 if B1>0 && B2<0
                     imgdata(i,j)=zij;
                 else
                     imgdata(i,j)=zmed;
                 end
             else
                 imgdata(i,j)=zmed;
             end
         elseif y==1
             B1=zij-zmin;
             B2=zij-zmax;
             if B1>0 && B2<0
                 imgdata(i,j)=zij;
             else
                 imgdata(i,j)=zmed;
             end
         else
             imgdata(i,j)=zmed;
         end
    end
             
end
NEWimg=imgdata;
% figure;
% imshow(NEWimg);
% title('����Ӧ��ֵ�˲���');

end

