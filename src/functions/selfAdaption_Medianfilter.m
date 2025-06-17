function [ NEWimg ] = selfAdaption_Medianfilter( imgdata )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
[rows,columns]=size(imgdata);
n=5;%设定窗口尺寸，默认是n*n的窗
nmax=15;%最大尺寸
flag=0;
%以下的start1和start2设定，以及min(j+(n-1)/2,columns)是为了能够滤除图像边框上的噪声。
%定义函数A为A过程。（需迭代调用所以另写函数）。B过程直接放在主程序elseif中。输出zmed的过程放在else中
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
% title('自适应中值滤波器');

end

