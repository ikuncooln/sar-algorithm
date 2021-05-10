function [GMG,LS,Dynamic_range,EVA,Mean,Var] = ImageEvaluation(img)
%图像质量评估函数
%input:二维图像
%output：评价值：GMG(灰度平均梯度值)，LS(拉普拉斯算子和)，Mean(均值），Var(方差）
%{
灰度平均梯度值方法(GMG)是分别将图像长度和宽度方向上的相邻像素灰度值做差后求平方和，再求均方根值，它能较好
地反映图像的对比度和纹理变化特征，其值越大表示图像越清晰，图像质量越好.
%}
%{
拉普拉斯算子(LS)和是对每一个像素在3×3 的邻域内采用拉普拉斯算子得到8 邻域微分值，然后在图像范围内
求和，一般图像越清晰，轮廓越鲜明，则每一像素附近的灰度值变化越大，LS值就越大.
%}
%{
动态范围（英语：dynamic range）是可变化信号（例如声音或光）最大值和最小值的比值。
也可以用以10为底的对数表示。
%}
%{
清晰度(EVA)是图像细节边缘变化的敏锐程度，在图像细节的边缘处，光学密度或亮度随位置的变化越敏锐
(变化快)、越剧烈(反差大)，则细节的边缘就越清晰，可辨程度越高。
%}
%{
等效视数(ENL)：等效视数度量了图像区分具有不同后向散射特性区域的能力，是衡量一幅SAR图像斑点噪声相对强度的一种指标，
等效视数越大，表明图像上斑点越弱。其计算方法为对SAR图像均匀区域的灰度均值与方差求商再平方。
%}
[M,N]=size(img);
img=abs(img);

Column=reshape(img,1,M*N);
Mean=mean(Column);
Var=var(Column);
clear Column

%% Gray Mean Grads，GMG 灰度平均梯度值
X=zeros(M-1,N-1);
for i=1:M-1
    for j=1:N-1
        X(i,j)=sqrt(((img(i,j+1)-img(i,j))^2+(img(i+1,j)-img(i,j))^2)/2);
    end
end
GMG=sum(sum(X))/(M-1)/(N-1);
disp('GMG计算完成');
clear X
%% Laplacian 拉普拉斯算子和(LS)
laplace=[-1 -1 -1;
         -1 1 -1;
         -1 -1 -1];
res=conv2(img,laplace,'valid');
LS=sum(sum(abs(res)))/(M-2)/(N-2);
values = sort(res(:),'ascend');
theshold1 = values(round(0.02*M*N));
theshold2 = values(round(0.98*M*N));
res(res < theshold1) = theshold1;
res(res > theshold2) = theshold2;
figure;imagesc(res); colormap('gray');title('拉普拉斯算子边缘检测1');

clear res
disp('LS计算完成');
%% 动态范围（英语：dynamic range）
Imin=min(min(img));
Imax=max(max(img));
if Imin==0
    disp('图像存在为0的点');
    Dynamic_range=0;
else
    Dynamic_range=10*log10(Imax/Imin);
end
disp('Dynamic range计算完成');
%% 
tmp1 = zeros(M-2,N-2);
tmp2 = zeros(M-2,N-2);
for x = 2:M-1
    for y = 2:N-1
        tmp1(x-1,y-1)= abs(img(x-1,y)-img(x,y))+abs(img(x+1,y)-img(x,y))+abs(img(x,y+1)-img(x,y))+abs(img(x,y-1)-img(x,y));
        tmp2(x-1,y-1)= abs(img(x-1,y-1)-img(x,y))+abs(img(x+1,y+1)-img(x,y))+abs(img(x-1,y+1)-img(x,y))+abs(img(x+1,y-1)-img(x,y));
    end
end
e=tmp1+tmp2;%e(x-1,y-1)存有f(x,y)的8领域加权差值和
EVA=sum(sum(e))/(M-2)/(N-2);
disp('EVA计算完成');
clear tmp1 tmp2
%% 等效视数(ENL)
% ENL = (Mean/Var).^2;
%怎么寻找均匀区域？
end

