% clc;
% clear;
% close all;
% 
% load s_wk_20480_16384.mat
img = abs(s);
%% 2%灰度增强
values = sort(img(:),'ascend');
theshold1 = values(round(0.02*length(values)));
theshold2 = values(round(0.98*length(values)));
img(img < theshold1) = theshold1;
img(img > theshold2) = theshold2;
figure;imagesc(img); colormap('gray');
img_uint8 = uint8((img-min(img(:)))/(max(img(:))-min(img(:)))*255);
clear img;
imwrite(img_uint8,'D:\研一下课程资料\SAR信号处理与运动补偿\第二次大作业\img_rd_full_nonwin.jpg');
imwrite(img_uint8,'D:\研一下课程资料\SAR信号处理与运动补偿\第二次大作业\img_rd_full_nonwin.bmp');