clc;
clear;
close all;

load s_rd_20480_16384_win.mat s
% [GMG,LS,Dynamic_range,EVA,Mean,Var] = ImageEvaluation(s); 
% img_filtered_matlab = medfilt2(abs(s),[5,5]);
% [GMG1,LS1,Dynamic_range1,EVA1,Mean1,Var1] = ImageEvaluation(img_filtered_matlab); 
img=abs(s);
clear s;
[M,N]=size(img);
values = sort(img(:),'ascend');
theshold1 = values(round(0.02*M*N));
theshold2 = values(round(0.98*M*N));
img(img < theshold1) = theshold1;
img(img > theshold2) = theshold2;
figure;imshow(img); colormap('gray');