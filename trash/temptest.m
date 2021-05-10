clc;
clear;
close all;

load s_rd_20480_16384_win.mat s
[GMG,LS,Dynamic_range,EVA,Mean,Var] = ImageEvaluation(s); 

% values = sort(img(:),'ascend');
% theshold1 = values(round(0.02*Nrg*Naz));
% theshold2 = values(round(0.98*Nrg*Naz));
% img(img < theshold1) = theshold1;
% img(img > theshold2) = theshold2;
% figure;imagesc(img); colormap('gray');