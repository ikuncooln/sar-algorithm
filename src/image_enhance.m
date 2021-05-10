% 图像增强探究
path = 'E:/zhaofei/repo/sar-algorithm/output/';

% filename1 = 'scene1_wk_nofilter_2048_2048.bmp';
% filename2 = 'scene1_wk_Leefilter_2048_2048.bmp';
% filename3 = 'scene1_wk_medfilter_2048_2048.bmp';
% filename4 = 'scene1_wk_adpmedfilter_2048_2048.bmp';

filename1 = 'scene2_wk_nofilter_4096_4096.bmp';
filename2 = 'scene3_wk_Leefilter_4096_4096.bmp';
filename3 = 'scene2_wk_medfilter_4096_4096.bmp';
filename4 = 'scene2_wk_adpmedfilter_4096_4096.bmp';
im =  imread([path, filename1]);
% figure;
% imshow(im);
% title('未滤波');
% 
% % 中值滤波
% im_median_filtered = medfilt2(im);
% imwrite(im_median_filtered, [path, filename3]);
% figure;
% imshow(im_median_filtered);
% title('中值滤波');
% 
% % 自适应中值滤波
% im_median_adp_filtered = medfilt2(im);
% imwrite(im_median_adp_filtered, [path, filename4]);
% figure;
% imshow(im_median_adp_filtered);
% title('自适应中值滤波');

% Lee滤波
im_Lee_filtered = Lee_filter(im);
imwrite(im_Lee_filtered, [path, filename2]);
figure;
imshow(im_Lee_filtered);
title('Lee滤波更新');

