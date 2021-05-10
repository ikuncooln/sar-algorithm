% Í¼ÏñÔöÇ¿Ì½¾¿
path = 'E:/zhaofei/repo/sar-algorithm/output/';

filename1 = 'scene_wk_nofilter_2048_2048.bmp';
filename2 = 'scene_wk_Leefilter_2048_2048.bmp';
im =  imread([path, filename1]);
figure;
imshow(im);

% LeeÂË²¨
im_Lee_filtered = Lee_filter(im);
imwrite(im_Lee_filtered, [path, filename2]);