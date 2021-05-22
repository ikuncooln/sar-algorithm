close all;
clear all;clc;
% parameters
azimuth_angle = 0.04;
wave_length = 0.03125;
chirp_rate = 200000000000000.0;
pulse_width = 2.4e-6;
range_sample_rate = 548571428.571429;
range_size = 16384;
azimuth_size = 20480;
near_range = 23306.25;
velocity = 154.195864;
prf = 533.330793;

% 4B real + 4B imag each pixel
% I,Q,I,Q
file_size = range_size * azimuth_size * 8;

%% 1. read data
x0 = 5500; y0 =5400;
% height = 20480; width = 16384;
height = 2000; width = 2000;
% data_file = 'E:/学校/研一下/SAR信号处理与运动补偿/h2/data_after_moco.dat';
data_file = 'D:\研一下课程资料\SAR信号处理与运动补偿\第二次大作业\data_after_moco.dat';

disp('读取数据中');
s0 = read_data( data_file, range_size, x0, y0, height, width);
disp('数据读取完成');

% figure;
% imagesc(real(s0));
% colormap('gray');
% title('原始信号实部');
% figure;
% subplot 211
% plot(real(s0(1,:)));
% xlabel('时间（采样点)');
% ylabel('幅度');
% subplot 212
% f=(-range_sample_rate/2:range_sample_rate/width:range_sample_rate/2-range_sample_rate/width);
% plot(f/1e6,fftshift(abs(fft(real(s0(1,:))))));
% xlabel('频率MHz');
% ylabel('幅度');
% axis tight
%% 2. convert the prameters to standard vaiables
c = 299792458;
lambda = wave_length;
f0 = c/lambda;
Kr = chirp_rate;
Vr = velocity;
Fr = range_sample_rate;
Fa = prf;
theta_rc_deg = 0;
delta_r = c/Fr/2;
center_R0 = near_range + (x0-1+width/2)*delta_r;
theta_bw = azimuth_angle;
Tr = pulse_width;
theta_rc = theta_rc_deg * pi / 180;
Naz = height;
Nrg = width;
PRF = Fa;
f_etac = 2 * Vr * sin(theta_rc) / lambda;

%% 3. 成像处理
disp('开始成像');
tic
s = RDA(s0, lambda, Kr, Vr, Fr, Fa, center_R0, theta_rc_deg,Tr);
% s = wKA( s0, lambda, Kr, Vr, Fr, Fa, center_R0, f_etac ,Tr);
% s = wKA1(s0,theta_bw,lambda,Kr,Tr,Fr,theta_rc,Nrg,Naz,near_range,Vr,PRF,0);
% s = CSA(s0,theta_bw,lambda,Kr,Tr,Fr,theta_rc,Nrg,Naz,near_range,Vr,PRF,1);
% s = RDA_smallSquintAngle(s0, lambda, Kr, Vr, Fr, PRF, theta_rc_deg ,near_range,Tr);
toc
clear s0;
img = abs(s);
% figure;imagesc(img); colormap('gray');

disp('成像完成');
%% 2%灰度增强
values = sort(img(:),'ascend');
theshold1 = values(round(0.02*Nrg*Naz));
theshold2 = values(round(0.98*Nrg*Naz));
img(img < theshold1) = theshold1;
img(img > theshold2) = theshold2;
figure;imagesc(img); colormap('gray');
img_uint8 = uint8((img-min(img(:)))/(max(img(:))-min(img(:)))*255);
clear img;
imwrite(img_uint8,'D:\研一下课程资料\SAR信号处理与运动补偿\第二次大作业\img_rd_4096.jpg');
% imwrite(img_uint8,'E:/zhaofei/repo/sar-algorithm/output/scene_rd.tiff');




% figure;
% img = abs(s);
% min_v = min(img(:));
% tmp = (img - min_v)/(max(img(:))-min_v)*255;
% imwrite(uint16(tmp), 'tmp.tif');
% img(img>255) = 255;
% imagesc(img);
% colormap('gray');
% a = abs(s);
% b = a(:);
% figure;
% plot(b);



%图像均衡化
% figure;
% imshow(histeq(abs(tmp)),[]);
% colormap('gray');
% 
% %% 自适应中值滤波
img_filtered_matlab = medfilt2(abs(s),[5,5]);
img_filtered=selfAdaption_Medianfilter(abs(s));

values = sort(img_filtered_matlab(:),'ascend');
theshold1 = values(round(0.02*Nrg*Naz));
theshold2 = values(round(0.98*Nrg*Naz));
img_filtered_matlab(img_filtered_matlab < theshold1) = theshold1;
img_filtered_matlab(img_filtered_matlab > theshold2) = theshold2;
figure;imagesc(img_filtered_matlab); colormap('gray');
img_filtered_matlab_uint8 = uint8((img_filtered_matlab-min(img_filtered_matlab(:)))/(max(img_filtered_matlab(:))-min(img_filtered_matlab(:)))*255);
imwrite(img_filtered_matlab_uint8,'D:\研一下课程资料\SAR信号处理与运动补偿\第二次大作业\img_rd_medfilt_4096.jpg');

values = sort(img_filtered(:),'ascend');
theshold1 = values(round(0.02*Nrg*Naz));
theshold2 = values(round(0.98*Nrg*Naz));
img_filtered(img_filtered < theshold1) = theshold1;
img_filtered(img_filtered > theshold2) = theshold2;
figure;imagesc(img_filtered); colormap('gray');
img_filtered_uint8 = uint8((img_filtered-min(img_filtered(:)))/(max(img_filtered(:))-min(img_filtered(:)))*255);
imwrite(img_filtered_uint8,'D:\研一下课程资料\SAR信号处理与运动补偿\第二次大作业\img_rd_adaptmedfilt_4096.jpg');


