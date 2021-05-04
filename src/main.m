close all;clear all;
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
x0 = 1; y0 = 1;
height = 2048; width = 2048;
data_file = 'E:/学校/研一下/SAR信号处理与运动补偿/h2/data_after_moco.dat';
s0 = read_data( data_file, range_size, x0, y0, height, width);

figure;
imagesc(real(s0));
colormap('gray');
figure;
plot(real(s0(1,:)));

%% 2. convert the prameters to standar vaiables
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
%% . focus
% s = rd_big_func(s0, f0, Kr, Vr, Fr, Fa, center_R0, theta_rc_deg );
s = omega_k_func( s0, f0, Kr, Vr, Fr, Fa, center_R0, 0 );
figure;
img = abs(s);
min_v = min(img(:));
tmp = (img - min_v)/(max(img(:))-min_v)*255;
imwrite(uint16(tmp), 'tmp.tif');

img(img>255) = 255;
imagesc(img);
colormap('gray');



a = abs(s);
b = a(:);
figure;
plot(b);








