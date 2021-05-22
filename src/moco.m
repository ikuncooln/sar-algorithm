close all; clear;

%% basic parameters
moco_file = 'E:\ѧУ\��һ��\SAR�źŴ������˶�����\h3\�˶���������\mocodata.dat';
file1 = 'E:\ѧУ\��һ��\SAR�źŴ������˶�����\h3\�˶���������\data_before_moco.dat';
file2 = 'E:\zhaofei\repo\sar-algorithm\output\moco\data_after_1th_phase_compensation.dat';
file3 = 'E:\zhaofei\repo\sar-algorithm\output\moco\data_after_range_resample.dat';
file4 = 'E:\zhaofei\repo\sar-algorithm\output\moco\data_after_azimuth_resample.dat';

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

%% ѡ��˴δ������������С
last_pulse_count = 0;
pulse_count = 2048;
MAX_MEM_GB = 1; % 1�ν���������1GB����

%% һ����λ����
first_order_phase_error_compensation(moco_file,...
    file1, file2, wave_length, near_range, range_sample_rate,...
    range_size, azimuth_size, pulse_count, last_pulse_count, MAX_MEM_GB);

%% �����ز�������/�����㶯����
range_resample(moco_file,...
    file2, file3, wave_length, near_range, range_sample_rate,...
    range_size, azimuth_size, pulse_count, last_pulse_count, MAX_MEM_GB);

%% ��λ�ز�������
azimuth_resample( moco_file, file3, file4,...
    range_size, azimuth_size, pulse_count, last_pulse_count, MAX_MEM_GB);

%% �������
% ת��main.m
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
% ����λ��ѡ��
x0 = 1; y0 =1;
% height = 20480; width = 16384;
height = 2048; width = 2048;

% parameters convert
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

disp('��ʼ����...');
s0 = read_data( file4, range_size, x0, y0, height, width);
disp('���ݶ�ȡ���');
figure;imagesc(real(s0)); colormap('gray');


s = CSA(s0,theta_bw,lambda,Kr,Tr,Fr,theta_rc,Nrg,Naz,near_range,Vr,PRF,0);
img = abs(s);
% figure;imagesc(img); colormap('gray');
disp('�������');

% 2%�Ҷ���ǿ
values = sort(img(:),'ascend');
theshold1 = values(round(0.02*Nrg*Naz));
theshold2 = values(round(0.98*Nrg*Naz));
img(img < theshold1) = theshold1;
img(img > theshold2) = theshold2;
img_uint8 = uint8((img-min(img(:)))/(max(img(:))-min(img(:)))*255);
figure;imagesc(img_uint8); colormap('gray');

figure;
plot(img(:));
