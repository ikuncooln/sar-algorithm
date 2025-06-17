% 读取运动补偿数据
clear; close all;clc;
% file_name = 'E:\学校\研一下\SAR信号处理与运动补偿\h3\运动补偿数据\mocodata.dat';

% file_name = 'D:\研一下课程资料\SAR信号处理与运动补偿\第三次大作业\mocodata.dat';
file_name =  'D:\Design_of_matlab\SAR\运动补偿数据\mocodata.dat';
% basic parameters
c = 299792458;
wave_length = 0.03125;
near_range = 23306.25;
range_size = 16384;
range_sample_rate = 548571428.571429;
azimuth_size = 20480;
azimuth_angle = 0.04;

% parameters convert
Nrg = range_size;
Fr = range_sample_rate;
Naz = azimuth_size;
delta_r = c/2/Fr; % 距离向采样间距
lambda = wave_length;
pulse_count = azimuth_size;
center_R0 = near_range + (range_size/2)*delta_r;


[ MOCO_UNIT_HEAD, MOCO_UNIT ] = read_mocodata( file_name, pulse_count );
% disp(['载机参考prf', num2str(MOCO_UNIT_HEAD.ref_prf)]);
PRF = MOCO_UNIT_HEAD.ref_prf;
azimuth_pos = MOCO_UNIT.forward;

%% analyze
% 载机飞行空间轨迹
M = length(MOCO_UNIT.forward);
xi = reshape(MOCO_UNIT.forward, [], 1);
yi = -reshape(MOCO_UNIT.cross, [], 1);
zi = reshape(MOCO_UNIT.height, [], 1);

%% 运补数据预处理
% figure;
% plot(diff(MOCO_UNIT.time), 'o');
% 应该已经做了运补数据插值了

% 拟合出理想航迹
N = length(xi);
p = polyfit(xi, yi, 1);
yi_ideal = polyval(p, xi);
href = mean(zi(:)); % 理想参考高度
zi_ideal = repmat(href, N, 1);


% 抽取样点分析
down_rate = 50;
idx = 1:down_rate:M;
xi = xi(idx); yi = yi(idx); zi = zi(idx);
yi_ideal = yi_ideal(idx); zi_ideal = zi_ideal(idx);

figure;
plot3(xi, yi, zi);
xlabel('forward'); ylabel('cross'); zlabel('height');
grid on;
hold on;
plot3(xi, yi_ideal, zi_ideal);
hold off;
legend('真实航迹','参考航迹');

%% 根据载机轨迹分析斜距误差空变特性
% 距离偏差
R0 = near_range + (0:Nrg-1)*delta_r;
R0 = R0(1:down_rate:Nrg);
cos_alpha = href ./ R0;
sin_alpha = sqrt(1-cos_alpha.^2);
[cos_alpha_mtx, delta_y] = meshgrid(cos_alpha, yi-yi_ideal);
[sin_alpha_mtx, delta_z] = meshgrid(sin_alpha, zi-href);

delta_R = delta_z .* cos_alpha_mtx;
clear delta_z cos_alpha_mtx;
delta_R = delta_R - delta_y .* sin_alpha_mtx;
clear delta_y sin_alpha_mtx

figure;
mesh(delta_R);
colormap('jet');
xlabel('距离向'); ylabel('方位向'); zlabel('斜距误差(米)')
title('运动误差引起的距离向和方位向各点的斜距误差');
%%
max_R_error = max(delta_R, [], 2);
min_R_error = min(delta_R, [], 2);
figure;
plot([max_R_error,min_R_error]);
legend('斜距误差最大值', '斜距误差最小值');grid on;
xlabel('方位向位置'); ylabel('距离误差（米）');

figure;
plot(min_R_error-max_R_error);grid on;
xlabel('方位向位置'); ylabel('运动误差导致的最大和最小距离误差的差（米）');
%%
figure;
plot(diff(azimuth_pos));
hold on;
plot(diff(MOCO_UNIT.ref_forward));
hold off;
grid on;
legend('实际差分','理想差分');
xlabel('方位采样点'); ylabel('方位坐标差分：米');
% title('运动误差引起的方位非均匀采样');

% 理想方位采样位置
azimuth_pos_ideal = linspace(azimuth_pos(1), azimuth_pos(end), pulse_count);
figure;
plot(azimuth_pos-azimuth_pos_ideal);grid on;
xlabel('方位采样点'); ylabel('偏移量（米）')
% title('运动误差引起的实际方位采样位置偏离理想位置的偏移量');
%%
figure;
plot(delta_R(round(Naz/down_rate/2), :));grid on;
xlabel('距离向'); ylabel('斜距误差（米）');
title('方位中间位置处的斜距误差沿着距离向的变化（米）');
%%


%  孔径时间内距离误差
theta_bw_deg = azimuth_angle*180/pi;
Ls = center_R0 * azimuth_angle;
delta_a = MOCO_UNIT_HEAD.ref_vel/PRF;
Num = Ls/delta_a;
Ls_range = linspace(-Ls/2 , Ls/2, 550);
center_R1 = sqrt((sqrt(center_R0^2-href^2) - (yi-yi_ideal)).^2 + zi.^2);
[Ls_range_grid, center_R1_grid] = meshgrid(Ls_range, center_R1);
delta_R1 = sqrt(center_R1_grid.^2 + Ls_range_grid.^2) - sqrt(center_R0^2+Ls_range_grid.^2);
figure;
mesh(delta_R1);
colormap('jet');
xlabel('合成孔径方向'); ylabel('方位向'); zlabel('斜距误差(米)')
% title('孔径时间内斜距误差变化');
% figure;
% plot([delta_R1(round(Naz/down_rate/2), :).',delta_R1(round(Naz/down_rate), :).']);grid on;axis tight;
% xlabel('合成孔径方向'); ylabel('斜距误差（米）');
% title('方位中间位置处的斜距误差沿着合成孔径方向的变化（米）');
%%
%  孔径时间内距离误差
theta_bw = azimuth_angle;

% 距离向中心一个合成孔径时间内方位相位误差变化
theta = linspace(-theta_bw/2, theta_bw/2, 550);
Ls_range = center_R0 * tan(theta);
center_R1 = sqrt((sqrt(center_R0^2-href^2) - (yi-yi_ideal)).^2 + zi.^2);
[Ls_range_grid, center_R1_grid] = meshgrid(Ls_range, center_R1);
delta_R1 = sqrt(center_R1_grid.^2 + Ls_range_grid.^2) - sqrt(center_R0^2+Ls_range_grid.^2);
phase_error1 = 4*pi*(delta_R1(end/2, :) - delta_R1(end/2, end/2))/lambda*180/pi;

% 距离向最远端一个合成孔径时间内方位相位误差变化
max_R0 = center_R0 + Nrg/2*delta_r;
Ls = max_R0  * azimuth_angle;
delta_a = MOCO_UNIT_HEAD.ref_vel/PRF;
Num = Ls/delta_a;
Ls_range = max_R0 * tan(theta);
center_R1 = sqrt((sqrt(max_R0 ^2-href^2) - (yi-yi_ideal)).^2 + zi.^2);
[Ls_range_grid, center_R1_grid] = meshgrid(Ls_range, center_R1);
delta_R1 = sqrt(center_R1_grid.^2 + Ls_range_grid.^2) - sqrt(max_R0 ^2+Ls_range_grid.^2);
phase_error2 = 4*pi*(delta_R1(end/2, :) - delta_R1(end/2, end/2))/lambda*180/pi;

figure;
plot(theta*180/pi,[phase_error1.', phase_error2.']);grid on;axis tight;
legend('中间距离','最远距离');
xlabel('照射角度'); ylabel('空变相位误差（°）');
title('空变相位误差随照射角度变化');
%% 根据载机轨迹分析方位相位误差的空变特性
% 方位相位误差
phase_error = 4*pi*delta_R / lambda;
phase_error_deg = phase_error * 180/pi;

%% 作业题要求

% 1.1 方位相位沿距离向变化
% figure;
% plot(phase_error(round(Naz/down_rate/2), :));
% xlabel('距离向'); ylabel('方位相位误差');

% 1.2 距离徙动量沿距离向变化

% 2.1 方位相位在测绘带中心一个合成孔径时间内沿方位向变化
figure;
plot(xi,href-zi);
title('泡面');

figure;
plot(yi);
% 2.2 距离徙动量在测绘带中心一个合成孔径时间内沿方位向变化
%%



%% 实际数据处理




%% 空变相位补偿
