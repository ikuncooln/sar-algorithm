clear all;
theta_bw_deg = 10;
theta_bw = 10*pi/180;
h = 8e3;
m1 = 5e3;
m2 = 1e4;
theta = linspace(-theta_bw/2, theta_bw/2, 550);
dx = 3;
dz = 3;
href = h;
lambda = 0.03;

% 距离向中心一个合成孔径时间内方位相位误差变化
center_R0 = sqrt(m1^2+h^2);
Ls_range = center_R0 * tan(theta);
center_R1 = sqrt((sqrt(center_R0 ^2-href^2) - dx)^2 + (dz+h)^2);
delta_R1 = sqrt(center_R1^2 + Ls_range.^2) - sqrt(center_R0 ^2+Ls_range.^2);
phase_error1 = 4*pi*(delta_R1 - delta_R1(end/2))/lambda*180/pi;

% 距离向最远端一个合成孔径时间内方位相位误差变化
max_R0 = sqrt(m2^2+h^2);
Ls_range = max_R0 * tan(theta);
center_R1 = sqrt((sqrt(max_R0 ^2-href^2) - dx)^2 + (dz+h)^2);
delta_R1 = sqrt(center_R1^2 + Ls_range.^2) - sqrt(max_R0 ^2+Ls_range.^2);
phase_error2 = 4*pi*(delta_R1 - delta_R1(end/2))/lambda*180/pi;

figure;
plot([phase_error1.', -phase_error2.']);
title('中心方位位置处斜距误差随孔径时间变化曲线');