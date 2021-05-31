% ��ȡ�˶���������
clear; close all;clc;
% file_name = 'E:\ѧУ\��һ��\SAR�źŴ������˶�����\h3\�˶���������\mocodata.dat';

% file_name = 'D:\��һ�¿γ�����\SAR�źŴ������˶�����\�����δ���ҵ\mocodata.dat';
file_name =  'D:\Design_of_matlab\SAR\�˶���������\mocodata.dat';
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
delta_r = c/2/Fr; % ������������
lambda = wave_length;
pulse_count = azimuth_size;
center_R0 = near_range + (range_size/2)*delta_r;


[ MOCO_UNIT_HEAD, MOCO_UNIT ] = read_mocodata( file_name, pulse_count );
% disp(['�ػ��ο�prf', num2str(MOCO_UNIT_HEAD.ref_prf)]);
PRF = MOCO_UNIT_HEAD.ref_prf;
azimuth_pos = MOCO_UNIT.forward;

%% analyze
% �ػ����пռ�켣
M = length(MOCO_UNIT.forward);
xi = reshape(MOCO_UNIT.forward, [], 1);
yi = -reshape(MOCO_UNIT.cross, [], 1);
zi = reshape(MOCO_UNIT.height, [], 1);

%% �˲�����Ԥ����
% figure;
% plot(diff(MOCO_UNIT.time), 'o');
% Ӧ���Ѿ������˲����ݲ�ֵ��

% ��ϳ����뺽��
N = length(xi);
p = polyfit(xi, yi, 1);
yi_ideal = polyval(p, xi);
href = mean(zi(:)); % ����ο��߶�
zi_ideal = repmat(href, N, 1);


% ��ȡ�������
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
legend('��ʵ����','�ο�����');

%% �����ػ��켣����б�����ձ�����
% ����ƫ��
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
xlabel('������'); ylabel('��λ��'); zlabel('б�����(��)')
title('�˶��������ľ�����ͷ�λ������б�����');
%%
max_R_error = max(delta_R, [], 2);
min_R_error = min(delta_R, [], 2);
figure;
plot([max_R_error,min_R_error]);
legend('б��������ֵ', 'б�������Сֵ');grid on;
xlabel('��λ��λ��'); ylabel('�������ף�');

figure;
plot(min_R_error-max_R_error);grid on;
xlabel('��λ��λ��'); ylabel('�˶����µ�������С�������Ĳ�ף�');
%%
figure;
plot(diff(azimuth_pos));
hold on;
plot(diff(MOCO_UNIT.ref_forward));
hold off;
grid on;
legend('ʵ�ʲ��','������');
xlabel('��λ������'); ylabel('��λ�����֣���');
% title('�˶��������ķ�λ�Ǿ��Ȳ���');

% ���뷽λ����λ��
azimuth_pos_ideal = linspace(azimuth_pos(1), azimuth_pos(end), pulse_count);
figure;
plot(azimuth_pos-azimuth_pos_ideal);grid on;
xlabel('��λ������'); ylabel('ƫ�������ף�')
% title('�˶���������ʵ�ʷ�λ����λ��ƫ������λ�õ�ƫ����');
%%
figure;
plot(delta_R(round(Naz/down_rate/2), :));grid on;
xlabel('������'); ylabel('б�����ף�');
title('��λ�м�λ�ô���б��������ž�����ı仯���ף�');
%%


%  �׾�ʱ���ھ������
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
xlabel('�ϳɿ׾�����'); ylabel('��λ��'); zlabel('б�����(��)')
% title('�׾�ʱ����б�����仯');
% figure;
% plot([delta_R1(round(Naz/down_rate/2), :).',delta_R1(round(Naz/down_rate), :).']);grid on;axis tight;
% xlabel('�ϳɿ׾�����'); ylabel('б�����ף�');
% title('��λ�м�λ�ô���б��������źϳɿ׾�����ı仯���ף�');
%%
%  �׾�ʱ���ھ������
theta_bw = azimuth_angle;

% ����������һ���ϳɿ׾�ʱ���ڷ�λ��λ���仯
theta = linspace(-theta_bw/2, theta_bw/2, 550);
Ls_range = center_R0 * tan(theta);
center_R1 = sqrt((sqrt(center_R0^2-href^2) - (yi-yi_ideal)).^2 + zi.^2);
[Ls_range_grid, center_R1_grid] = meshgrid(Ls_range, center_R1);
delta_R1 = sqrt(center_R1_grid.^2 + Ls_range_grid.^2) - sqrt(center_R0^2+Ls_range_grid.^2);
phase_error1 = 4*pi*(delta_R1(end/2, :) - delta_R1(end/2, end/2))/lambda*180/pi;

% ��������Զ��һ���ϳɿ׾�ʱ���ڷ�λ��λ���仯
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
legend('�м����','��Զ����');
xlabel('����Ƕ�'); ylabel('�ձ���λ���㣩');
title('�ձ���λ���������Ƕȱ仯');
%% �����ػ��켣������λ��λ���Ŀձ�����
% ��λ��λ���
phase_error = 4*pi*delta_R / lambda;
phase_error_deg = phase_error * 180/pi;

%% ��ҵ��Ҫ��

% 1.1 ��λ��λ�ؾ�����仯
% figure;
% plot(phase_error(round(Naz/down_rate/2), :));
% xlabel('������'); ylabel('��λ��λ���');

% 1.2 �����㶯���ؾ�����仯

% 2.1 ��λ��λ�ڲ�������һ���ϳɿ׾�ʱ�����ط�λ��仯
figure;
plot(xi,href-zi);
title('����');

figure;
plot(yi);
% 2.2 �����㶯���ڲ�������һ���ϳɿ׾�ʱ�����ط�λ��仯
%%



%% ʵ�����ݴ���




%% �ձ���λ����
