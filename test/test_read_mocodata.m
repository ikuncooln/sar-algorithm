% 读取运动补偿数据
clc;
clear;
close all;

file_name = 'D:\研一下课程资料\SAR信号处理与运动补偿\第三次大作业\mocodata.dat';
pulse_count = 20480;

[ MOCO_UNIT_HEAD, MOCO_UNIT ] = read_mocodata( file_name, pulse_count );
disp(['载机参考prf', num2str(MOCO_UNIT_HEAD.ref_prf)]);
% 载机飞行空间轨迹
figure;
plot3(MOCO_UNIT.forward, MOCO_UNIT.cross, MOCO_UNIT.height);
xlabel('forward'); ylabel('cross'); zlabel('height');
% 载机理想飞行轨迹
hold on
MOCO_UNIT.ref_corss = zeros(1,length(MOCO_UNIT.ref_forward));%？
plot3(MOCO_UNIT.ref_forward,MOCO_UNIT.ref_corss , MOCO_UNIT.ref_height);
grid on;




% 运动误差的空变特性(斜距误差随距离变化）

