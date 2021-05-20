% 读取运动补偿数据

file_name = 'E:\学校\研一下\SAR信号处理与运动补偿\h3\运动补偿数据\mocodata.dat';
pulse_count = 20480;

[ MOCO_UNIT_HEAD, MOCO_UNIT ] = read_mocodata( file_name, pulse_count );
disp(['载机参考prf', num2str(MOCO_UNIT_HEAD.ref_prf)]);
% 载机飞行空间轨迹
figure;
plot3(MOCO_UNIT.forward, MOCO_UNIT.cross, MOCO_UNIT.height);
xlabel('forward'); ylabel('cross'); zlabel('height');
grid on;