% ��ȡ�˶���������

file_name = 'E:\ѧУ\��һ��\SAR�źŴ������˶�����\h3\�˶���������\mocodata.dat';
pulse_count = 20480;

[ MOCO_UNIT_HEAD, MOCO_UNIT ] = read_mocodata( file_name, pulse_count );
disp(['�ػ��ο�prf', num2str(MOCO_UNIT_HEAD.ref_prf)]);
% �ػ����пռ�켣
figure;
plot3(MOCO_UNIT.forward, MOCO_UNIT.cross, MOCO_UNIT.height);
xlabel('forward'); ylabel('cross'); zlabel('height');
grid on;