% ��ȡ�˶���������
clc;
clear;
close all;

file_name = 'D:\��һ�¿γ�����\SAR�źŴ������˶�����\�����δ���ҵ\mocodata.dat';
pulse_count = 20480;

[ MOCO_UNIT_HEAD, MOCO_UNIT ] = read_mocodata( file_name, pulse_count );
disp(['�ػ��ο�prf', num2str(MOCO_UNIT_HEAD.ref_prf)]);
% �ػ����пռ�켣
figure;
plot3(MOCO_UNIT.forward, MOCO_UNIT.cross, MOCO_UNIT.height);
xlabel('forward'); ylabel('cross'); zlabel('height');
% �ػ�������й켣
hold on
MOCO_UNIT.ref_corss = zeros(1,length(MOCO_UNIT.ref_forward));%��
plot3(MOCO_UNIT.ref_forward,MOCO_UNIT.ref_corss , MOCO_UNIT.ref_height);
grid on;




% �˶����Ŀձ�����(б����������仯��

