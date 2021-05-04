% Test for generate_point_data() function
close all;clear all;

%% 1. ������� (�ο� p142, table 6.1)
center_Rc = 20e3;  % ������б��
Vr = 150;   % ��Ч�״��ٶ�
Tr = 2.5e-6;    % ��������ʱ��
Kr = 20e12; % �����Ƶ��
f0 = 5.3e9; % �״﹤��Ƶ��
BW_dop = 80;    % �����մ���
Fr = 60e6;  % ���������
Fa = 100;   % ��λ������
Naz = 256;  % ��λ�����������������������
Nrg = 256;  % ��������������������߲���������
theta_rc_deg = 21.9; % ��б�ӽ�21.9��
c = 299792458;    % ����

% derived params
lambda = c / f0;
theta_rc = theta_rc_deg * pi / 180;
Vs = Vr;
Vg = Vr;
Np = Tr * Fr;   % �������г��ȣ�����������
alpha_os_r = Fr / (Kr*Tr);
alpha_os_a = Fa / BW_dop;
%% 2. ����ԭʼ�״�����
NUM_TARGETS = 3;    % �����Ŀ����Ϊ3
rs = [0, 0, 30];    % ��Ŀ����������
as = [-20, 0, -10]; % Ŀ����Է�λ�����
parameters = struct(...
    'center_Rc', center_Rc,...          % ������б��
    'theta_rc_deg', theta_rc_deg,...    % б�ӽ�
    'Nrg', Nrg,...                      % �������������
    'Naz', Naz,...                      % ��λ���������
    'Vr', Vr,...                        % �ػ��ٶ�
    'f0', f0,...                        % �ز�Ƶ��
    'Tr', Tr,...                        % ����������
    'Kr', Kr,...                        % ���������Ƶ��
    'BW_dop', BW_dop,...                % �����մ���
    'alpha_os_r', alpha_os_r,...        % �������������
    'alpha_os_a', alpha_os_a,...        % ��λ���������
    'NUM_TARGETS', NUM_TARGETS,...      % ��Ŀ������
    'rs', rs,...                        % ��Ŀ����������꣨m��
    'as', as...                         % ��Ŀ�귽λ�����꣨m��
);

[ s0, f_etac, delta_r, delta_a, center_R0, center_Rc ] = generate_point_data(parameters);
figure; % �����״�ԭʼ�����ź�
imagesc(abs(s0));xlabel('������ʱ�䣨�����㣩');ylabel('��λ��ʱ�䣨�����㣩');
title('��Ŀ���״�ԭʼ�����źŷ��ȣ�ʱ��');

%% 3. ����
s = rd_big_func(  s0, f0, Kr, Vr, Fr, Fa, center_R0, theta_rc_deg );
x = ((-Nrg / 2) : (Nrg / 2 - 1)) / Fr * c / 2;
y = ((-Naz / 2 : Naz / 2 - 1)) / Fa * Vg ;
figure; % ���Ƶ�б�ӽ�����¾���ѹ���ҷ�λѹ�����ź�
subplot(121);imagesc(real(s));xlabel('������ʱ�䣨�����㣩');ylabel('��λ��ʱ�䣨�����㣩');title('(a)ʵ��');
subplot(122);imagesc(x, y, abs(s));xlabel('������ʱ�䣨�����㣩');title('(b)����');set(gca, 'YDir', 'normal');
suptitle('21.9��б�ӽǾ���ѹ���ҷ�λѹ������źţ�ʱ��');

%% 4. ����
