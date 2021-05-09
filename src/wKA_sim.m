close all;clear all;
%% 1. ������� (�ο� p142, table 6.1)
center_Rc = 20e3;  % ������б��
Vr = 150;       % ��Ч�״��ٶ�
Tr = 2.5e-6;    % ��������ʱ��
Kr = 40e12;     % �����Ƶ��
f0 = 5.3e9;     % �״﹤��Ƶ��
BW_dop = 80;    % �����մ���
Fr = 120e6;  % ���������
Fa = 100;   % ��λ������
Naz = 256;  % ��λ�����������������������
Nrg = 3072;  % ��������������������߲���������
theta_rc_deg = 5; % ��б�ӽ�5��
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
NUM_TARGETS = 7;    % �����Ŀ����
gap = 50 * 8;
% gap = 25;
rs = [-3, -2, -1, 0, 1, 2, 3]*gap;
as = [-3, -2, -1, 0, 1, 2, 3]*gap*tan(theta_rc);
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

%% 3. ������
Rref = center_R0;
s2 = wKA( s0, lambda, Kr, Vr, Fr, Fa, Rref, f_etac, Tr );


%% show
Vg = Vr;
[Naz, Nrg] = size(s2);
x = ((-Nrg / 2) : (Nrg / 2 - 1)) / Fr * c / 2;
y = ((-Naz / 2 : Naz / 2 - 1)) / Fa * Vg;
figure;
imagesc(x, y, abs(s2));
xlabel('������m��');ylabel('��λ��m��');
title('����������Ŀ��');set(gca, 'YDir', 'normal');

%% 4. ��Ŀ�����
% ����ÿ�������λ�õ�����ֵ
ns = round(rs/delta_r) + (Nrg/2 + 1);
ms = round(as/delta_a) + (Naz/2 + 1);
len = 16;

figure; % �Ŵ���ʾÿ����Ŀ��
for i = 1:NUM_TARGETS
    target = s2(ms(i)-len/2:ms(i)+len/2-1, ns(i)-len/2:ns(i)+len/2-1);
    subplot(ceil(NUM_TARGETS/3), 3, i);
    imagesc(abs(target));
    xlabel('�����򣨲����㣩');ylabel('��λ�򣨲����㣩');ylabel('��λ�򣨲����㣩');title(['Ŀ��', num2str(i)]);
end

% ������������Ŀ��A������һ����Ŀ�꣩
p = 1;
target = s2(ms(p)-len/2:ms(p)+len/2-1, ns(p)-len/2:ns(p)+len/2-1);
[image_upsample,signal_r,quality_r,signal_a,quality_a] = f_point_analyse(target,delta_r,delta_a);

BW_r= abs(Kr*Tr);
La = 0.886 * 2 * Vs * cos(theta_rc) / BW_dop;   % ���߿׾�����
IRW_r_theory = c/2/BW_r*0.886*1.18;
IRW_a_theory = La/2*Vg/Vs*1.185;
disp(['���������۷ֱ���:',num2str(IRW_r_theory),'m']);
disp(['��λ�����۷ֱ���:',num2str(IRW_a_theory),'m']);

