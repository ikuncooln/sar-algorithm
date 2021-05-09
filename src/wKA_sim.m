close all;clear all;
%% 1. ������� (�ο� p142, table 6.1)
mode = 1; % 0: ��б�ӽǲ�����1��Сб�ӽǲ���

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
if mode == 0
    theta_rc_deg = 21.9; % ��б�ӽ�21.9��
else
    theta_rc_deg = 3.5; % Сб�ӽ�21.9��
end
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

figure; % ���ƴ�б�ӽ�����µ������״�ԭʼ�����ź�
subplot(221);imagesc(real(s0));ylabel('��λ��ʱ�䣨�����㣩');title('(a)ʵ��');
subplot(222);imagesc(imag(s0));title('(b)�鲿');
subplot(223);imagesc(abs(s0));xlabel('������ʱ�䣨�����㣩');ylabel('��λ��ʱ�䣨�����㣩');title('(c)����');
subplot(224);imagesc(angle(s0));xlabel('������ʱ�䣨�����㣩');title('(d)��λ');
suptitle([num2str(theta_rc_deg), '��б�ӽ�����µ�', num2str(NUM_TARGETS),'���״�ԭʼ�����źţ�ʱ��']);

%% 3. ������
% step1: do something

% step2: try something

% step3: complete

%% 4. ��Ŀ�����



