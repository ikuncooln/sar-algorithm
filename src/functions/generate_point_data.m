function [  s0, f_etac, delta_r, delta_a, center_R0, center_etac  ] = generate_point_data( parameters )
%GENERATE_POINT_DATA generate SAR target raw data
%   parameters is a struct (including kinds of simulation parameters as its
%   properties.
%   Tips:
%       Please give a complete parameters struct.
%       If parameters.NUM_TARGETS equals to zero, this function will
%       display a scene, guiding you to set proper targets.
%       
%   ���磺
%         NUM_TARGETS = 3;    % �����Ŀ����Ϊ3
%         rs = [0, 0, 20];    % ��Ŀ����������
%         as = [-20, 0, -10]; % Ŀ����Է�λ�����
%         parameters = struct(...
%             'center_Rc', center_Rc,...          % ������б��
%             'theta_rc_deg', theta_rc_deg,...    % б�ӽ�
%             'Nrg', Nrg,...                      % �������������
%             'Naz', Naz,...                      % ��λ���������
%             'Vr', Vr,...                        % �ػ��ٶ�
%             'f0', f0,...                        % �ز�Ƶ��
%             'Tr', Tr,...                        % ����������
%             'Kr', Kr,...                        % ���������Ƶ��
%             'BW_dop', BW_dop,...                % �����մ���
%             'alpha_os_r', alpha_os_r,...        % �������������
%             'alpha_os_a', alpha_os_a,...        % ��λ���������
%             'NUM_TARGETS', NUM_TARGETS,...      % ��Ŀ������
%             'rs', rs,...                        % ��Ŀ����������꣨m��
%             'as', as...                         % ��Ŀ�귽λ�����꣨m��
%         );
%         [ s0, f_etac, delta_r, delta_a, center_R0, center_Rc ] = generate_point_data(parameters);
%
%   ����ֵ:
%   s0 ����ز��źž���
%   f_etac �źŶ���������Ƶ��
%   delta_r ���ݾ����������ࣨm��
%   delta_a ���ݷ�λ������ࣨm��
%   center_R0 �����ĵ����б��
%   center_Rc �����Ĳ������Ĵ�Խʱ��б��

% parse the parameters
center_Rc = parameters.center_Rc;   % �����ĵ㲨����Խʱ��б��
theta_rc = parameters.theta_rc_deg * pi / 180;
Nrg = parameters.Nrg;
Naz = parameters.Naz;
Vr = parameters.Vr;
f0 = parameters.f0;
Tr = parameters.Tr;
Kr = parameters.Kr;
BW_dop = parameters.BW_dop;
alpha_os_r = parameters.alpha_os_r;
alpha_os_a = parameters.alpha_os_a;

% targets info
NUM_TARGETS = parameters.NUM_TARGETS;
rs = parameters.rs;
as = parameters.as;

% default parameters
c = 299792458;    % light speed

% derived parameters
Vs = Vr;
Vg = Vr;
lambda = c / f0;
center_R0 = center_Rc * cos(theta_rc);          % �����ĵ����б��
center_etac = - center_Rc * sin(theta_rc) / Vr; % �����ĵ��Ӧ����Բ������Ĵ�Խʱ��;
La = 0.886 * 2 * Vs * cos(theta_rc) / BW_dop;   % ���߿׾�����
beta_bw = 0.886 * lambda / La;                  % we suppose
Fr = Tr * Kr * alpha_os_r;  % ���������
Fa = BW_dop * alpha_os_a;   % �����ظ�Ƶ��
delta_r = c/2/Fr;           % �����������ࣨSAR�źſռ䣩
delta_a = Vr / Fa;          % ��λ��������
% not used by generating the signal, and they will be returned directly for
% the user
f_etac = 2 * Vr * sin(theta_rc) / lambda; % ����������Ƶ��

% �۲�ʱ��������ȷ��
tau0 = 2 * center_Rc / c;
eta0 = center_etac;
tau = ((-Nrg / 2) : (Nrg / 2 - 1)) / Fr + tau0;
eta = ((-Naz / 2 : Naz / 2 - 1)) / Fa + eta0;

%##########################################################################
% �ڴ�̽���������õĵ�Ŀ������ȡֵ��Χ
% SARԭʼ�źſռ�ľ�����ͷ�λ��ʱ����
tau_w = tau(end) - tau(1);
eta_h = eta(end) - eta(1);
% �����б�����귶Χ
% ��linspace�����⣬�����Ż����ο�extra/better_generate_time_freq_axis.m
r_w = linspace(-(tau_w/2-Tr/2), (tau_w/2-Tr/2), Nrg) * c / 2 * cos(theta_rc); % -Tp/2����Ϊ������������
% б�ӽ�ʹ�ò�ͬб�ദĿ�������ķ�λλ�ò�ͬ(1)
a_h_top = eta_h / 2 * Vg + r_w * tan(theta_rc);
a_h_bottom = -eta_h / 2 * Vg + r_w * tan(theta_rc);
% % �ؿ�eta�ᣬ�Ա������ܳ����������е㣨���ڲ�����2��֮���������Ϊeta�ỹҪ�㹻��ȷ��ԭʼ�ź��ܱ���ȫ���գ�
% % ��Ϊuser�Ѿ�������Naz�����ǾͲ�Ҫ���Ըı���
% Naz = round((a_h_top(end)-a_h_bottom(1)) / interval_a);
% Naz = ceil(Naz/2)*2;
% % disp(['Naz2:', num2str(Naz)]);
% eta = ((-Naz / 2 : Naz / 2 - 1)) / Fa + eta0;
% б��Ĳ�ͬʹ�ò�ͬб�ദĿ��ʱ��r-a�������켣��λ�򳤶Ȳ�ͬ(2)
a_h_top = a_h_top - (r_w + center_R0) / cos(theta_rc) * beta_bw / 2;
a_h_bottom = a_h_bottom + (r_w + center_R0) / cos(theta_rc) * beta_bw / 2;

% ��ʾ�ܷ��õ�Ŀ���ķ�Χ
figure;
plot(r_w, a_h_top, 'r-');
hold on;
plot(r_w, a_h_bottom, 'r-');
line([r_w(1); r_w(1)], [a_h_bottom(1); a_h_top(1)], 'color', 'r');
line([r_w(end); r_w(end)], [a_h_bottom(end); a_h_top(end)], 'color', 'r');
xlabel('��Ծ�����б����r��m��');ylabel('��λ��m��');
grid on;

if NUM_TARGETS == 0
    % just want to know the coordiante range of target we can set
    suptitle('б��-��λƽ����Ŀ�������ó������귶Χʾ��ͼ');
    s0 = [];
    hold off;
    return;
else
    suptitle('б��-��λƽ����Ŀ������ó���ʾ��ͼ');
    plot(rs, as, '*');
    hold off;
end

%##########################################################################
% �����Ŀ�������б��
R0s = center_R0 + rs;
% �����Ŀ���ľ����������ʱ��
eta_0s = as / Vr;
% �����Ŀ���ľ��Զ���������ʱ��
eta_cs = eta_0s + center_etac - rs * tan(theta_rc) / Vr;

%% 3. �����״�ԭʼ����
% ʱ����������
[tau_mtx, eta_mtx] = meshgrid(tau, eta);

% ����Ŀ���˲ʱб�ࣨ�淽λʱ��仯��
% ע��etaY-eta_0s(i)ΪĿ�����Ը����������ʱ�̵ķ�λ��ʱ��eta
R_eta = zeros(NUM_TARGETS, Naz, Nrg);
for i = 1:NUM_TARGETS
    R_eta(i, :, :) = sqrt((R0s(i)^2 + Vr^2 * (eta_mtx - eta_0s(i)).^2 ));
end

A0 = 1;
s0 = zeros(Naz, Nrg);
for i = 1:NUM_TARGETS
    % ����
    w_r = (abs(tau_mtx - 2 * reshape(R_eta(i, :, :), Naz, Nrg) / c) <= Tr / 2);
    w_a = sinc(0.886 / beta_bw * atan(Vg * (eta_mtx - eta_cs(i)) / R0s(i))).^2;
    % ��λ
    theta1 = -4 * pi * f0 * reshape(R_eta(i, :, :), Naz, Nrg) / c;
    theta2 = pi * Kr * (tau_mtx - 2 * reshape(R_eta(i, :, :), Naz, Nrg) / c).^2;
    % �źŶ���ۼ�
    s0 = s0 + A0 * w_r .* w_a .* exp(1j*theta1) .* exp(1j*theta2);
end
end

