close all;clear all;
%% 1. 仿真参数 (参考 p142, table 6.1)
mode = 1; % 0: 大斜视角参数，1：小斜视角参数

center_Rc = 20e3;  % 景中心斜距
Vr = 150;   % 等效雷达速度
Tr = 2.5e-6;    % 发射脉冲时宽
Kr = 20e12; % 距离调频率
f0 = 5.3e9; % 雷达工作频率
BW_dop = 80;    % 多普勒带宽
Fr = 60e6;  % 距离采样率
Fa = 100;   % 方位采样率
Naz = 256;  % 方位向采样点数（距离线条数）
Nrg = 256;  % 距离向采样点数（距离线采样点数）
if mode == 0
    theta_rc_deg = 21.9; % 大斜视角21.9度
else
    theta_rc_deg = 3.5; % 小斜视角21.9度
end
c = 299792458;    % 光速

% derived params
lambda = c / f0;
theta_rc = theta_rc_deg * pi / 180;
Vs = Vr;
Vg = Vr;
Np = Tr * Fr;   % 脉冲序列长度（采样点数）
alpha_os_r = Fr / (Kr*Tr);
alpha_os_a = Fa / BW_dop;

%% 2. 生成原始雷达数据
NUM_TARGETS = 3;    % 仿真的目标数为3
rs = [0, 0, 30];    % 各目标距离向距离
as = [-20, 0, -10]; % 目标相对方位向距离
parameters = struct(...
    'center_Rc', center_Rc,...          % 景中心斜距
    'theta_rc_deg', theta_rc_deg,...    % 斜视角
    'Nrg', Nrg,...                      % 距离向采样点数
    'Naz', Naz,...                      % 方位向采样点数
    'Vr', Vr,...                        % 载机速度
    'f0', f0,...                        % 载波频率
    'Tr', Tr,...                        % 发射脉冲宽度
    'Kr', Kr,...                        % 发射脉冲调频率
    'BW_dop', BW_dop,...                % 多普勒带宽
    'alpha_os_r', alpha_os_r,...        % 距离向过采样率
    'alpha_os_a', alpha_os_a,...        % 方位向过采样率
    'NUM_TARGETS', NUM_TARGETS,...      % 点目标数量
    'rs', rs,...                        % 点目标距离向坐标（m）
    'as', as...                         % 点目标方位向坐标（m）
);

[ s0, f_etac, delta_r, delta_a, center_R0, center_Rc ] = generate_point_data(parameters);

figure; % 绘制大斜视角情况下的三点雷达原始仿真信号
subplot(221);imagesc(real(s0));ylabel('方位向时间（采样点）');title('(a)实部');
subplot(222);imagesc(imag(s0));title('(b)虚部');
subplot(223);imagesc(abs(s0));xlabel('距离向时间（采样点）');ylabel('方位向时间（采样点）');title('(c)幅度');
subplot(224);imagesc(angle(s0));xlabel('距离向时间（采样点）');title('(d)相位');
suptitle([num2str(theta_rc_deg), '度斜视角情况下的', num2str(NUM_TARGETS),'点雷达原始仿真信号（时域）']);

%% 3. 成像处理
% step1: do something

% step2: try something

% step3: complete

%% 4. 点目标分析



