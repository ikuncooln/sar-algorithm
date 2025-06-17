% Test for generate_point_data() function

close all;clear;



%% 1. 仿真参数 (参考 p142, table 6.1)

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

theta_rc_deg = 21.9; % 低斜视角21.9度

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

parameters = struct('center_Rc', center_Rc, 'theta_rc_deg', theta_rc_deg, 'Nrg', Nrg, 'Naz', Naz, 'Vr', Vr, 'f0', f0, 'Tr', Tr, 'Kr', Kr, 'BW_dop', BW_dop, 'alpha_os_r', alpha_os_r, 'alpha_os_a', alpha_os_a, 'NUM_TARGETS', NUM_TARGETS, 'rs', rs, 'as', as);



[ s0, f_etac, delta_r, delta_a, center_R0, center_Rc ] = generate_point_data(parameters);

figure; % 绘制雷达原始仿真信号

imagesc(abs(s0));xlabel('距离向时间（采样点）');ylabel('方位向时间（采样点）');

title('点目标雷达原始仿真信号幅度（时域）');



%% 3. 处理

s = RDA(s0, lambda, Kr, Vr, Fr, Fa, center_R0, theta_rc_deg, Tr);

x = ((-Nrg / 2) : (Nrg / 2 - 1)) / Fr * c / 2;

y = ((-Naz / 2 : Naz / 2 - 1)) / Fa * Vg ;

figure; % 绘制低斜视角情况下距离压缩且方位压缩后信号

subplot(121);imagesc(real(s));xlabel('距离向时间（采样点）');ylabel('方位向时间（采样点）');title('(a)实部');

subplot(122);imagesc(x, y, abs(s));xlabel('距离向时间（采样点）');title('(b)幅度');set(gca, 'YDir', 'normal');

sgtitle('21.9度斜视角距离压缩且方位压缩后的信号（时域）');



%% 4. 分析

