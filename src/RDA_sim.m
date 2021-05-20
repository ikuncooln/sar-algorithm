close all;clear all;
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
theta_rc_deg = 21.9; % 大斜视角21.9度
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
as = [-20, 0, -7.94]; % 目标相对方位向距离
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


%% 距离压缩(采用方式3匹配滤波）
f_tau = ifftshift((-Nrg/2:Nrg/2-1)*Fr/Nrg); % 距离向频率轴
Hrc = exp(1j*pi*f_tau.^2/Kr);  % Matched filter in Frequency domain
a_os_r = Fr/abs(Kr*Tr);
N_BW_r = round(Nrg/a_os_r);            % Kr*Tr包含的点数
window_r = ifftshift(kaiser(N_BW_r,2.5)');    % Kaiser窗
window_r = repmat([window_r(1:ceil(N_BW_r/2)),zeros(1,Nrg-N_BW_r),window_r(ceil(N_BW_r/2)+1:N_BW_r)],Naz,1);
Hrc = repmat(Hrc,Naz,1);
Hrc = Hrc.*window_r;
s0_tmp = fft(s0.').';          %距离频域方位时域 fft默认按列
%注意这里不用fftshift
Src = s0_tmp.*Hrc;             %匹配滤波
s_rc = ifft(Src.').';

%% 方位向傅里叶变换
Srd = fft(s_rc);
figure; % 绘制距离压缩后的结果
subplot(221);imagesc(real(s_rc));ylabel('方位向时间（采样点）');title('(a)二次距离压缩前实部（时域）');
subplot(222);imagesc(abs(s_rc));title('(b)二次距离压缩前幅度（时域）');
subplot(223);imagesc(real(Srd));xlabel('距离向时间（采样点）');ylabel('方位向频率（采样点）');title('(c)二次距离压缩前实部（距离多普勒域）');set(gca, 'YDir', 'normal');
subplot(224);imagesc(abs(Srd));xlabel('距离向时间（采样点）');title('(d)二次距离压缩前幅度（距离多普勒域）');set(gca, 'YDir', 'normal');
suptitle('21.9度斜视角距离压缩后信号（时域与距离多普勒域）');
%% 二次距离压缩(二维频域进行）
S2df  = fft(Srd.').';          %前面为了观察频域做了傅里叶逆变换，实际可以省略。
f_eta = (ifftshift((-Naz/2:Naz/2-1)*Fa/Naz)).';
f_eta = f_eta + round((f_etac - f_eta)/Fa)*Fa;
tau0 = 2*center_R0/cos(theta_rc)/c;
tau = (-Nrg/2:Nrg/2-1)/Fr+tau0;  % 距离时间轴
R0 = tau*c/2*cos(theta_rc);
[R0_grid,f_eta_grid] = meshgrid(R0,f_eta);

%计算Range-Doppler域中的距离徙动因子
D = sqrt(1-lambda^2.*f_eta.^2/4/Vr^2);
%计算二次压缩调频率
Ksrc = 2*Vr^2*f0^3.*D.^3/c/center_R0./f_eta.^2;
f_tau_mtx = repmat(f_tau,Naz,1);
%计算二次压缩滤波器
Hsrc = exp(-1j*pi*f_tau_mtx.^2./repmat(Ksrc,1,Nrg));
Ssrc = S2df.*Hsrc;              %二维频域中实现二次压缩
s_src = ifft(Ssrc.').';

figure; % 绘制距离多普勒域里距离压缩后和二次距离压缩后的结果
subplot(221);imagesc(real(Srd));ylabel('方位向频率（采样点）');title('(a)二次距离压缩前实部');set(gca, 'YDir', 'normal');
subplot(222);imagesc(abs(Srd));title('(b)二次距离压缩前幅度');set(gca, 'YDir', 'normal');
subplot(223);imagesc(real(s_src));xlabel('距离向时间（采样点）');ylabel('方位向频率（采样点）');title('(c)二次距离压缩后实部');set(gca, 'YDir', 'normal');
subplot(224);imagesc(abs(s_src));xlabel('距离向时间（采样点）');title('(d)二次距离压缩后幅度');set(gca, 'YDir', 'normal');
suptitle('21.9度斜视角二次距离压缩前后信号对比（距离多普勒域）');

%% 距离徙动校正（RCMC)
D_grid = repmat(D,1,Nrg);
RCM = R0_grid./D_grid-R0_grid;
RCM = RCM - (1/cos(theta_rc)-1)*R0_grid;
RCM = RCM / delta_r; %将距离徙动量转换为距离单元偏移量。
% 计算插值核系数表
x_tmp = repmat(-4:3, 16, 1);
offset_tmp = (1:16)/16;
x_tmp = x_tmp + repmat(offset_tmp.', 1, 8);
hx = sinc(x_tmp);
x_tmp16 = x_tmp .* 16;
x_tmp16 = round(x_tmp16 + 16 * 8 / 2);
kwin = repmat(kaiser(16*8, 2.5).', 16, 1);
hx = kwin(x_tmp16) .* hx;
hx = hx ./ repmat(sum(hx, 2), 1, 8);

% 插值校正
Srcmc = zeros(Naz, Nrg);  % 存放距离徙动校正后的回波信号
for i = 1:Naz
    for j = 1:Nrg
        offset_int = ceil(RCM(i,j));
        offset_frac = round((offset_int - RCM(i,j)) * 16);
        if offset_frac == 0
            Srcmc(i,j) = s_src(i,ceil(mod(j+offset_int-0.1,Nrg)));   % 利用信号数据S1的周期性假定
        else
            Srcmc(i,j) = s_src(i, ceil(mod((j+offset_int-4:j+offset_int+3)-0.1,Nrg))) * hx(offset_frac,:).';
        end
        
    end
end
figure; % 绘制距离多普勒域里的距离压徙动校正前的结果
subplot(121);imagesc(real(s_src));xlabel('距离向时间（采样点）');ylabel('方位向频率（采样点）');title('(a)实部');set(gca, 'YDir', 'normal');
subplot(122);imagesc(abs(s_src));xlabel('距离向时间（采样点）');title('(b)幅度');set(gca, 'YDir', 'normal');
suptitle('21.9度斜视角距离徙动校正前信号（距离多普勒域）');

figure; % 绘制距离多普勒域里的距离压徙动校正后的结果
subplot(121);imagesc(real(Srcmc));xlabel('距离向时间（采样点）');ylabel('方位向频率（采样点）');title('(a)实部');set(gca, 'YDir', 'normal');
subplot(122);imagesc(abs(Srcmc));xlabel('距离向时间（采样点）');title('(b)幅度');set(gca, 'YDir', 'normal');
suptitle('21.9度斜视角距离徙动校正后信号（距离多普勒域）');

%% 方位向压缩
% Srcmc=s_src;%不做距离徙动校正
Haz = exp(1j*4*pi.*R0_grid.*D_grid *f0 /c);% 注意此处方位压缩多补偿了个4*pi*R0*f0/c的相位
Srd_ac = Srcmc.*Haz;
%%
eta0 = -center_R0 / cos(theta_rc)*sin(theta_rc)/Vr; %景中心点对应的相对波束中心穿越时刻；
Srd_ac = Srd_ac.*exp(-1j*2*pi*f_eta_grid*eta0);
img_rd = ifft(Srd_ac);

figure; % 绘制低斜视角情况下距离压缩且方位压缩后信号
subplot(121);imagesc(real(img_rd));xlabel('距离向时间（采样点）');ylabel('方位向时间（采样点）');title('(a)实部');
subplot(122);imagesc(abs(img_rd));xlabel('距离向时间（采样点）');title('(b)幅度');
suptitle('21.9度斜视角距离压缩且方位压缩后的信号（时域）');

%% 4. 点目标分析
% 计算每个点出现位置的索引值
delta_r=delta_r * cos(theta_rc);
ns = round(rs/(delta_r)) + (Nrg/2 + 1);
ms = round(as/delta_a) + (Naz/2 + 1);
len = 16;
p = 1;
target = img_rd(ms(p)-len/2:ms(p)+len/2-1, ns(p)-len/2:ns(p)+len/2-1);
[image_upsample,signal_r,quality_r,signal_a,quality_a] = f_point_analyse(target,delta_r,delta_a);
BW_r= abs(Kr*Tr);
La = 0.886 * 2 * Vs * cos(theta_rc) / BW_dop;   % 天线孔径长度
IRW_r_theory = c/2/BW_r*0.886*1.18;
IRW_a_theory = La/2*Vg/Vs*1.185;
disp(['距离向理论分辨率:',num2str(IRW_r_theory),'m']);
disp(['方位向理论分辨率:',num2str(IRW_a_theory),'m']);
