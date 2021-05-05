function [img_rd] = RDA_smallSquintAngle(s0, lambda, Kr, Vr, Fr, PRF, theta_rc_deg ,near_range,Tr)
%   Range Doppler Algorithm
% input:
%   s0,基带回波数据，每一列为等距离门。
%   lambda,波长。
%   Kr,调频率。
%   Vr,载机速度。
%   Fr,距离向采样频率。
%   PRF,脉冲重复频率。
%   near_range,第一个采样点距离。
%   theta_rc_deg,斜视角，单位为角度（°）。
% output:
%   img_rd: RDA聚焦后的复图像。

%% 1. 仿真参数 (参考 p142, table 6.1)


Fa = PRF;   % 方位采样率
theta_rc = theta_rc_deg*pi/180;
[Naz,Nrg] = size(s0);
c = 299792458;
f_etac = 2 * Vr * sin(theta_rc) / lambda;   % 多普勒中心频率
delta_r = c/2/Fr;    % 距离向采样间距
center_R0=near_range+Nrg/2/Fr*c/2;
tau0 = 2 * center_R0 / cos(theta_rc) / c;
tau = ((-Nrg / 2) : (Nrg / 2 - 1)) / Fr + tau0;
eta0 = - center_R0 / cos(theta_rc) * sin(theta_rc) / Vr;  % 景中心点对应的相对波束中心穿越时刻;

%% 3. 距离压缩

% 方式3：根据脉冲频谱特性在频域生成滤波器
f_tau = ifftshift((-Nrg/2:Nrg/2-1) * Fr / Nrg); % 生成距离向频率轴
Hrc3_tmp = exp(1j * pi .* f_tau.^2 / Kr);
W3 = abs(f_tau) <= (abs(Kr)*Tr/2); % 这里加这个窗限制的目的在于：（复）信号带宽并不一定等于采样带宽
Hrc3 = W3 .* Hrc3_tmp;    % 滤波器频谱的有效值分布宽度也应等于信号带宽，而不仅仅由f的取值宽度（Fs）决定
Hrc3 = repmat(Hrc3, Naz, 1);

% 在频域匹配滤波
S0 = fft(s0.').';
Src = S0 .* Hrc3 .* repmat(ifftshift(kaiser(Nrg, 2.5).'), Naz, 1);   % 选择方式1或2或3进行距离压缩
s_rc = ifft(Src.').';

%% 4. 方位向傅里叶变换
Srd = fft(s_rc);
% figure; % 绘制距离多普勒域里的距离压缩后的结果
% subplot(121);imagesc(real(Srd));xlabel('距离向时间（采样点）');ylabel('方位向频率（采样点）');title('(a)实部');set(gca, 'YDir', 'normal');
% subplot(122);imagesc(abs(Srd));xlabel('距离向时间（采样点）');title('(b)幅度');set(gca, 'YDir', 'normal');
% suptitle('3.5度斜视角距离压缩后信号（距离多普勒域）');

%% 5. 距离徙动校正
f_eta = (ifftshift((-Naz/2 : Naz/2-1) * Fa / Naz)).';
f_eta = f_eta + round((f_etac - f_eta) / Fa) * Fa;
R0 = tau * c / 2 * cos(theta_rc);
[R0_grid, f_eta_grid] = meshgrid(R0, f_eta);
% 计算距离徙动量矩阵
RCM = lambda^2 * R0_grid .* f_eta_grid.^2 / 8 / Vr^2;
% 此处要非常小心：因为校正时认为Srcmc和Srd矩阵元素是一一对应，但它们所表示的
% 距离线天然的就存在距离差，所以后面校正时要考虑到这一部分已经偏移了，不要重复
% 例：尽管Srcmc(1,:)处对应距离徙动为0，也不能直接Srcmc(1,:)=Srd(1,:)，这是因为
% Srd(1,:)处对应的距离为：Srcmc(1,:)对应的斜距+R_etac-R0c。因此有如下处理：
RCM = RCM - (1/cos(theta_rc)-1)*R0_grid;
RCM = RCM / delta_r;     % 将距离徙动量转换为距离单元偏移量

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
            Srcmc(i,j) = Srd(i,ceil(mod(j+offset_int-0.1,Nrg)));   % 利用信号数据S1的周期性假定
        else
            Srcmc(i,j) = Srd(i, ceil(mod((j+offset_int-4:j+offset_int+3)-0.1,Nrg))) * hx(offset_frac,:).';
        end
        
    end
end

% figure; % 绘制距离多普勒域里的距离压徙动校正后的结果
% subplot(121);imagesc(real(Srcmc));xlabel('距离向时间（采样点）');ylabel('方位向频率（采样点）');title('(a)实部');set(gca, 'YDir', 'normal');
% subplot(122);imagesc(abs(Srcmc));xlabel('距离向时间（采样点）');title('(b)幅度');set(gca, 'YDir', 'normal');
% suptitle('3.5度斜视角距离徙动校正后信号（距离多普勒域）');

%% 6. 方位压缩
Ka = 2 * Vr^2 / lambda ./ R0_grid;
Haz = exp(-1j*pi*f_eta_grid.^2./Ka);
Srd_ac = Srcmc .* Haz;

%% 7. 得到时域SAR图像
% 方位向乘以线性相位，以使得轴的中心对应中心目标的零多普勒时刻（有利于最终图像的显示）
% 或者另一种解释：最终的时间轴是eta0为中心，而我们的目标点显示以0为中心，所以我们可以将目标们搬到eta0附近，以便显示
Srd_ac = Srd_ac .* exp(-1j*2*pi*f_eta_grid*eta0);
img_rd = ifft(Srd_ac);

% x = ((-Nrg / 2) : (Nrg / 2 - 1)) / Fr * c / 2;
% y = ((-Naz / 2 : Naz / 2 - 1)) / Fa * Vr ;
% zero_pos = round((0 - eta0) * Fa + Naz / 2 + 1);
% y = y - y(zero_pos);

% figure; % 绘制低斜视角情况下距离压缩且方位压缩后信号
% subplot(121);imagesc(real(s_ac));xlabel('距离向时间（采样点）');ylabel('方位向时间（采样点）');title('(a)实部');
% subplot(122);imagesc(x, y, abs(s_ac));xlabel('距离向时间（采样点）');title('(b)幅度');set(gca, 'YDir', 'normal');
% suptitle('3.5度斜视角距离压缩且方位压缩后的信号（时域）');


end

