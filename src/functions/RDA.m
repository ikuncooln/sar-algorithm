function [ img_rd] = RDA( s0, lambda, Kr, Vr, Fr, PRF,center_R0,theta_rc_deg,flag )
%   Range Doppler Algorithm
% input:
%   s0,基带回波数据，每一列为等距离门。
%   lambda,波长。
%   Kr,调频率。
%   Vr,载机速度。
%   Fr,距离向采样频率。
%   PRF,脉冲重复频率。
%   center_R0,第一个采样点距离。
%   theta_rc_deg,斜视角，单位为角度（°）。
%   flag,是否画出中间步骤图。
% output:
%   img_rd: RDA聚焦后的复图像。

%% 确定其它参数
c = 299792458;
Fa = PRF;                      % 方位向采样率
theta_rc = theta_rc_deg*pi/180;
[Naz,Nrg] = size(s0);
f0=c/lambda;
delta_r = c/2/Fr;              % 距离向采样间距。
f_etac = 2*Vr*sin(theta_rc)/lambda;% 多普勒中心频率

%% 距离压缩(采用方式3匹配滤波）
f_tau = ifftshift((-Nrg/2:Nrg/2-1)*Fr/Nrg); % 距离向频率轴
Hrc = exp(1j*pi*f_tau.^2/Kr);  % Matched filter in Frequency domain
Hrc = repmat(Hrc,Naz,1);
s0_tmp = fft(s0.').';          %距离频域方位时域 fft默认按列
%注意这里不用fftshift
Src = s0_tmp.*Hrc;             %匹配滤波
s_rc = ifft(Src.').';

if flag == 1 
    figure;
    imagesc(abs(s_rc));
    title('距离压缩后图像');
    xlabel('距离向时间（采样点）');
    ylabel('方位向时间（采样点）');
%     colormap('gray');
end 

clear s0  Hrc s0_tmp Src
%% 方位向傅里叶变换
Srd = fft(s_rc);
clear s_rc
%% 二次距离压缩(二维频域进行）
S2df  = fft(Srd.').';          %前面为了观察频域做了傅里叶逆变换，实际可以省略。
clear Srd
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
clear f_tau_mtx f_tau
Ssrc = S2df.*Hsrc;              %二维频域中实现二次压缩
clear  Hsrc   S2df  f_eta
s_src = ifft(Ssrc.').';
clear Ssrc
if flag == 1 
    figure;
    imagesc(abs(s_src));
    title('二次距离压缩后图像');
    xlabel('距离向时间（采样点）');
    ylabel('方位向频域（采样点）');
%     colormap('gray');
end 


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
hwait=waitbar(0,'请等待>>>>>>>>');
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
    if mod(i,200)== 0
    waitbar(i/Naz,hwait,['插值处理中: ', num2str(i/Naz*100), '%']);
    end
end
close(hwait);
 
if flag == 1
    figure;
    imagesc(abs(Srcmc));
    title('距离徙动校正后图像');
    xlabel('距离向时间（采样点）');
    ylabel('方位向频域（采样点）');
%     colormap('gray');
end

clear RCM kwin
%% 方位向压缩
% Srcmc=s_src;%不做距离徙动校正
Haz = exp(1j*4*pi.*R0_grid.*D_grid *f0 /c);% 注意此处方位压缩多补偿了个4*pi*R0*f0/c的相位
Srd_ac = Srcmc.*Haz;
clear Srcmc R0_grid D_grid  Haz
%%
eta0 = -center_R0 / cos(theta_rc)*sin(theta_rc)/Vr; %景中心点对应的相对波束中心穿越时刻；
Srd_ac = Srd_ac.*exp(-1j*2*pi*f_eta_grid*eta0);
clear f_eta_grid 
img_rd = ifft(Srd_ac);
clear Srd_ac
end
