%% Chirp Scaling Algorithm
% 直接点运行，跑出的结果对应书上的图
clear;
close all;
clc;
%% 表4.1C波段星载参数
% 距离向参数
R_nc = 850e3;   % 景中心斜距
Height = 800e3; % 高度
Tr = 40e-6;     % 发射脉冲时宽
Kr = -0.5e12;    % 距离脉冲调频率
BW_r = abs(Kr*Tr);   % 信号带宽
Fr = 24e6;      % 距离采样率
width_r = 50e3; % 斜距条带宽度
% 方位向参数
Vr = 7100;      % 等效雷达速度
Vs = Vr*1.06;   % 卫星切向速度
Vg = Vr*0.94;   % 波束覆盖区沿地面的移动速度
f0 = 5.3e9;     % 雷达工作频率
c = 3e8;        % 光速
lambda = c/f0;   % 雷达工作波长
Ls = 4.8e3;     % 合成孔径长度
La = 10;        % 天线长度
BW_dop = 1338;  % 多普勒带宽
Fa = 1700;      % 方位采样率
theta_r_c = 8*pi/180;    %等效斜视角
% R_nc = (1000e3+20e3)/cos(theta_r_c);
% R_nc = 1000e3/cos(theta_r_c);
% 其他参数
Naz = 1024;    % 方位向采样点数
Nrg = 1024;    % 距离向采样点数

A0 = 1;        % 复常数
% 导出参数
theta_sq_c = theta_r_c;
beta_bw = 0.886*lambda/La;                      % 雷达3dB波束
R0 = R_nc*cos(theta_r_c);	                   % 与R_nc相对应的最近斜距，记为R0
% R_ref = R0-20e3; % 参考距离
R_ref = R0;
Ta = 0.886*lambda*R_nc/(La*Vg*cos(theta_r_c));  % 目标照射时间
Ka = 2*Vr^2*cos(theta_r_c)^3/(lambda*R0);       % 方位向调频率
f_nc = 2*Vr*sin(theta_r_c)/lambda;              % 多普勒中心频率
nc = (-R0*tan(theta_r_c))/Vr;                  % 景中心的波束中心穿越时刻
Nr = round(Tr*Fr);                             % 线性调频信号采样点数
Na = round(Ta*Fa);                             % 方位向采样点数?
a_os_r = Fr/BW_r;                              % 距离向过采样因子
a_os_a = Fa/BW_dop;                            % 方位向过采样因子
Mamb = floor(f_nc/Fa);                         % 多普勒模糊

% Nrg = ceil(Tr*Fr);
% Naz = ceil(Ta*Fa);
%% 目标点相关参数
delta_r1 = -1500;         % 目标1和景中心的距离向距离差
delta_a1 = delta_r1*tan(theta_r_c);            % 目标1和景中心的方位向距离差
% delta_r2 = 0;            % 目标2和景中心的距离向距离差     
delta_r2 = -20e3;            % 目标2和景中心的距离向距离差    
% delta_a2 = 0;            % 目标2和景中心的方位向距离差
delta_a2 = delta_r2*tan(theta_r_c);    
delta_r3 = 1500;          % 目标3和景中心的距离向距离差     
delta_a3 = delta_r3*tan(theta_r_c);            % 目标3和景中心的方位向距离差    

% 目标1
x1 = R0 + delta_r1;      % 目标1的距离向距离
y1 = delta_a1;           % 目标1的方位向距离
% 目标2
x2 = R0 + delta_r2;      % 目标2的距离向距离
y2 = delta_a2;           % 目标2的方位向距离
% 目标3
x3 = R0 + delta_r3;      % 目标3的距离向距离
y3 = delta_a3;           % 目标3的方位向距离
% 构造向量
x_range = [x1,x2,x3];
y_azimuth = [y1,y2,y3];
% 计算波束中心穿越时刻
nc_1 = (y1-x1*tan(theta_r_c))/Vr;           % 目标1的波束中心穿越时刻。
nc_2 = (y2-x2*tan(theta_r_c))/Vr;           % 目标2的波束中心穿越时刻。
nc_3 = (y3-x3*tan(theta_r_c))/Vr;           % 目标3的波束中心穿越时刻。
nc_target = [nc_1,nc_2,nc_3];               % 定义该数组，便于处理。
%% 定义距离向和方位向的时间轴和频率轴
% tr = 2*R_nc/c+(-Nrg/2:(Nrg/2-1))/Fr;           % 距离时间轴，中心为2*R_nc/c
tr = 2*x2/cos(theta_r_c)/c+(-Nrg/2:(Nrg/2-1))/Fr;
ta = nc+(-Naz/2:Naz/2-1)/Fa;                   % 方位时间轴，将y=0（景中心）的零多普勒时刻为方位向时间零点 
% ta = nc_2+(-Naz/2:Naz/2-1)/Fa;                   % 方位时间轴，将y=0（景中心）的零多普勒时刻为方位向时间零点 
% 生成距离（方位）时间（频率）矩阵
tr_mtx = ones(Naz,1)*tr;                       % 距离时间轴矩阵，大小：Naz*Nrg
ta_mtx = ta.'*ones(1,Nrg);                     % 方位时间轴矩阵，大小：Naz*Nrg
% 距离频率轴
fr_mtx = 0:Fr/Nrg:(Nrg-1)*Fr/Nrg;   % 距离向已被解调到基带              
fr_mtx = floor((Fr/2-fr_mtx)/Fr)*Fr+fr_mtx; 
fr_mtx = ones(Naz,1)*fr_mtx;
% 方位频率轴
fa_mtx = 0:Fa/Naz:(Naz-1)*Fa/Naz;     % 要把这些频率值+整数倍的Fa，映射到(f_nc-Fa/2,f_nc+Fa/2]区间里              
fa_mtx = floor((f_nc+Fa/2-fa_mtx)/Fa)*Fa+fa_mtx; 
fa_mtx = fa_mtx.'*ones(1,Nrg);
%% 雷达原始数据
s0_tt = zeros(Naz,Nrg);                                                                   % 时域信号
% for k = 1:3
for k = 2
    R_n = sqrt(x_range(k)^2+(Vr*ta_mtx-y_azimuth(k)).^2);                                 % 瞬时斜距
    wr = (abs(tr_mtx-2*R_n/c)<=Tr/2);                                                     % 距离包络
    wa = sinc(0.886*(atan(Vg*(ta_mtx-nc_target(k))/x_range(k)))/beta_bw).^2;              % 方位包络
    s0_tt = s0_tt + A0*wr.*wa.*exp(-1j*4*pi*f0*R_n/c).*exp(1j*pi*Kr*(tr_mtx-2*R_n/c).^2); % 叠加信号
end
%% 距离多普勒域  
Srd_ft = fft(s0_tt,Naz,1);  % 方位向FFT
figure;
subplot(3,2,[1,3,5]);
imagesc(real(Srd_ft));
set(gca,'YDir','normal');
title('(a)原始信号的距离多普勒频谱');
xlabel('距离（采样点）');
ylabel('方位频率（采样点）');
sample = [300,600,800];  % 选取3个方位频率点观察，具体位置随参数变化，要根据距离多普勒域观察
subplot(322);
plot(real(Srd_ft(sample(1),:)));
xlim([1,Nrg]);
ylabel('归一化幅度');
title(['(b)第',num2str(sample(1)),'个方位频率采样']);
subplot(324);
plot(real(Srd_ft(sample(2),:)));
xlim([1,Nrg]);
ylabel('归一化幅度');
title(['(c)第',num2str(sample(2)),'个方位频率采样']);
subplot(326);
plot(real(Srd_ft(sample(3),:)));
xlim([1,Nrg]);
xlabel('距离（采样点）');
ylabel('归一化幅度');
title(['(d)第',num2str(sample(3)),'个方位频率采样']);
%% 一致和补余距离徙动的幅度
D_mtx = sqrt(1-c^2*fa_mtx.^2/(4*Vr^2*f0^2));     % 徙动参数
tr_RCMC = tr*cos(theta_r_c);                     % 在新的距离线长度下的时间轴，注意中心值发生了变化
R0_mtx = ones(Naz,1)*(c/2)*tr_RCMC;              % 随距离线变化的最近距离
% D_ref = D_mtx(ceil(mod(round(mod(f_nc,Fa)*Naz/Fa)-0.1,Naz)),1);    % 参考方位频率（多普勒中心频率）处的徙动参数
D_ref = cos(theta_r_c);
RCM_bulk = 2*(R_ref./D_mtx-R_ref./D_ref)*Fr/c;                                     % 一致距离徙动
RCM_diff = 2*(R0_mtx./D_mtx - R0_mtx./D_ref - R_ref./D_mtx + R_ref./D_ref)*Fr/c;   % 补余距离徙动
% fa_mtx = floor((f_nc+Fa/2-fa_mtx)/Fa)*Fa+fa_mtx; 
dif = diff(fa_mtx(:,1));
p_break = find(abs(dif)==max(abs(dif)));
sample_r = Nrg/2+1;
figure;
plot(RCM_bulk(1:p_break,sample_r),1:p_break,'r',RCM_diff(1:p_break,sample_r),1:p_break,'b--',...
    RCM_bulk(p_break+1:Naz,sample_r),p_break+1:Naz,'r',RCM_diff(p_break+1:Naz,sample_r),p_break+1:Naz,'b--');
legend('一致距离徙动','补余距离徙动','location','Bestoutside');
ylim([1,Naz]);
set(gca,'YDir','reverse');
title('一致和补余距离徙动的幅度');
xlabel('距离徙动（采样点）');
ylabel('方位频率（采样点）');
%% 只进行距离压缩和一致距离徙动校正
s2_ff = fft(Srd_ft,Nrg,2);                       % 距离向FFT，变换到二维频域
D_mtx = sqrt(1-c^2*fa_mtx.^2/(4*Vr^2*f0^2));     % 徙动参数
Km = Kr./(1-Kr*c*R_ref.*fa_mtx.^2./(2*Vr^2*f0^3*D_mtx.^3));       % 随距离变化的距离向调频率
% D_ref = D_mtx(ceil(mod(round(mod(f_nc,Fa)*Naz/Fa)-0.1,Naz)),1);    % 参考方位频率（多普勒中心频率）处的徙动参数
D_ref = cos(theta_r_c);
H_range = exp(1j*pi*D_mtx.*fr_mtx.^2./(Km.*D_ref));            % 消除距离调制项
H_RCM_bulk = exp(1j*4*pi*R_ref*fr_mtx.*(1./D_mtx-1/D_ref)/c);         % 一致距离徙动校正
N_BW_r = round(Nrg/a_os_r);    
window = ifftshift(kaiser(N_BW_r,2.5)');    % Kaiser窗
window = (ones(Naz,1)*[window(1:ceil(N_BW_r/2)),zeros(1,Nrg-N_BW_r),window(ceil(N_BW_r/2)+1:N_BW_r)]);
s3_ff = s2_ff.*H_range.*H_RCM_bulk.*window;    % 滤波     
s4_ft = ifft(s3_ff,Nrg,2);                     % 距离向IFFT
N_cut = 32;
figure;
subplot(3,2,[1,3,5]);
imagesc(abs(s4_ft));
xlim([Nrg/2+1-N_cut/2,Nrg/2+N_cut/2]);
set(gca,'YDir','normal');
title('(a)距离压缩后的距离多普勒频谱');
xlabel('距离（采样点）');
ylabel('方位频率（采样点）');
% 升采样分析
Num = N_cut*8;      % 升采样点数
% s4_ff = fft(s4_ft,Nrg,2);
% s4_ft_up = ifft([s4_ff(:,1:Nrg/2),zeros(Naz,Num-Nrg),s4_ff(:,Nrg/2+1:Nrg)],Num,2);
% figure;imagesc(abs(s4_ft_up));
% set(gca,'YDir','normal');
% title('升采样-距离压缩与一致RCMC后');
% xlabel('距离时间（采样点）');
% ylabel('方位频率（采样点）');
% 选取第20、60和240个方位频率采样点的距离线升采样分析
line1_f = fft(s4_ft(sample(1),Nrg/2+1-N_cut/2:Nrg/2+N_cut/2)); 
line1_up = ifft([line1_f(1:N_cut/2),zeros(1,Num-N_cut),line1_f(N_cut/2+1:N_cut)]);
line2_f = fft(s4_ft(sample(2),Nrg/2+1-N_cut/2:Nrg/2+N_cut/2));
line2_up = ifft([line2_f(1:N_cut/2),zeros(1,Num-N_cut),line2_f(N_cut/2+1:N_cut)]);
line3_f = fft(s4_ft(sample(3),Nrg/2+1-N_cut/2:Nrg/2+N_cut/2));
line3_up = ifft([line3_f(1:N_cut/2),zeros(1,Num-N_cut),line3_f(N_cut/2+1:N_cut)]);
subplot(322);
plot(0:N_cut/Num:N_cut-N_cut/Num,abs(line1_up));
xlim([0,N_cut]);
[M,I] = max(abs(line1_up));
text((I-1)*N_cut/Num+4,M,['峰值位于第',num2str((I-1)*N_cut/Num),'个采样处']);
ylabel('幅度');
title(['(b)第',num2str(sample(1)),'个方位频率采样']);
subplot(324);
plot(0:N_cut/Num:N_cut-N_cut/Num,abs(line2_up));
xlim([0,N_cut]);
[M,I] = max(abs(line2_up));
text((I-1)*N_cut/Num+4,M,['峰值位于第',num2str((I-1)*N_cut/Num),'个采样处']);
ylabel('幅度');
title(['(c)第',num2str(sample(2)),'个方位频率采样']);
subplot(326);
plot(0:N_cut/Num:N_cut-N_cut/Num,abs(line3_up));
xlim([0,N_cut]);
[M,I] = max(abs(line3_up));
text((I-1)*N_cut/Num+4,M,['峰值位于第',num2str((I-1)*N_cut/Num),'个采样处']);
xlabel('距离（采样点）');
ylabel('幅度');
title(['(d)第',num2str(sample(3)),'个方位频率采样']);
%% 进行变标、距离压缩和一致距离徙动校正
% 变标方程
D_mtx = sqrt(1-c^2*fa_mtx.^2/(4*Vr^2*f0^2));                       % 徙动参数
% D_ref = D_mtx(ceil(mod(round(mod(f_nc,Fa)*Naz/Fa)-0.1,Naz)),1);    % 参考方位频率（多普勒中心频率）处的徙动参数
D_ref = cos(theta_r_c);
tr_mtx_new = tr_mtx - 2*R_ref./(c*D_mtx);                             % 新的距离时间
Km = Kr./(1-Kr*c*R_ref.*fa_mtx.^2./(2*Vr^2*f0^3*D_mtx.^3));           % 假设距离多普勒域的Km不随距离改变
ssc_ft = exp(1j*pi*Km.*(D_ref./D_mtx-1).*tr_mtx_new.^2);              % 变标方程
S1_ft = ssc_ft.*Srd_ft;                                            % 与变标方程相乘
s2_ff = fft(S1_ft,Nrg,2);                                          % 距离向FFT
% 距离压缩和一致距离徙动校正
Km = Kr./(1-Kr*c*R_ref.*fa_mtx.^2./(2*Vr^2*f0^3*D_mtx.^3));       
H_range = exp(1j*pi*D_mtx.*fr_mtx.^2./(Km.*D_ref));            % 消除距离调制项
H_RCM_bulk = exp(1j*4*pi*R_ref*fr_mtx.*(1./D_mtx-1/D_ref)/c);         % 一致距离徙动校正
N_BW_r = round(BW_r/Fr*Nrg);            % Kr*Tr包含的点数
window = ifftshift(kaiser(N_BW_r,2.5)');    % Kaiser窗
window = (ones(Naz,1)*[window(1:ceil(N_BW_r/2)),zeros(1,Nrg-N_BW_r),window(ceil(N_BW_r/2)+1:N_BW_r)]);
% window = ifftshift(ones(Naz,1)*kaiser(Nrg,2.5)');         % 距离平滑窗
s3_ff = s2_ff.*H_range.*H_RCM_bulk.*window;    % 滤波     
s4_ft = ifft(s3_ff,Nrg,2);                     % 距离向IFFT

N_cut = 32;
figure;
subplot(3,2,[1,3,5]);
imagesc(abs(s4_ft));
xlim([Nrg/2+1-N_cut/2,Nrg/2+N_cut/2]);
set(gca,'YDir','normal');
title('(a)距离压缩后的距离多普勒频谱');
xlabel('距离（采样点）');
ylabel('方位频率（采样点）');
% 升采样分析
Num = N_cut*8;      % 升采样点数
% s4_ff = fft(s4_ft,Nrg,2);
% s4_ft_up = ifft([s4_ff(:,1:Nrg/2),zeros(Naz,Num-Nrg),s4_ff(:,Nrg/2+1:Nrg)],Num,2);
% figure;imagesc(abs(s4_ft_up));
% set(gca,'YDir','normal');
% title('升采样-距离压缩与一致RCMC后');
% xlabel('距离时间（采样点）');
% ylabel('方位频率（采样点）');
% 选取第20、60和240个方位频率采样点的距离线升采样分析
line1_f = fft(s4_ft(sample(1),Nrg/2+1-N_cut/2:Nrg/2+N_cut/2)); 
line1_up = ifft([line1_f(1:N_cut/2),zeros(1,Num-N_cut),line1_f(N_cut/2+1:N_cut)]);
line2_f = fft(s4_ft(sample(2),Nrg/2+1-N_cut/2:Nrg/2+N_cut/2));
line2_up = ifft([line2_f(1:N_cut/2),zeros(1,Num-N_cut),line2_f(N_cut/2+1:N_cut)]);
line3_f = fft(s4_ft(sample(3),Nrg/2+1-N_cut/2:Nrg/2+N_cut/2));
line3_up = ifft([line3_f(1:N_cut/2),zeros(1,Num-N_cut),line3_f(N_cut/2+1:N_cut)]);
subplot(322);
plot(0:N_cut/Num:N_cut-N_cut/Num,abs(line1_up));
xlim([0,N_cut]);
[M,I] = max(abs(line1_up));
text((I-1)*N_cut/Num+4,M,['峰值位于第',num2str((I-1)*N_cut/Num),'个采样处']);
ylabel('幅度');
title(['(b)第',num2str(sample(1)),'个方位频率采样']);
subplot(324);
plot(0:N_cut/Num:N_cut-N_cut/Num,abs(line2_up));
xlim([0,N_cut]);
[M,I] = max(abs(line2_up));
text((I-1)*N_cut/Num+4,M,['峰值位于第',num2str((I-1)*N_cut/Num),'个采样处']);
ylabel('幅度');
title(['(c)第',num2str(sample(2)),'个方位频率采样']);
subplot(326);
plot(0:N_cut/Num:N_cut-N_cut/Num,abs(line3_up));
xlim([0,N_cut]);
[M,I] = max(abs(line3_up));
text((I-1)*N_cut/Num+4,M,['峰值位于第',num2str((I-1)*N_cut/Num),'个采样处']);
xlabel('距离（采样点）');
ylabel('幅度');
title(['(d)第',num2str(sample(3)),'个方位频率采样']);
%% 方位处理
tr_RCMC = tr*cos(theta_r_c);                     % 在新的距离线长度下的时间轴，注意中心值发生了变化
% tr_RCMC = 2*R0/c + (-Nrg/2:(Nrg/2-1))/Fr;
% D_ref = D_mtx(ceil(mod(round(mod(f_nc,Fa)*Naz/Fa)-0.1,Naz)),1);
R0_mtx = ones(Naz,1)*(c/2)*tr_RCMC;              % 随距离线变化的最近距离
Km_mtx = Kr./(1-Kr*c*R0_mtx.*fa_mtx.^2./(2*Vr^2*f0^3*D_mtx.^3));   % 随距离变化的距离向调频率
H_az = exp(1j*4*pi*f0*R0_mtx.*D_mtx/c);           % 方位向匹配滤波
H_add = exp(-1j*4*pi*Km_mtx.*(1-D_mtx./D_ref).*(R0_mtx./D_mtx-R_ref./D_mtx).^2./c^2);   % 附加相位校正
H_offset = exp(-1j*2*pi*(fa_mtx*nc));            % 把方位时间中心改为0
% H_offset = exp(-1j*2*pi*(fa_mtx*nc_2));            % 把方位时间中心改为0
% window = kaiser(Naz,2.5)*ones(1,Nrg);         % 方位平滑窗
s5_ft = s4_ft.*H_az.*H_add.*H_offset;          % 滤波
s5_tt = ifft(s5_ft,Naz,1);           % 方位向IFFT
figure;
imagesc(abs(s5_tt));
title('CSA处理结果');
xlabel('距离时间（采样点）');
%% 点目标放大
area = -8:1:7;       % 设置切片的大小
% 目标2
pos_2_a = round(Naz/2 + 1 + y2*Fa/Vr);          % 目标2方位向的位置
% pos_2_a = round(Naz/2 + 1);          % 目标2方位向的位置
% pos_2_r = round(Nrg/2 + 1 + 2*(x2-R0)*Fr/c/cos(theta_r_c));    % 目标2距离向的位置
pos_2_r = round(Nrg/2 + 1);    % 目标2距离向的位置
amplify_2 = s5_tt(ceil(mod(pos_2_a+area-0.1,Naz)),ceil(mod(pos_2_r+area-0.1,Nrg)));
figure;imagesc(abs(amplify_2));
title('点目标2放大（幅度）');
xlabel('距离时间（采样点）');
ylabel('方位时间（采样点）');
target = amplify_2;        % 选择目标进行分析
[image_upsample,signal_r,quality_r,signal_a,quality_a] = f_point_analyse(target,c/2/Fr,Vg/Fa,32,1);
IRW_r_theory = c/2/BW_r*0.886*1.18;
IRW_a_theory = La/2*Vg/Vs*1.185;


% [image_upsample,az_cut,rg_cut,IRW_a,IRW_r,PSLR_a,PSLR_r] = f_point_analyze_trial(target,16,1);
% Len = length(area);
% Num = size(image_upsample,1);
% image_upsample_dB = 20*log10(abs(image_upsample)/max(max(abs(image_upsample))));
% angle_upsample = angle(image_upsample);
% figure;
% subplot(321);imagesc(0:Len/Num:Len-Len/Num,0:Len/Num:Len-Len/Num,abs(image_upsample));title('(a)放大后的点目标');xlabel('距离向（采样点）');ylabel('方位向（采样点）');
% subplot(322);contour(0:Len/Num:Len-Len/Num,0:Len/Num:Len-Len/Num,abs(image_upsample),30);title('(b)放大后的点目标等值线图');set(gca,'YDir','reverse');xlabel('距离向（采样点）');ylabel('方位向（采样点）');
% subplot(323);plot(0:Len/Num:Len-Len/Num,image_upsample_dB(az_cut,:));xlim([0,Len]);ylim([-35,0]);title('(c)距离剖面图');xlabel('距离时间（采样点）');ylabel('幅度(dB)');
% subplot(324);plot(0:Len/Num:Len-Len/Num,image_upsample_dB(:,rg_cut));xlim([0,Len]);ylim([-35,0]);title('(d)方位剖面图');xlabel('方位时间（采样点）');ylabel('幅度(dB)');
% subplot(325);plot(0:Len/Num:Len-Len/Num,angle_upsample(az_cut,:)*180/pi);xlim([0,Len]);title('(e)距离相位');xlabel('距离时间（采样点）');ylabel('相位（°）');
% subplot(326);plot(0:Len/Num:Len-Len/Num,angle_upsample(:,rg_cut)*180/pi);xlim([0,Len]);title('(f)方位相位');xlabel('方位时间（采样点）');ylabel('相位（°）');