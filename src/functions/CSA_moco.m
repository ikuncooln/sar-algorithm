function img_cs = CSA_moco(s0,theta_bw,lambda,Kr,Tr,Fr,theta_rc,Nrg,Naz,near_range,Vr,PRF,flag,...
    moco_file, ref_range, subaperture_num, range_size, azimuth_size, range_start, last_pulse_count)
%   Chirp Scaling算法，带有二阶运动补偿（距离空变相位）和方位子孔径处理（方位空变相位）
%   s0 输入信号数据（二维复数组）
%   theta_bw 天线方位向波束宽度(rad)
%   lambda 波长(m)
%   Kr 距离向调频率(Hz/s)
%   Tr 发射脉冲宽度(s)
%   Fr 距离向采样率(Hz)
%   theta_rc 斜视角(rad)
%   Nrg 距离向采样点数
%   Naz 方位向采样点数
%   near_range 第一个采样点距离(m)
%   Vr 载机速度(m/s)
%   PRF 脉冲重复频率(Hz)
%   flag 为1/0表示画/不画出中间步骤的图
%   moco_file 运动补偿元数据文件名（即平台轨迹信息）
%   ref_range 一阶相位运动补偿时选择的参考距离
%   subaperture_num 方位向空变补偿时子孔径数量
%   range_size 每个回波信号总的采样点数
%   azimuth_size 用于理想轨迹拟合的总的回波个数
%   range_start 起始距离向采样点数
%   last_pulse_count 上次已经处理过的脉冲个数（用于追加式处理）
%% 其他参数
c=299792458;
% 距离向
delta_r = c/2/Fr;
center_Rc = near_range + Nrg/2*delta_r;
R_ref = center_Rc*cos(theta_rc);
BW_r = abs(Kr)*Tr;
gama_wr = 1.18;
pr = 0.886*gama_wr/BW_r;
start = near_range*2/c;
% Nr = ceil(Fr*Tr/2)*2;
% 方位向
f0 = c/lambda;
Fa = PRF;
eta_c = (-R_ref*tan(theta_rc))/Vr;
f_etac = 2*Vr*sin(theta_rc)/lambda;
delta_a = Vr/Fa;
gama_wa = 1.185;
La = 0.886*lambda/theta_bw;
pa = La/2*gama_wa;
D_ref = cos(theta_rc);
Vg = Vr;
% Ta = 0.886*lambda*center_Rc/(La*Vg*cos(theta_rc));
% Na = ceil(Fa*Ta/2)*2;
%% 距离多普勒域 变标
if(flag == 1)
    figure;subplot(221);
    imagesc(real(s0));
    xlabel('距离向（采样点）');ylabel('方位向（采样点）');title('(a)原始信号实部');
end
Srd = fft(s0,Naz,1);
clear s0;
if(flag == 1)
    subplot(222);
    imagesc(abs(Srd));
    xlabel('距离向（采样点）');ylabel('方位频率（采样点）');title('(b)原始信号的距离多普勒域');
end
f_eta_mtx = (0:Fa/Naz:(Naz-1)*Fa/Naz).';                   
f_eta_mtx = (round((f_etac-f_eta_mtx)/Fa)*Fa+f_eta_mtx); 
D = sqrt(1-c^2*f_eta_mtx.^2/(4*Vr^2*f0^2));                         % 徙动参数
Km = Kr./(1-Kr*c*R_ref.*f_eta_mtx.^2./(2*Vr^2*f0^3*D.^3));             % 假设距离多普勒域的Km不随距离改变
clear f_eta_mtx;
% ssc_ft = exp(1j*pi*Km.*(D_ref./D_mtx-1).*tr_mtx_new.^2);             % 变标方程
ssc = exp(1j*pi*repmat(Km.*(D_ref./D-1),1,Nrg).*(repmat(start +(0:(Nrg-1))/Fr,Naz,1) - repmat(2*R_ref./(c*D),1,Nrg)).^2);             % 变标方程
clear Km D;
S1 = ssc.*Srd;           % 与变标方程相乘
clear Srd ssc;
%% 二维频域 距离处理
S2 = fft(S1,Nrg,2);         % 距离向FFT
clear S1;
f_eta_mtx = (0:Fa/Naz:(Naz-1)*Fa/Naz).';                   
f_eta_mtx = (round((f_etac-f_eta_mtx)/Fa)*Fa+f_eta_mtx); 
D = sqrt(1-c^2*f_eta_mtx.^2/(4*Vr^2*f0^2));                   % 徙动参数
Km = Kr./(1-Kr*c*R_ref.*f_eta_mtx.^2./(2*Vr^2*f0^3*D.^3));             % 假设距离多普勒域的Km不随距离改变
clear f_eta_mtx;
f_tau_mtx = ifftshift((-Nrg/2:Nrg/2-1)*Fr/Nrg);
% H_range_bulk = exp(1j*pi*(D_mtx.*f_tau_mtx.^2./(Km.*D_ref)+(4*R_ref*f_tau_mtx.*(1./D-1/D_ref)/c)));  % 滤波器，最占内存的1步 
H_range_bulk = exp(1j*pi*(repmat(D./(Km.*D_ref),1,Nrg).*repmat(f_tau_mtx.^2,Naz,1)+(4*R_ref*repmat(f_tau_mtx,Naz,1).*repmat((1./D-1/D_ref)/c,1,Nrg))));  % 滤波器
clear Km D f_tau_mtx;
N_BW_r = round(BW_r/Fr*Nrg);            % Kr*Tr包含的点数
window_r = ifftshift(kaiser(N_BW_r,2.5).');
S3 = S2.*H_range_bulk;
clear S2 H_range_bulk
S3 = S3.*repmat([window_r(1:ceil(N_BW_r/2)),zeros(1,Nrg-N_BW_r),window_r(ceil(N_BW_r/2)+1:N_BW_r)],Naz,1);           % 滤波   
clear window_r;

%% 二维时域 距离空变相位补偿
s4 = ifft2(S3);
clear S3;
% delta_R 一个矩阵，每一行代表1到Nrg这些距离点相对于参考点的距离误差
delta_R = range_space_variant( moco_file,...
    near_range, ref_range, Fr, range_size, azimuth_size,... 
    Nrg, range_start, Naz, last_pulse_count);
s4 = s4 .* exp(1j*4*pi*delta_R/lambda);
clear delta_R;

%% 距离多普勒域 方位处理
S4 = fft(s4);
clear s4;
S4 = azimuth_space_variant( moco_file, S4, lambda, f_etac, Vr,...
    subaperture_num, near_range, Fr, PRF, azimuth_size, last_pulse_count);
%%
%下面需要距离多普勒域的S4
if(flag == 1)
    subplot(223);
    imagesc(abs(S4));
    xlabel('距离向（采样点）');ylabel('方位频率（采样点）');title('(c)距离处理后的距离多普勒域');
end
R0_mtx = (c/2)*(start +(0:(Nrg-1))/Fr)*cos(theta_rc);  % 随距离线变化的最近距离
f_eta_mtx = (0:Fa/Naz:(Naz-1)*Fa/Naz).';                   
f_eta_mtx = (round((f_etac-f_eta_mtx)/Fa)*Fa+f_eta_mtx); 
D = sqrt(1-c^2*f_eta_mtx.^2/(4*Vr^2*f0^2));                   % 徙动参数 
Km_mtx = Kr./(1-Kr*c*repmat(R0_mtx,Naz,1).*repmat(f_eta_mtx.^2./(2*Vr^2*f0^3*D.^3),1,Nrg));    % 随距离变化的距离向调频率
% H_az_add_offset = exp(1j*4*pi*(f0*R0_mtx.*D_mtx/c-Km_mtx.*(1-D_mtx./D_ref).*(R0_mtx./D_mtx-R0./D_mtx).^2./c^2)-1j*2*pi*fa_dop*nc);
S4 = S4.*exp(1j*4*pi*(f0*repmat(R0_mtx,Naz,1).*repmat(D,1,Nrg)/c));
S5 = S4.*exp(1j*4*pi*(-Km_mtx.*repmat((1-D./D_ref),1,Nrg)...
    .*(repmat(R0_mtx,Naz,1)./repmat(D,1,Nrg)-repmat(R_ref./D,1,Nrg)).^2./c^2));
clear S4 Km_mtx R0_mtx D
S5 = S5.*exp(-1j*2*pi*repmat(f_eta_mtx,1,Nrg)*eta_c);
% H_az_add_offset = exp(1j*4*pi*(-Km_mtx.*repmat((1-D./D_ref),1,Nrg)...
%     .*(repmat(R0_mtx,Naz,1)./repmat(D,1,Nrg)-repmat(R_ref./D,1,Nrg)).^2./c^2)-1j*2*pi*repmat(f_eta_mtx,1,Nrg)*eta_c);   % 滤波器，最占内存的1步
clear f_eta_mtx 
% S5 = S4.*H_az_add_offset;                  % 滤波
img_cs = ifft(S5,Naz,1);                % 方位向IFFT
if(flag == 1)
    subplot(224);
    imagesc(abs(img_cs));
    xlabel('距离向（采样点）');ylabel('方位向（采样点）');title('(d)CSA成像结果');
end
end