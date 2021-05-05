function img_cs = CSA(s0,theta_bw,lambda,Kr,Tr,Fr,theta_rc,Nrg,Naz,near_range,Vr,PRF,flag)
%   Chirp Scaling算法成像
%   s0是输入信号数据（二维复数组）
%   theta_bw是天线方位向波束宽度(rad)
%   lambda是波长(m)
%   Kr是距离向调频率(Hz/s)
%   Tr是发射脉冲宽度(s)
%   Fr是距离向采样率(Hz)
%   theta_rc是斜视角(rad)
%   Nrg是距离向采样点数
%   Naz是方位向采样点数
%   near_range是第一个采样点距离(m)
%   Vr是载机速度(m/s)
%   PRF是脉冲重复频率(Hz)
%   flag为1/0表示画/不画出中间步骤的图
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
f_eta_mtx = (round((f_etac-f_eta_mtx)/Fa)*Fa+f_eta_mtx)*ones(1,Nrg); 
D = sqrt(1-c^2*f_eta_mtx.^2/(4*Vr^2*f0^2));                         % 徙动参数
Km = Kr./(1-Kr*c*R_ref.*f_eta_mtx.^2./(2*Vr^2*f0^3*D.^3));             % 假设距离多普勒域的Km不随距离改变
clear f_eta_mtx;
tau_mtx_new = ones(Naz,1)*(start +(0:(Nrg-1))/Fr) - 2*R_ref./(c*D);  % 新的距离时间
ssc = exp(1j*pi*Km.*(D_ref./D-1).*tau_mtx_new.^2);             % 变标方程
clear tau_mtx_new Km D;
S1 = ssc.*Srd;           % 与变标方程相乘
clear Srd ssc;
%% 二维频域 距离处理
S2 = fft(S1,Nrg,2);         % 距离向FFT
clear S1;
f_eta_mtx = (0:Fa/Naz:(Naz-1)*Fa/Naz).';                   
f_eta_mtx = (round((f_etac-f_eta_mtx)/Fa)*Fa+f_eta_mtx)*ones(1,Nrg); 
D = sqrt(1-c^2*f_eta_mtx.^2/(4*Vr^2*f0^2));                   % 徙动参数
Km = Kr./(1-Kr*c*R_ref.*f_eta_mtx.^2./(2*Vr^2*f0^3*D.^3));             % 假设距离多普勒域的Km不随距离改变
clear f_eta_mtx;
f_tau_mtx = ones(Naz,1)*ifftshift((-Nrg/2:Nrg/2-1)*Fr/Nrg);
H_range_bulk = exp(1j*pi*(D.*f_tau_mtx.^2./(Km.*D_ref)+(4*R_ref*f_tau_mtx.*(1./D-1/D_ref)/c)));  % 滤波器
clear Km D f_tau_mtx;
N_BW_r = round(BW_r/Fr*Nrg);            % Kr*Tr包含的点数
window_r = ifftshift(kaiser(N_BW_r,2.5).');
window_r = (ones(Naz,1)*[window_r(1:ceil(N_BW_r/2)),zeros(1,Nrg-N_BW_r),window_r(ceil(N_BW_r/2)+1:N_BW_r)]);  % 距离平滑窗
S3 = S2.*H_range_bulk.*window_r;           % 滤波   
clear S2 H_range_bulk window_r;
%% 距离多普勒域 方位处理
S4 = ifft(S3,Nrg,2);       % 距离向IFFT
clear S3;
if(flag == 1)
    subplot(223);
    imagesc(abs(S4));
    xlabel('距离向（采样点）');ylabel('方位频率（采样点）');title('(c)距离处理后的距离多普勒域');
end
R0_mtx = (c/2)*ones(Naz,1)*(start +(0:(Nrg-1))/Fr)*cos(theta_rc);  % 随距离线变化的最近距离
f_eta_mtx = (0:Fa/Naz:(Naz-1)*Fa/Naz).';                   
f_eta_mtx = (round((f_etac-f_eta_mtx)/Fa)*Fa+f_eta_mtx)*ones(1,Nrg); 
D = sqrt(1-c^2*f_eta_mtx.^2/(4*Vr^2*f0^2));                   % 徙动参数 
Km_mtx = Kr./(1-Kr*c*R0_mtx.*f_eta_mtx.^2./(2*Vr^2*f0^3*D.^3));    % 随距离变化的距离向调频率
H_az_add_offset = exp(1j*4*pi*(f0*R0_mtx.*D/c-Km_mtx.*(1-D./D_ref).*(R0_mtx./D-R_ref./D).^2./c^2)-1j*2*pi*f_eta_mtx*eta_c);   % 滤波器，最占内存的1步
clear f_eta_mtx Km_mtx R0_mtx D;
S5 = S4.*H_az_add_offset;                  % 滤波
clear S4 H_az_add_offset;
img_cs = ifft(S5,Naz,1);                % 方位向IFFT
if(flag == 1)
    subplot(224);
    imagesc(abs(img_cs));
    xlabel('距离向（采样点）');ylabel('方位向（采样点）');title('(d)CSA成像结果');
end
end