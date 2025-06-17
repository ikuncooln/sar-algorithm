function img_wk = wKA1(s0,theta_bw,lambda,Kr,Tr,Fr,theta_rc,Nrg,Naz,near_range,Vr,PRF,flag)
%   wK算法成像
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
R_nc = near_range + Nrg/2*delta_r;
R_ref = R_nc*cos(theta_rc);
BW_r = abs(Kr)*Tr;
gama_wr = 1.18;
pr = 0.886*gama_wr/BW_r;
start = near_range*2/c;
a_os_r = Fr/BW_r;
% Nr = ceil(Fr*Tr/2)*2;
% 方位向
f0 = c/lambda;
Fa = PRF;
nc = (-R_ref*tan(theta_rc))/Vr;
f_nc = 2*Vr*sin(theta_rc)/lambda;
delta_a = Vr/Fa;
gama_wa = 1.185;
La = 0.886*lambda/theta_bw;
pa = La/2*gama_wa;
D_ref = cos(theta_rc);
Vg = Vr;
% Ta = 0.886*lambda*center_Rc/(La*Vg*cos(theta_rc));
% Na = ceil(Fa*Ta/2)*2;
% 距离频率轴
fr_mtx = ifftshift((-Nrg/2:Nrg/2-1)*Fr/Nrg);
% 方位频率轴
fa_mtx = (0:Fa/Naz:(Naz-1)*Fa/Naz).';     % 要把这些频率值+整数倍的Fa，映射到[f_nc-Fa/2,f_nc+Fa/2)区间里
fa_mtx = round((f_nc-fa_mtx)/Fa)*Fa+fa_mtx;
%% 二维FFT
s0_ff = fft2(s0,Naz,Nrg);  % 二维FFT
clear s0
%% 参考函数相乘（一致压缩）
% H_ref = exp(1j.*(4*pi*R_ref/c*sqrt((f0+repmat(fr_mtx,Naz,1)).^2-c^2.*repmat(fa_mtx,1,Nrg).^2/(4*Vr^2))+pi*repmat(fr_mtx,Naz,1).^2/Kr));
s0_ff = s0_ff.*exp(1j.*(4*pi*R_ref/c*sqrt((f0+repmat(fr_mtx,Naz,1)).^2-c^2.*repmat(fa_mtx,1,Nrg).^2/(4*Vr^2))+pi*repmat(fr_mtx,Naz,1).^2/Kr));
clear H_ref;
N_BW_r = round(Nrg/a_os_r);            % Kr*Tr包含的点数
window_r = ifftshift(kaiser(N_BW_r,2.5)');    % Kaiser窗
window_r = repmat([window_r(1:ceil(N_BW_r/2)),zeros(1,Nrg-N_BW_r),window_r(ceil(N_BW_r/2)+1:N_BW_r)],Naz,1);
s1_ff = s0_ff.*window_r;       % 参考函数相乘
clear s0_ff window_r
%% Stolt映射
% H_shift = exp(-1j*2*pi*(repmat(fr_mtx,Naz,1)*(-Nrg/2/Fr)+fa_mtx*(-Naz/2/Fa)));  % 将距离和方位时间轴的零点由中心挪到开始
H_shift = exp(-1j*2*pi*(repmat(fr_mtx,Naz,1)*(2*R_nc/c-Nrg/2/Fr))); 
% H_shift = 1;
s1_ff_shift = s1_ff.*H_shift;                             % 消除各种偏置，否则最后目标点出现的位置不对
clear s1_ff H_shift
% 映射原则：fr_mtx_stolt+整数倍的Fr，直到fr_origin在[-Fr/2,Fr/2]区间内
% 如果频率映射关系错误，斜视角非零时会出现重影、非参考距离处的距离多普勒域散焦、最终结果有多余横线、点目标二维频谱中间列不连续等一系列问题
fr_mtx_stolt = 0:Fr/Nrg:(Nrg-1)*Fr/Nrg;
fr_mtx_stolt = repmat(fr_mtx_stolt,Naz,1);
% 反推fr_mtx_stolt所在的区间，注意stolt映射主要表现为加偏置，区间长度近似不变。（映射后的区间比映射前略长，因此映射回去变短）
% fr_mtx_stolt = floor(((sqrt((f0-Fr/2)^2-c^2*fa_mtx.^2/(4*Vr^2))+sqrt((f0+Fr/2)^2-c^2*fa_mtx.^2/(4*Vr^2)))/2-f0+Fr/2-fr_mtx_stolt)/Fr)*Fr+fr_mtx_stolt;
fr_mtx_stolt_mid = (sqrt((f0-Fr/2)^2-c^2*repmat(fa_mtx,1,Nrg).^2/(4*Vr^2))+sqrt((f0+Fr/2)^2-c^2*repmat(fa_mtx,1,Nrg).^2/(4*Vr^2)))/2-f0;
fr_mtx_stolt = round((fr_mtx_stolt_mid-fr_mtx_stolt)/Fr)*Fr+fr_mtx_stolt;
clear fr_mtx_stolt_mid
fr_origin = sqrt((f0+fr_mtx_stolt).^2+c^2*repmat(fa_mtx,1,Nrg).^2/(4*Vr^2))-f0; % 从fr_mtx_stolt反推fr
clear fr_mtx_stolt
% fr_origin = fr_mtx_stolt;  % 调试用，验证插值
R = 8;                        % 插值核长度
window = kaiser(R,2.5)';      % Kaiser窗
r_TBD = fr_origin*Nrg/Fr+1;   % 转换量级，求得fr在s1_ff_new里对应的坐标（注意+1）
clear fr_origin
point = ceil(r_TBD);          % 取整
s2_ff = zeros(Naz,Nrg);
% 插值关系式：s2_ff(a,fr_mtx_stolt) = s1_ff(a,fr_origin)
for a = 1:Naz
    for r = 1:Nrg
        points = point(a,r) + (-R/2:1:R/2-1);         % 参与插值运算的点的坐标
        kernel = window.*sinc(r_TBD(a,r)-points);     % 参与插值运算的核函数
        kernel = kernel./sum(kernel);                 % 核函数归一化
        points = ceil(mod(points-0.1,Nrg));           % 假设周期循环
        s2_ff(a,r) = s1_ff_shift(a,points)*kernel';     % 计算插值结果
    end
end
clear s1_ff_shift r_TBD point
H_shift = exp(1j*2*pi*(repmat(fr_mtx,Naz,1)*(Nrg/2/Fr)+repmat(fa_mtx,1,Nrg)*(nc)));
s2_ff = s2_ff./H_shift;      % 和s2_tt_2 = circshift(circshift(s2_tt,Nrg/2,2),Naz/2,1);是等效的
clear H_shift
figure;imagesc(abs(ifft2(s2_ff)));
img_wk = ifft2(s2_ff);             % 完整wKA的结果（没有cos(theta_r_c)尺度）
end