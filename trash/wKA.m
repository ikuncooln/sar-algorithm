function img_wk = wKA(s0,theta_bw,lambda,Kr,Tr,Fr,theta_rc,Nrg,Naz,near_range,Vr,PRF,flag)
%   wK�㷨����
%   s0�������ź����ݣ���ά�����飩
%   theta_bw�����߷�λ�������(rad)
%   lambda�ǲ���(m)
%   Kr�Ǿ������Ƶ��(Hz/s)
%   Tr�Ƿ���������(s)
%   Fr�Ǿ����������(Hz)
%   theta_rc��б�ӽ�(rad)
%   Nrg�Ǿ������������
%   Naz�Ƿ�λ���������
%   near_range�ǵ�һ�����������(m)
%   Vr���ػ��ٶ�(m/s)
%   PRF�������ظ�Ƶ��(Hz)
%   flagΪ1/0��ʾ��/�������м䲽���ͼ
%% ��������
c=299792458;
% ������
delta_r = c/2/Fr;
R_nc = near_range + Nrg/2*delta_r;
R_ref = R_nc*cos(theta_rc);
BW_r = abs(Kr)*Tr;
gama_wr = 1.18;
pr = 0.886*gama_wr/BW_r;
start = near_range*2/c;
a_os_r = Fr/BW_r;
% Nr = ceil(Fr*Tr/2)*2;
% ��λ��
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
% ����Ƶ����
fr_mtx = ifftshift((-Nrg/2:Nrg/2-1)*Fr/Nrg);
% ��λƵ����
fa_mtx = (0:Fa/Naz:(Naz-1)*Fa/Naz).';     % Ҫ����ЩƵ��ֵ+��������Fa��ӳ�䵽[f_nc-Fa/2,f_nc+Fa/2)������
fa_mtx = round((f_nc-fa_mtx)/Fa)*Fa+fa_mtx;
%% ��άFFT
s0_ff = fft2(s0,Naz,Nrg);  % ��άFFT
clear s0
%% �ο�������ˣ�һ��ѹ����
% H_ref = exp(1j.*(4*pi*R_ref/c*sqrt((f0+repmat(fr_mtx,Naz,1)).^2-c^2.*repmat(fa_mtx,1,Nrg).^2/(4*Vr^2))+pi*repmat(fr_mtx,Naz,1).^2/Kr));
s0_ff = s0_ff.*exp(1j.*(4*pi*R_ref/c*sqrt((f0+repmat(fr_mtx,Naz,1)).^2-c^2.*repmat(fa_mtx,1,Nrg).^2/(4*Vr^2))+pi*repmat(fr_mtx,Naz,1).^2/Kr));
clear H_ref;
N_BW_r = round(Nrg/a_os_r);            % Kr*Tr�����ĵ���
window_r = ifftshift(kaiser(N_BW_r,2.5)');    % Kaiser��
window_r = repmat([window_r(1:ceil(N_BW_r/2)),zeros(1,Nrg-N_BW_r),window_r(ceil(N_BW_r/2)+1:N_BW_r)],Naz,1);
s1_ff = s0_ff.*window_r;       % �ο��������
clear s0_ff window_r
%% Stoltӳ��
% H_shift = exp(-1j*2*pi*(repmat(fr_mtx,Naz,1)*(-Nrg/2/Fr)+fa_mtx*(-Naz/2/Fa)));  % ������ͷ�λʱ��������������Ų����ʼ
H_shift = exp(-1j*2*pi*(repmat(fr_mtx,Naz,1)*(2*R_nc/c-Nrg/2/Fr))); 
% H_shift = 1;
s1_ff_shift = s1_ff.*H_shift;                             % ��������ƫ�ã��������Ŀ�����ֵ�λ�ò���
clear s1_ff H_shift
% ӳ��ԭ��fr_mtx_stolt+��������Fr��ֱ��fr_origin��[-Fr/2,Fr/2]������
% ���Ƶ��ӳ���ϵ����б�ӽǷ���ʱ�������Ӱ���ǲο����봦�ľ����������ɢ�������ս���ж�����ߡ���Ŀ���άƵ���м��в�������һϵ������
fr_mtx_stolt = 0:Fr/Nrg:(Nrg-1)*Fr/Nrg;
fr_mtx_stolt = repmat(fr_mtx_stolt,Naz,1);
% ����fr_mtx_stolt���ڵ����䣬ע��stoltӳ����Ҫ����Ϊ��ƫ�ã����䳤�Ƚ��Ʋ��䡣��ӳ���������ӳ��ǰ�Գ������ӳ���ȥ��̣�
% fr_mtx_stolt = floor(((sqrt((f0-Fr/2)^2-c^2*fa_mtx.^2/(4*Vr^2))+sqrt((f0+Fr/2)^2-c^2*fa_mtx.^2/(4*Vr^2)))/2-f0+Fr/2-fr_mtx_stolt)/Fr)*Fr+fr_mtx_stolt;
fr_mtx_stolt_mid = (sqrt((f0-Fr/2)^2-c^2*repmat(fa_mtx,1,Nrg).^2/(4*Vr^2))+sqrt((f0+Fr/2)^2-c^2*repmat(fa_mtx,1,Nrg).^2/(4*Vr^2)))/2-f0;
fr_mtx_stolt = round((fr_mtx_stolt_mid-fr_mtx_stolt)/Fr)*Fr+fr_mtx_stolt;
clear fr_mtx_stolt_mid
fr_origin = sqrt((f0+fr_mtx_stolt).^2+c^2*repmat(fa_mtx,1,Nrg).^2/(4*Vr^2))-f0; % ��fr_mtx_stolt����fr
clear fr_mtx_stolt
% fr_origin = fr_mtx_stolt;  % �����ã���֤��ֵ
R = 8;                        % ��ֵ�˳���
window = kaiser(R,2.5)';      % Kaiser��
r_TBD = fr_origin*Nrg/Fr+1;   % ת�����������fr��s1_ff_new���Ӧ�����꣨ע��+1��
clear fr_origin
point = ceil(r_TBD);          % ȡ��
s2_ff = zeros(Naz,Nrg);
% ��ֵ��ϵʽ��s2_ff(a,fr_mtx_stolt) = s1_ff(a,fr_origin)
for a = 1:Naz
    for r = 1:Nrg
        points = point(a,r) + (-R/2:1:R/2-1);         % �����ֵ����ĵ������
        kernel = window.*sinc(r_TBD(a,r)-points);     % �����ֵ����ĺ˺���
        kernel = kernel./sum(kernel);                 % �˺�����һ��
        points = ceil(mod(points-0.1,Nrg));           % ��������ѭ��
        s2_ff(a,r) = s1_ff_shift(a,points)*kernel';     % �����ֵ���
    end
end
clear s1_ff_shift r_TBD point
H_shift = exp(1j*2*pi*(repmat(fr_mtx,Naz,1)*(Nrg/2/Fr)+repmat(fa_mtx,1,Nrg)*(nc)));
s2_ff = s2_ff./H_shift;      % ��s2_tt_2 = circshift(circshift(s2_tt,Nrg/2,2),Naz/2,1);�ǵ�Ч��
clear H_shift
figure;imagesc(abs(ifft2(s2_ff)));
img_wk = ifft2(s2_ff);             % ����wKA�Ľ����û��cos(theta_r_c)�߶ȣ�
end