function img_cs = CSA(s0,theta_bw,lambda,Kr,Tr,Fr,theta_rc,Nrg,Naz,near_range,Vr,PRF,flag)
%   Chirp Scaling�㷨����
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
center_Rc = near_range + Nrg/2*delta_r;
R_ref = center_Rc*cos(theta_rc);
BW_r = abs(Kr)*Tr;
gama_wr = 1.18;
pr = 0.886*gama_wr/BW_r;
start = near_range*2/c;
% Nr = ceil(Fr*Tr/2)*2;
% ��λ��
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
%% ����������� ���
if(flag == 1)
    figure;subplot(221);
    imagesc(real(s0));
    xlabel('�����򣨲����㣩');ylabel('��λ�򣨲����㣩');title('(a)ԭʼ�ź�ʵ��');
end
Srd = fft(s0,Naz,1);
clear s0;
if(flag == 1)
    subplot(222);
    imagesc(abs(Srd));
    xlabel('�����򣨲����㣩');ylabel('��λƵ�ʣ������㣩');title('(b)ԭʼ�źŵľ����������');
end
f_eta_mtx = (0:Fa/Naz:(Naz-1)*Fa/Naz).';                   
f_eta_mtx = (round((f_etac-f_eta_mtx)/Fa)*Fa+f_eta_mtx); 
D = sqrt(1-c^2*f_eta_mtx.^2/(4*Vr^2*f0^2));                         % �㶯����
Km = Kr./(1-Kr*c*R_ref.*f_eta_mtx.^2./(2*Vr^2*f0^3*D.^3));             % ���������������Km�������ı�
clear f_eta_mtx;
% ssc_ft = exp(1j*pi*Km.*(D_ref./D_mtx-1).*tr_mtx_new.^2);             % ��귽��
ssc = exp(1j*pi*repmat(Km.*(D_ref./D-1),1,Nrg).*(repmat(start +(0:(Nrg-1))/Fr,Naz,1) - repmat(2*R_ref./(c*D),1,Nrg)).^2);             % ��귽��
clear Km D;
S1 = ssc.*Srd;           % ���귽�����
clear Srd ssc;
%% ��άƵ�� ���봦��
S2 = fft(S1,Nrg,2);         % ������FFT
clear S1;
f_eta_mtx = (0:Fa/Naz:(Naz-1)*Fa/Naz).';                   
f_eta_mtx = (round((f_etac-f_eta_mtx)/Fa)*Fa+f_eta_mtx); 
D = sqrt(1-c^2*f_eta_mtx.^2/(4*Vr^2*f0^2));                   % �㶯����
Km = Kr./(1-Kr*c*R_ref.*f_eta_mtx.^2./(2*Vr^2*f0^3*D.^3));             % ���������������Km�������ı�
clear f_eta_mtx;
f_tau_mtx = ifftshift((-Nrg/2:Nrg/2-1)*Fr/Nrg);
% H_range_bulk = exp(1j*pi*(D_mtx.*f_tau_mtx.^2./(Km.*D_ref)+(4*R_ref*f_tau_mtx.*(1./D-1/D_ref)/c)));  % �˲�������ռ�ڴ��1�� 
H_range_bulk = exp(1j*pi*(repmat(D./(Km.*D_ref),1,Nrg).*repmat(f_tau_mtx.^2,Naz,1)+(4*R_ref*repmat(f_tau_mtx,Naz,1).*repmat((1./D-1/D_ref)/c,1,Nrg))));  % �˲���
clear Km D f_tau_mtx;
N_BW_r = round(BW_r/Fr*Nrg);            % Kr*Tr�����ĵ���
window_r = ifftshift(kaiser(N_BW_r,2.5).');
S3 = S2.*H_range_bulk;
clear S2 H_range_bulk
S3 = S3.*repmat([window_r(1:ceil(N_BW_r/2)),zeros(1,Nrg-N_BW_r),window_r(ceil(N_BW_r/2)+1:N_BW_r)],Naz,1);           % �˲�   
clear window_r;
%% ����������� ��λ����
S4 = ifft(S3,Nrg,2);       % ������IFFT
clear S3;
if(flag == 1)
    subplot(223);
    imagesc(abs(S4));
    xlabel('�����򣨲����㣩');ylabel('��λƵ�ʣ������㣩');title('(c)���봦���ľ����������');
end
R0_mtx = (c/2)*(start +(0:(Nrg-1))/Fr)*cos(theta_rc);  % ������߱仯���������
f_eta_mtx = (0:Fa/Naz:(Naz-1)*Fa/Naz).';                   
f_eta_mtx = (round((f_etac-f_eta_mtx)/Fa)*Fa+f_eta_mtx); 
D = sqrt(1-c^2*f_eta_mtx.^2/(4*Vr^2*f0^2));                   % �㶯���� 
Km_mtx = Kr./(1-Kr*c*repmat(R0_mtx,Naz,1).*repmat(f_eta_mtx.^2./(2*Vr^2*f0^3*D.^3),1,Nrg));    % �����仯�ľ������Ƶ��
% H_az_add_offset = exp(1j*4*pi*(f0*R0_mtx.*D_mtx/c-Km_mtx.*(1-D_mtx./D_ref).*(R0_mtx./D_mtx-R0./D_mtx).^2./c^2)-1j*2*pi*fa_dop*nc);
S4 = S4.*exp(1j*4*pi*(f0*repmat(R0_mtx,Naz,1).*repmat(D,1,Nrg)/c));
S5 = S4.*exp(1j*4*pi*(-Km_mtx.*repmat((1-D./D_ref),1,Nrg)...
    .*(repmat(R0_mtx,Naz,1)./repmat(D,1,Nrg)-repmat(R_ref./D,1,Nrg)).^2./c^2));
clear S4 Km_mtx R0_mtx D
S5 = S5.*exp(-1j*2*pi*repmat(f_eta_mtx,1,Nrg)*eta_c);
% H_az_add_offset = exp(1j*4*pi*(-Km_mtx.*repmat((1-D./D_ref),1,Nrg)...
%     .*(repmat(R0_mtx,Naz,1)./repmat(D,1,Nrg)-repmat(R_ref./D,1,Nrg)).^2./c^2)-1j*2*pi*repmat(f_eta_mtx,1,Nrg)*eta_c);   % �˲�������ռ�ڴ��1��
clear f_eta_mtx 
% S5 = S4.*H_az_add_offset;                  % �˲�
img_cs = ifft(S5,Naz,1);                % ��λ��IFFT
if(flag == 1)
    subplot(224);
    imagesc(abs(img_cs));
    xlabel('�����򣨲����㣩');ylabel('��λ�򣨲����㣩');title('(d)CSA������');
end
end