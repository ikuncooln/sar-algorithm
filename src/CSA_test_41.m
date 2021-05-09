%% Chirp Scaling Algorithm
% ֱ�ӵ����У��ܳ��Ľ����Ӧ���ϵ�ͼ
clear;
close all;
clc;
%% ��4.1C�������ز���
% ���������
R_nc = 850e3;   % ������б��
Height = 800e3; % �߶�
Tr = 40e-6;     % ��������ʱ��
Kr = -0.5e12;    % ���������Ƶ��
BW_r = abs(Kr*Tr);   % �źŴ���
Fr = 24e6;      % ���������
width_r = 50e3; % б���������
% ��λ�����
Vr = 7100;      % ��Ч�״��ٶ�
Vs = Vr*1.06;   % ���������ٶ�
Vg = Vr*0.94;   % �����������ص�����ƶ��ٶ�
f0 = 5.3e9;     % �״﹤��Ƶ��
c = 3e8;        % ����
lambda = c/f0;   % �״﹤������
Ls = 4.8e3;     % �ϳɿ׾�����
La = 10;        % ���߳���
BW_dop = 1338;  % �����մ���
Fa = 1700;      % ��λ������
theta_r_c = 8*pi/180;    %��Чб�ӽ�
% R_nc = (1000e3+20e3)/cos(theta_r_c);
% R_nc = 1000e3/cos(theta_r_c);
% ��������
Naz = 1024;    % ��λ���������
Nrg = 1024;    % �������������

A0 = 1;        % ������
% ��������
theta_sq_c = theta_r_c;
beta_bw = 0.886*lambda/La;                      % �״�3dB����
R0 = R_nc*cos(theta_r_c);	                   % ��R_nc���Ӧ�����б�࣬��ΪR0
% R_ref = R0-20e3; % �ο�����
R_ref = R0;
Ta = 0.886*lambda*R_nc/(La*Vg*cos(theta_r_c));  % Ŀ������ʱ��
Ka = 2*Vr^2*cos(theta_r_c)^3/(lambda*R0);       % ��λ���Ƶ��
f_nc = 2*Vr*sin(theta_r_c)/lambda;              % ����������Ƶ��
nc = (-R0*tan(theta_r_c))/Vr;                  % �����ĵĲ������Ĵ�Խʱ��
Nr = round(Tr*Fr);                             % ���Ե�Ƶ�źŲ�������
Na = round(Ta*Fa);                             % ��λ���������?
a_os_r = Fr/BW_r;                              % ���������������
a_os_a = Fa/BW_dop;                            % ��λ�����������
Mamb = floor(f_nc/Fa);                         % ������ģ��

% Nrg = ceil(Tr*Fr);
% Naz = ceil(Ta*Fa);
%% Ŀ�����ز���
delta_r1 = -1500;         % Ŀ��1�;����ĵľ���������
delta_a1 = delta_r1*tan(theta_r_c);            % Ŀ��1�;����ĵķ�λ������
% delta_r2 = 0;            % Ŀ��2�;����ĵľ���������     
delta_r2 = -20e3;            % Ŀ��2�;����ĵľ���������    
% delta_a2 = 0;            % Ŀ��2�;����ĵķ�λ������
delta_a2 = delta_r2*tan(theta_r_c);    
delta_r3 = 1500;          % Ŀ��3�;����ĵľ���������     
delta_a3 = delta_r3*tan(theta_r_c);            % Ŀ��3�;����ĵķ�λ������    

% Ŀ��1
x1 = R0 + delta_r1;      % Ŀ��1�ľ��������
y1 = delta_a1;           % Ŀ��1�ķ�λ�����
% Ŀ��2
x2 = R0 + delta_r2;      % Ŀ��2�ľ��������
y2 = delta_a2;           % Ŀ��2�ķ�λ�����
% Ŀ��3
x3 = R0 + delta_r3;      % Ŀ��3�ľ��������
y3 = delta_a3;           % Ŀ��3�ķ�λ�����
% ��������
x_range = [x1,x2,x3];
y_azimuth = [y1,y2,y3];
% ���㲨�����Ĵ�Խʱ��
nc_1 = (y1-x1*tan(theta_r_c))/Vr;           % Ŀ��1�Ĳ������Ĵ�Խʱ�̡�
nc_2 = (y2-x2*tan(theta_r_c))/Vr;           % Ŀ��2�Ĳ������Ĵ�Խʱ�̡�
nc_3 = (y3-x3*tan(theta_r_c))/Vr;           % Ŀ��3�Ĳ������Ĵ�Խʱ�̡�
nc_target = [nc_1,nc_2,nc_3];               % ��������飬���ڴ���
%% ���������ͷ�λ���ʱ�����Ƶ����
% tr = 2*R_nc/c+(-Nrg/2:(Nrg/2-1))/Fr;           % ����ʱ���ᣬ����Ϊ2*R_nc/c
tr = 2*x2/cos(theta_r_c)/c+(-Nrg/2:(Nrg/2-1))/Fr;
ta = nc+(-Naz/2:Naz/2-1)/Fa;                   % ��λʱ���ᣬ��y=0�������ģ����������ʱ��Ϊ��λ��ʱ����� 
% ta = nc_2+(-Naz/2:Naz/2-1)/Fa;                   % ��λʱ���ᣬ��y=0�������ģ����������ʱ��Ϊ��λ��ʱ����� 
% ���ɾ��루��λ��ʱ�䣨Ƶ�ʣ�����
tr_mtx = ones(Naz,1)*tr;                       % ����ʱ������󣬴�С��Naz*Nrg
ta_mtx = ta.'*ones(1,Nrg);                     % ��λʱ������󣬴�С��Naz*Nrg
% ����Ƶ����
fr_mtx = 0:Fr/Nrg:(Nrg-1)*Fr/Nrg;   % �������ѱ����������              
fr_mtx = floor((Fr/2-fr_mtx)/Fr)*Fr+fr_mtx; 
fr_mtx = ones(Naz,1)*fr_mtx;
% ��λƵ����
fa_mtx = 0:Fa/Naz:(Naz-1)*Fa/Naz;     % Ҫ����ЩƵ��ֵ+��������Fa��ӳ�䵽(f_nc-Fa/2,f_nc+Fa/2]������              
fa_mtx = floor((f_nc+Fa/2-fa_mtx)/Fa)*Fa+fa_mtx; 
fa_mtx = fa_mtx.'*ones(1,Nrg);
%% �״�ԭʼ����
s0_tt = zeros(Naz,Nrg);                                                                   % ʱ���ź�
% for k = 1:3
for k = 2
    R_n = sqrt(x_range(k)^2+(Vr*ta_mtx-y_azimuth(k)).^2);                                 % ˲ʱб��
    wr = (abs(tr_mtx-2*R_n/c)<=Tr/2);                                                     % �������
    wa = sinc(0.886*(atan(Vg*(ta_mtx-nc_target(k))/x_range(k)))/beta_bw).^2;              % ��λ����
    s0_tt = s0_tt + A0*wr.*wa.*exp(-1j*4*pi*f0*R_n/c).*exp(1j*pi*Kr*(tr_mtx-2*R_n/c).^2); % �����ź�
end
%% �����������  
Srd_ft = fft(s0_tt,Naz,1);  % ��λ��FFT
figure;
subplot(3,2,[1,3,5]);
imagesc(real(Srd_ft));
set(gca,'YDir','normal');
title('(a)ԭʼ�źŵľ��������Ƶ��');
xlabel('���루�����㣩');
ylabel('��λƵ�ʣ������㣩');
sample = [300,600,800];  % ѡȡ3����λƵ�ʵ�۲죬����λ��������仯��Ҫ���ݾ����������۲�
subplot(322);
plot(real(Srd_ft(sample(1),:)));
xlim([1,Nrg]);
ylabel('��һ������');
title(['(b)��',num2str(sample(1)),'����λƵ�ʲ���']);
subplot(324);
plot(real(Srd_ft(sample(2),:)));
xlim([1,Nrg]);
ylabel('��һ������');
title(['(c)��',num2str(sample(2)),'����λƵ�ʲ���']);
subplot(326);
plot(real(Srd_ft(sample(3),:)));
xlim([1,Nrg]);
xlabel('���루�����㣩');
ylabel('��һ������');
title(['(d)��',num2str(sample(3)),'����λƵ�ʲ���']);
%% һ�ºͲ�������㶯�ķ���
D_mtx = sqrt(1-c^2*fa_mtx.^2/(4*Vr^2*f0^2));     % �㶯����
tr_RCMC = tr*cos(theta_r_c);                     % ���µľ����߳����µ�ʱ���ᣬע������ֵ�����˱仯
R0_mtx = ones(Naz,1)*(c/2)*tr_RCMC;              % ������߱仯���������
% D_ref = D_mtx(ceil(mod(round(mod(f_nc,Fa)*Naz/Fa)-0.1,Naz)),1);    % �ο���λƵ�ʣ�����������Ƶ�ʣ������㶯����
D_ref = cos(theta_r_c);
RCM_bulk = 2*(R_ref./D_mtx-R_ref./D_ref)*Fr/c;                                     % һ�¾����㶯
RCM_diff = 2*(R0_mtx./D_mtx - R0_mtx./D_ref - R_ref./D_mtx + R_ref./D_ref)*Fr/c;   % ��������㶯
% fa_mtx = floor((f_nc+Fa/2-fa_mtx)/Fa)*Fa+fa_mtx; 
dif = diff(fa_mtx(:,1));
p_break = find(abs(dif)==max(abs(dif)));
sample_r = Nrg/2+1;
figure;
plot(RCM_bulk(1:p_break,sample_r),1:p_break,'r',RCM_diff(1:p_break,sample_r),1:p_break,'b--',...
    RCM_bulk(p_break+1:Naz,sample_r),p_break+1:Naz,'r',RCM_diff(p_break+1:Naz,sample_r),p_break+1:Naz,'b--');
legend('һ�¾����㶯','��������㶯','location','Bestoutside');
ylim([1,Naz]);
set(gca,'YDir','reverse');
title('һ�ºͲ�������㶯�ķ���');
xlabel('�����㶯�������㣩');
ylabel('��λƵ�ʣ������㣩');
%% ֻ���о���ѹ����һ�¾����㶯У��
s2_ff = fft(Srd_ft,Nrg,2);                       % ������FFT���任����άƵ��
D_mtx = sqrt(1-c^2*fa_mtx.^2/(4*Vr^2*f0^2));     % �㶯����
Km = Kr./(1-Kr*c*R_ref.*fa_mtx.^2./(2*Vr^2*f0^3*D_mtx.^3));       % �����仯�ľ������Ƶ��
% D_ref = D_mtx(ceil(mod(round(mod(f_nc,Fa)*Naz/Fa)-0.1,Naz)),1);    % �ο���λƵ�ʣ�����������Ƶ�ʣ������㶯����
D_ref = cos(theta_r_c);
H_range = exp(1j*pi*D_mtx.*fr_mtx.^2./(Km.*D_ref));            % �������������
H_RCM_bulk = exp(1j*4*pi*R_ref*fr_mtx.*(1./D_mtx-1/D_ref)/c);         % һ�¾����㶯У��
N_BW_r = round(Nrg/a_os_r);    
window = ifftshift(kaiser(N_BW_r,2.5)');    % Kaiser��
window = (ones(Naz,1)*[window(1:ceil(N_BW_r/2)),zeros(1,Nrg-N_BW_r),window(ceil(N_BW_r/2)+1:N_BW_r)]);
s3_ff = s2_ff.*H_range.*H_RCM_bulk.*window;    % �˲�     
s4_ft = ifft(s3_ff,Nrg,2);                     % ������IFFT
N_cut = 32;
figure;
subplot(3,2,[1,3,5]);
imagesc(abs(s4_ft));
xlim([Nrg/2+1-N_cut/2,Nrg/2+N_cut/2]);
set(gca,'YDir','normal');
title('(a)����ѹ����ľ��������Ƶ��');
xlabel('���루�����㣩');
ylabel('��λƵ�ʣ������㣩');
% ����������
Num = N_cut*8;      % ����������
% s4_ff = fft(s4_ft,Nrg,2);
% s4_ft_up = ifft([s4_ff(:,1:Nrg/2),zeros(Naz,Num-Nrg),s4_ff(:,Nrg/2+1:Nrg)],Num,2);
% figure;imagesc(abs(s4_ft_up));
% set(gca,'YDir','normal');
% title('������-����ѹ����һ��RCMC��');
% xlabel('����ʱ�䣨�����㣩');
% ylabel('��λƵ�ʣ������㣩');
% ѡȡ��20��60��240����λƵ�ʲ�����ľ���������������
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
text((I-1)*N_cut/Num+4,M,['��ֵλ�ڵ�',num2str((I-1)*N_cut/Num),'��������']);
ylabel('����');
title(['(b)��',num2str(sample(1)),'����λƵ�ʲ���']);
subplot(324);
plot(0:N_cut/Num:N_cut-N_cut/Num,abs(line2_up));
xlim([0,N_cut]);
[M,I] = max(abs(line2_up));
text((I-1)*N_cut/Num+4,M,['��ֵλ�ڵ�',num2str((I-1)*N_cut/Num),'��������']);
ylabel('����');
title(['(c)��',num2str(sample(2)),'����λƵ�ʲ���']);
subplot(326);
plot(0:N_cut/Num:N_cut-N_cut/Num,abs(line3_up));
xlim([0,N_cut]);
[M,I] = max(abs(line3_up));
text((I-1)*N_cut/Num+4,M,['��ֵλ�ڵ�',num2str((I-1)*N_cut/Num),'��������']);
xlabel('���루�����㣩');
ylabel('����');
title(['(d)��',num2str(sample(3)),'����λƵ�ʲ���']);
%% ���б�ꡢ����ѹ����һ�¾����㶯У��
% ��귽��
D_mtx = sqrt(1-c^2*fa_mtx.^2/(4*Vr^2*f0^2));                       % �㶯����
% D_ref = D_mtx(ceil(mod(round(mod(f_nc,Fa)*Naz/Fa)-0.1,Naz)),1);    % �ο���λƵ�ʣ�����������Ƶ�ʣ������㶯����
D_ref = cos(theta_r_c);
tr_mtx_new = tr_mtx - 2*R_ref./(c*D_mtx);                             % �µľ���ʱ��
Km = Kr./(1-Kr*c*R_ref.*fa_mtx.^2./(2*Vr^2*f0^3*D_mtx.^3));           % ���������������Km�������ı�
ssc_ft = exp(1j*pi*Km.*(D_ref./D_mtx-1).*tr_mtx_new.^2);              % ��귽��
S1_ft = ssc_ft.*Srd_ft;                                            % ���귽�����
s2_ff = fft(S1_ft,Nrg,2);                                          % ������FFT
% ����ѹ����һ�¾����㶯У��
Km = Kr./(1-Kr*c*R_ref.*fa_mtx.^2./(2*Vr^2*f0^3*D_mtx.^3));       
H_range = exp(1j*pi*D_mtx.*fr_mtx.^2./(Km.*D_ref));            % �������������
H_RCM_bulk = exp(1j*4*pi*R_ref*fr_mtx.*(1./D_mtx-1/D_ref)/c);         % һ�¾����㶯У��
N_BW_r = round(BW_r/Fr*Nrg);            % Kr*Tr�����ĵ���
window = ifftshift(kaiser(N_BW_r,2.5)');    % Kaiser��
window = (ones(Naz,1)*[window(1:ceil(N_BW_r/2)),zeros(1,Nrg-N_BW_r),window(ceil(N_BW_r/2)+1:N_BW_r)]);
% window = ifftshift(ones(Naz,1)*kaiser(Nrg,2.5)');         % ����ƽ����
s3_ff = s2_ff.*H_range.*H_RCM_bulk.*window;    % �˲�     
s4_ft = ifft(s3_ff,Nrg,2);                     % ������IFFT

N_cut = 32;
figure;
subplot(3,2,[1,3,5]);
imagesc(abs(s4_ft));
xlim([Nrg/2+1-N_cut/2,Nrg/2+N_cut/2]);
set(gca,'YDir','normal');
title('(a)����ѹ����ľ��������Ƶ��');
xlabel('���루�����㣩');
ylabel('��λƵ�ʣ������㣩');
% ����������
Num = N_cut*8;      % ����������
% s4_ff = fft(s4_ft,Nrg,2);
% s4_ft_up = ifft([s4_ff(:,1:Nrg/2),zeros(Naz,Num-Nrg),s4_ff(:,Nrg/2+1:Nrg)],Num,2);
% figure;imagesc(abs(s4_ft_up));
% set(gca,'YDir','normal');
% title('������-����ѹ����һ��RCMC��');
% xlabel('����ʱ�䣨�����㣩');
% ylabel('��λƵ�ʣ������㣩');
% ѡȡ��20��60��240����λƵ�ʲ�����ľ���������������
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
text((I-1)*N_cut/Num+4,M,['��ֵλ�ڵ�',num2str((I-1)*N_cut/Num),'��������']);
ylabel('����');
title(['(b)��',num2str(sample(1)),'����λƵ�ʲ���']);
subplot(324);
plot(0:N_cut/Num:N_cut-N_cut/Num,abs(line2_up));
xlim([0,N_cut]);
[M,I] = max(abs(line2_up));
text((I-1)*N_cut/Num+4,M,['��ֵλ�ڵ�',num2str((I-1)*N_cut/Num),'��������']);
ylabel('����');
title(['(c)��',num2str(sample(2)),'����λƵ�ʲ���']);
subplot(326);
plot(0:N_cut/Num:N_cut-N_cut/Num,abs(line3_up));
xlim([0,N_cut]);
[M,I] = max(abs(line3_up));
text((I-1)*N_cut/Num+4,M,['��ֵλ�ڵ�',num2str((I-1)*N_cut/Num),'��������']);
xlabel('���루�����㣩');
ylabel('����');
title(['(d)��',num2str(sample(3)),'����λƵ�ʲ���']);
%% ��λ����
tr_RCMC = tr*cos(theta_r_c);                     % ���µľ����߳����µ�ʱ���ᣬע������ֵ�����˱仯
% tr_RCMC = 2*R0/c + (-Nrg/2:(Nrg/2-1))/Fr;
% D_ref = D_mtx(ceil(mod(round(mod(f_nc,Fa)*Naz/Fa)-0.1,Naz)),1);
R0_mtx = ones(Naz,1)*(c/2)*tr_RCMC;              % ������߱仯���������
Km_mtx = Kr./(1-Kr*c*R0_mtx.*fa_mtx.^2./(2*Vr^2*f0^3*D_mtx.^3));   % �����仯�ľ������Ƶ��
H_az = exp(1j*4*pi*f0*R0_mtx.*D_mtx/c);           % ��λ��ƥ���˲�
H_add = exp(-1j*4*pi*Km_mtx.*(1-D_mtx./D_ref).*(R0_mtx./D_mtx-R_ref./D_mtx).^2./c^2);   % ������λУ��
H_offset = exp(-1j*2*pi*(fa_mtx*nc));            % �ѷ�λʱ�����ĸ�Ϊ0
% H_offset = exp(-1j*2*pi*(fa_mtx*nc_2));            % �ѷ�λʱ�����ĸ�Ϊ0
% window = kaiser(Naz,2.5)*ones(1,Nrg);         % ��λƽ����
s5_ft = s4_ft.*H_az.*H_add.*H_offset;          % �˲�
s5_tt = ifft(s5_ft,Naz,1);           % ��λ��IFFT
figure;
imagesc(abs(s5_tt));
title('CSA������');
xlabel('����ʱ�䣨�����㣩');
%% ��Ŀ��Ŵ�
area = -8:1:7;       % ������Ƭ�Ĵ�С
% Ŀ��2
pos_2_a = round(Naz/2 + 1 + y2*Fa/Vr);          % Ŀ��2��λ���λ��
% pos_2_a = round(Naz/2 + 1);          % Ŀ��2��λ���λ��
% pos_2_r = round(Nrg/2 + 1 + 2*(x2-R0)*Fr/c/cos(theta_r_c));    % Ŀ��2�������λ��
pos_2_r = round(Nrg/2 + 1);    % Ŀ��2�������λ��
amplify_2 = s5_tt(ceil(mod(pos_2_a+area-0.1,Naz)),ceil(mod(pos_2_r+area-0.1,Nrg)));
figure;imagesc(abs(amplify_2));
title('��Ŀ��2�Ŵ󣨷��ȣ�');
xlabel('����ʱ�䣨�����㣩');
ylabel('��λʱ�䣨�����㣩');
target = amplify_2;        % ѡ��Ŀ����з���
[image_upsample,signal_r,quality_r,signal_a,quality_a] = f_point_analyse(target,c/2/Fr,Vg/Fa,32,1);
IRW_r_theory = c/2/BW_r*0.886*1.18;
IRW_a_theory = La/2*Vg/Vs*1.185;


% [image_upsample,az_cut,rg_cut,IRW_a,IRW_r,PSLR_a,PSLR_r] = f_point_analyze_trial(target,16,1);
% Len = length(area);
% Num = size(image_upsample,1);
% image_upsample_dB = 20*log10(abs(image_upsample)/max(max(abs(image_upsample))));
% angle_upsample = angle(image_upsample);
% figure;
% subplot(321);imagesc(0:Len/Num:Len-Len/Num,0:Len/Num:Len-Len/Num,abs(image_upsample));title('(a)�Ŵ��ĵ�Ŀ��');xlabel('�����򣨲����㣩');ylabel('��λ�򣨲����㣩');
% subplot(322);contour(0:Len/Num:Len-Len/Num,0:Len/Num:Len-Len/Num,abs(image_upsample),30);title('(b)�Ŵ��ĵ�Ŀ���ֵ��ͼ');set(gca,'YDir','reverse');xlabel('�����򣨲����㣩');ylabel('��λ�򣨲����㣩');
% subplot(323);plot(0:Len/Num:Len-Len/Num,image_upsample_dB(az_cut,:));xlim([0,Len]);ylim([-35,0]);title('(c)��������ͼ');xlabel('����ʱ�䣨�����㣩');ylabel('����(dB)');
% subplot(324);plot(0:Len/Num:Len-Len/Num,image_upsample_dB(:,rg_cut));xlim([0,Len]);ylim([-35,0]);title('(d)��λ����ͼ');xlabel('��λʱ�䣨�����㣩');ylabel('����(dB)');
% subplot(325);plot(0:Len/Num:Len-Len/Num,angle_upsample(az_cut,:)*180/pi);xlim([0,Len]);title('(e)������λ');xlabel('����ʱ�䣨�����㣩');ylabel('��λ���㣩');
% subplot(326);plot(0:Len/Num:Len-Len/Num,angle_upsample(:,rg_cut)*180/pi);xlim([0,Len]);title('(f)��λ��λ');xlabel('��λʱ�䣨�����㣩');ylabel('��λ���㣩');