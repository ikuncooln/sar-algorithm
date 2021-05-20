close all;clear all;
%% 1. ������� (�ο� p142, table 6.1)
center_Rc = 20e3;  % ������б��
Vr = 150;   % ��Ч�״��ٶ�
Tr = 2.5e-6;    % ��������ʱ��
Kr = 20e12; % �����Ƶ��
f0 = 5.3e9; % �״﹤��Ƶ��
BW_dop = 80;    % �����մ���
Fr = 60e6;  % ���������
Fa = 100;   % ��λ������
Naz = 256;  % ��λ�����������������������
Nrg = 256;  % ��������������������߲���������
theta_rc_deg = 21.9; % ��б�ӽ�21.9��
c = 299792458;    % ����

% derived params
lambda = c / f0;
theta_rc = theta_rc_deg * pi / 180;
Vs = Vr;
Vg = Vr;
Np = Tr * Fr;   % �������г��ȣ�����������
alpha_os_r = Fr / (Kr*Tr);
alpha_os_a = Fa / BW_dop;

%% 2. ����ԭʼ�״�����
NUM_TARGETS = 3;    % �����Ŀ����Ϊ3
rs = [0, 0, 30];    % ��Ŀ����������
as = [-20, 0, -7.94]; % Ŀ����Է�λ�����
parameters = struct(...
    'center_Rc', center_Rc,...          % ������б��
    'theta_rc_deg', theta_rc_deg,...    % б�ӽ�
    'Nrg', Nrg,...                      % �������������
    'Naz', Naz,...                      % ��λ���������
    'Vr', Vr,...                        % �ػ��ٶ�
    'f0', f0,...                        % �ز�Ƶ��
    'Tr', Tr,...                        % ����������
    'Kr', Kr,...                        % ���������Ƶ��
    'BW_dop', BW_dop,...                % �����մ���
    'alpha_os_r', alpha_os_r,...        % �������������
    'alpha_os_a', alpha_os_a,...        % ��λ���������
    'NUM_TARGETS', NUM_TARGETS,...      % ��Ŀ������
    'rs', rs,...                        % ��Ŀ����������꣨m��
    'as', as...                         % ��Ŀ�귽λ�����꣨m��
);

[ s0, f_etac, delta_r, delta_a, center_R0, center_Rc ] = generate_point_data(parameters);

figure; % ���ƴ�б�ӽ�����µ������״�ԭʼ�����ź�
subplot(221);imagesc(real(s0));ylabel('��λ��ʱ�䣨�����㣩');title('(a)ʵ��');
subplot(222);imagesc(imag(s0));title('(b)�鲿');
subplot(223);imagesc(abs(s0));xlabel('������ʱ�䣨�����㣩');ylabel('��λ��ʱ�䣨�����㣩');title('(c)����');
subplot(224);imagesc(angle(s0));xlabel('������ʱ�䣨�����㣩');title('(d)��λ');
suptitle([num2str(theta_rc_deg), '��б�ӽ�����µ�', num2str(NUM_TARGETS),'���״�ԭʼ�����źţ�ʱ��']);


%% ����ѹ��(���÷�ʽ3ƥ���˲���
f_tau = ifftshift((-Nrg/2:Nrg/2-1)*Fr/Nrg); % ������Ƶ����
Hrc = exp(1j*pi*f_tau.^2/Kr);  % Matched filter in Frequency domain
a_os_r = Fr/abs(Kr*Tr);
N_BW_r = round(Nrg/a_os_r);            % Kr*Tr�����ĵ���
window_r = ifftshift(kaiser(N_BW_r,2.5)');    % Kaiser��
window_r = repmat([window_r(1:ceil(N_BW_r/2)),zeros(1,Nrg-N_BW_r),window_r(ceil(N_BW_r/2)+1:N_BW_r)],Naz,1);
Hrc = repmat(Hrc,Naz,1);
Hrc = Hrc.*window_r;
s0_tmp = fft(s0.').';          %����Ƶ��λʱ�� fftĬ�ϰ���
%ע�����ﲻ��fftshift
Src = s0_tmp.*Hrc;             %ƥ���˲�
s_rc = ifft(Src.').';

%% ��λ����Ҷ�任
Srd = fft(s_rc);
figure; % ���ƾ���ѹ����Ľ��
subplot(221);imagesc(real(s_rc));ylabel('��λ��ʱ�䣨�����㣩');title('(a)���ξ���ѹ��ǰʵ����ʱ��');
subplot(222);imagesc(abs(s_rc));title('(b)���ξ���ѹ��ǰ���ȣ�ʱ��');
subplot(223);imagesc(real(Srd));xlabel('������ʱ�䣨�����㣩');ylabel('��λ��Ƶ�ʣ������㣩');title('(c)���ξ���ѹ��ǰʵ���������������');set(gca, 'YDir', 'normal');
subplot(224);imagesc(abs(Srd));xlabel('������ʱ�䣨�����㣩');title('(d)���ξ���ѹ��ǰ���ȣ������������');set(gca, 'YDir', 'normal');
suptitle('21.9��б�ӽǾ���ѹ�����źţ�ʱ��������������');
%% ���ξ���ѹ��(��άƵ����У�
S2df  = fft(Srd.').';          %ǰ��Ϊ�˹۲�Ƶ�����˸���Ҷ��任��ʵ�ʿ���ʡ�ԡ�
f_eta = (ifftshift((-Naz/2:Naz/2-1)*Fa/Naz)).';
f_eta = f_eta + round((f_etac - f_eta)/Fa)*Fa;
tau0 = 2*center_R0/cos(theta_rc)/c;
tau = (-Nrg/2:Nrg/2-1)/Fr+tau0;  % ����ʱ����
R0 = tau*c/2*cos(theta_rc);
[R0_grid,f_eta_grid] = meshgrid(R0,f_eta);

%����Range-Doppler���еľ����㶯����
D = sqrt(1-lambda^2.*f_eta.^2/4/Vr^2);
%�������ѹ����Ƶ��
Ksrc = 2*Vr^2*f0^3.*D.^3/c/center_R0./f_eta.^2;
f_tau_mtx = repmat(f_tau,Naz,1);
%�������ѹ���˲���
Hsrc = exp(-1j*pi*f_tau_mtx.^2./repmat(Ksrc,1,Nrg));
Ssrc = S2df.*Hsrc;              %��άƵ����ʵ�ֶ���ѹ��
s_src = ifft(Ssrc.').';

figure; % ���ƾ���������������ѹ����Ͷ��ξ���ѹ����Ľ��
subplot(221);imagesc(real(Srd));ylabel('��λ��Ƶ�ʣ������㣩');title('(a)���ξ���ѹ��ǰʵ��');set(gca, 'YDir', 'normal');
subplot(222);imagesc(abs(Srd));title('(b)���ξ���ѹ��ǰ����');set(gca, 'YDir', 'normal');
subplot(223);imagesc(real(s_src));xlabel('������ʱ�䣨�����㣩');ylabel('��λ��Ƶ�ʣ������㣩');title('(c)���ξ���ѹ����ʵ��');set(gca, 'YDir', 'normal');
subplot(224);imagesc(abs(s_src));xlabel('������ʱ�䣨�����㣩');title('(d)���ξ���ѹ�������');set(gca, 'YDir', 'normal');
suptitle('21.9��б�ӽǶ��ξ���ѹ��ǰ���źŶԱȣ������������');

%% �����㶯У����RCMC)
D_grid = repmat(D,1,Nrg);
RCM = R0_grid./D_grid-R0_grid;
RCM = RCM - (1/cos(theta_rc)-1)*R0_grid;
RCM = RCM / delta_r; %�������㶯��ת��Ϊ���뵥Ԫƫ������
% �����ֵ��ϵ����
x_tmp = repmat(-4:3, 16, 1);
offset_tmp = (1:16)/16;
x_tmp = x_tmp + repmat(offset_tmp.', 1, 8);
hx = sinc(x_tmp);
x_tmp16 = x_tmp .* 16;
x_tmp16 = round(x_tmp16 + 16 * 8 / 2);
kwin = repmat(kaiser(16*8, 2.5).', 16, 1);
hx = kwin(x_tmp16) .* hx;
hx = hx ./ repmat(sum(hx, 2), 1, 8);

% ��ֵУ��
Srcmc = zeros(Naz, Nrg);  % ��ž����㶯У����Ļز��ź�
for i = 1:Naz
    for j = 1:Nrg
        offset_int = ceil(RCM(i,j));
        offset_frac = round((offset_int - RCM(i,j)) * 16);
        if offset_frac == 0
            Srcmc(i,j) = s_src(i,ceil(mod(j+offset_int-0.1,Nrg)));   % �����ź�����S1�������Լٶ�
        else
            Srcmc(i,j) = s_src(i, ceil(mod((j+offset_int-4:j+offset_int+3)-0.1,Nrg))) * hx(offset_frac,:).';
        end
        
    end
end
figure; % ���ƾ������������ľ���ѹ�㶯У��ǰ�Ľ��
subplot(121);imagesc(real(s_src));xlabel('������ʱ�䣨�����㣩');ylabel('��λ��Ƶ�ʣ������㣩');title('(a)ʵ��');set(gca, 'YDir', 'normal');
subplot(122);imagesc(abs(s_src));xlabel('������ʱ�䣨�����㣩');title('(b)����');set(gca, 'YDir', 'normal');
suptitle('21.9��б�ӽǾ����㶯У��ǰ�źţ������������');

figure; % ���ƾ������������ľ���ѹ�㶯У����Ľ��
subplot(121);imagesc(real(Srcmc));xlabel('������ʱ�䣨�����㣩');ylabel('��λ��Ƶ�ʣ������㣩');title('(a)ʵ��');set(gca, 'YDir', 'normal');
subplot(122);imagesc(abs(Srcmc));xlabel('������ʱ�䣨�����㣩');title('(b)����');set(gca, 'YDir', 'normal');
suptitle('21.9��б�ӽǾ����㶯У�����źţ������������');

%% ��λ��ѹ��
% Srcmc=s_src;%���������㶯У��
Haz = exp(1j*4*pi.*R0_grid.*D_grid *f0 /c);% ע��˴���λѹ���ಹ���˸�4*pi*R0*f0/c����λ
Srd_ac = Srcmc.*Haz;
%%
eta0 = -center_R0 / cos(theta_rc)*sin(theta_rc)/Vr; %�����ĵ��Ӧ����Բ������Ĵ�Խʱ�̣�
Srd_ac = Srd_ac.*exp(-1j*2*pi*f_eta_grid*eta0);
img_rd = ifft(Srd_ac);

figure; % ���Ƶ�б�ӽ�����¾���ѹ���ҷ�λѹ�����ź�
subplot(121);imagesc(real(img_rd));xlabel('������ʱ�䣨�����㣩');ylabel('��λ��ʱ�䣨�����㣩');title('(a)ʵ��');
subplot(122);imagesc(abs(img_rd));xlabel('������ʱ�䣨�����㣩');title('(b)����');
suptitle('21.9��б�ӽǾ���ѹ���ҷ�λѹ������źţ�ʱ��');

%% 4. ��Ŀ�����
% ����ÿ�������λ�õ�����ֵ
delta_r=delta_r * cos(theta_rc);
ns = round(rs/(delta_r)) + (Nrg/2 + 1);
ms = round(as/delta_a) + (Naz/2 + 1);
len = 16;
p = 1;
target = img_rd(ms(p)-len/2:ms(p)+len/2-1, ns(p)-len/2:ns(p)+len/2-1);
[image_upsample,signal_r,quality_r,signal_a,quality_a] = f_point_analyse(target,delta_r,delta_a);
BW_r= abs(Kr*Tr);
La = 0.886 * 2 * Vs * cos(theta_rc) / BW_dop;   % ���߿׾�����
IRW_r_theory = c/2/BW_r*0.886*1.18;
IRW_a_theory = La/2*Vg/Vs*1.185;
disp(['���������۷ֱ���:',num2str(IRW_r_theory),'m']);
disp(['��λ�����۷ֱ���:',num2str(IRW_a_theory),'m']);
