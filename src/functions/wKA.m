function [ img_wk ] = wKA( s0, lambda, Kr, Vr, Fr, PRF, Rref, f_etac )
%wKA Foucus your SAR data by Omega-K Algorithm
%   
%   s2 the focused image
%   s1 the focused image without stolt map
%   ...
c = 299792458;
f0 = c / lambda;
Fa = PRF;

%% 1. �ο��������
% ����ο�����
[Naz, Nrg] = size(s0);
f_tau = ifftshift(-Nrg/2:Nrg/2-1) * Fr / Nrg;
f_eta = (ifftshift(-Naz/2:Naz/2-1) * Fa / Naz).';
% ��Ƶ��[-Fa/2, Fa/2]ӳ�����ʵ�ʣ�����ǰ����Ӧ��Ƶ��
f_eta = f_eta + round((f_etac - f_eta) / Fa) * Fa;
[f_tau_grid, f_eta_grid] = meshgrid(f_tau, f_eta);

clear f_tau; clear f_eta;

theta_ref = 4*pi*Rref / c * sqrt((f0+f_tau_grid).^2 ...
- c^2*f_eta_grid.^2/(4*Vr^2)) + pi*f_tau_grid.^2/Kr;
Href = exp(1j * theta_ref);
clear theta_ref;
S2df = fft2(s0);
clear s0;
% ��Ϊ�ο������������������λ2*pi*2*Rref/c/D��������Ҫ���䲹��������Ȼʱ��ͷ�����ƽ��
D_fetac_Vref = sqrt(1-c^2*f_etac^2/(4*Vr^2*f0^2));
% һ�²ο�������˵��µľ���ʱ��ͷ�λʱ��ƫ�ƣ���Ҫ��theta_refչ��ʱ�й���f_tau��f_eta�������
tau_shift = (2*Rref/c/D_fetac_Vref);
eta_shift = Rref * c * f_etac / (2 * Vr^2 * f0 * D_fetac_Vref);
clear D_fetac_Vref;

Hshift1 = exp(-1j*2*pi*f_tau_grid*tau_shift) ...
    .*exp(1j*2*pi*f_eta_grid*eta_shift);
clear tau_shift; clear eta_shift;
S_RFM = S2df .* Href .* Hshift1;
clear Href; clear Hshift1; clear S2df; 
% s1 = ifft2(S_RFM);
% figure;
% imagesc(abs(s1)); title('��Stolt��ֵ��ѹ��Ŀ��');

%% 2. Stoltӳ��
% ����ӳ��Ƶ��ƫ����
f_tau1_0 = sqrt((f0 + 0)^2 - c^2*f_eta_grid.^2/(4*Vr^2)) - f0; % ӳ������������Ƶ��
% ��Ƶ��[-Fr/2, Fr/2]ӳ�����ʵ�ʣ�����ǰ����Ӧ��Ƶ�� 
f_tau1_grid = f_tau_grid + round((f_tau1_0 - f_tau_grid)/Fr)*Fr;
OFFSET = sqrt((f_tau1_grid + f0).^2 + c^2*f_eta_grid.^2/(4*Vr^2)) - f0 - f_tau_grid;
OFFSET = OFFSET / (Fr/Nrg);

clear f_eta_grid; clear f_tau1_grid;

% ���������ھ����㶯У�����в�ֵ
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

% פ����λԭ��Ҫ��Chirp�ź�s(t)�ĳ���ʱ��Ϊ[-T/2, T/2]���������ǲ����Ƴ����ϵ�
% stoltӳ�乫ʽ��������ֱ�Ӱ�[-T/2, T/2]��ȡ���ź�ȥ��DFT����ô�õ����ź�Ƶ�����
% s(t)�����������λ�����������Ҫ��ʱ����ѭ����λ����ֱ����Ƶ�����������λ����
% ����ο�ʹ��Matlab��FFT�任ֵ��ע���һЩ����.md
Hshift2 = exp(-1j*2*pi*(f_tau_grid*(-Nrg/2/Fr)));
clear f_tau_grid;
S_RFM = S_RFM .* Hshift2;
% ��ֵӳ��
Sstolt = zeros(Naz, Nrg);  % ���Stoltӳ�����ź�
hwait=waitbar(0,'��ȴ�>>>>>>>>');
for i = 1:Naz
    for j = 1:Nrg
        offset_int = ceil(OFFSET(i,j));
        offset_frac = round((offset_int - OFFSET(i,j)) * 16);
        if offset_frac == 0
            Sstolt(i,j) = S_RFM(i,ceil(mod(j+offset_int-0.1,Nrg)));   % �����ź�����S1�������Լٶ�
        else
            Sstolt(i,j) = S_RFM(i, ceil(mod((j+offset_int-4:j+offset_int+3)-0.1,Nrg))) * hx(offset_frac,:).';
        end
        
    end
    if mod(i,200)== 0
    waitbar(i/Naz,hwait,['Stoltӳ����: ', num2str(i/Naz*100), '%']);
    end
end
close(hwait);
clear OFFSET;

Sstolt = Sstolt ./ Hshift2;  % Ŀ����Ѿ�����λ��R0,ȥ��ǰ�����λЧ��
clear Hshift2;
img_wk = ifft2(Sstolt);
end

