function [  ] = range_resample( moco_file, ...
    in_file, out_file, wave_length, near_range, range_sample_rate,...
    range_size, azimuth_size, pulse_count, last_pulse_count, MAX_MEM_GB)
%RANGE_RESAMPLE �������ز���
%   moco_file �˶�����Ԫ�����ļ�������ƽ̨�켣��Ϣ��
%   in_file ������SAR�����ļ���
%   out_file ����������ļ���
%   wave_length ���Ĳ���
%   near_range ��������б��
%   range_sample_rate �����������
%   range_size ÿ���ز��ź��ܵĲ�������
%   azimuth_size ��������켣��ϵ��ܵĻز�����
%   pulse_count ���δ����������
%   last_pulse_count �ϴ��Ѿ���������������������׷��ʽ����
%   MAX_MEM_GB �ڴ��С���ƣ�GB��
%   
%   see fopec_rr() for more efficiently!
c = 299792458;
% parameters convert
MAX_MEM_GB = MAX_MEM_GB * 1024*1024*1024;   % 1�ν���������*GB����
lambda = wave_length;
f0 = c/lambda;
Fr = range_sample_rate;
Nrg = range_size;
delta_r = c/2/Fr;                           % ������������
batch_size = floor(MAX_MEM_GB / 128 / range_size);    % ÿ�δ������������
iter = ceil(pulse_count/batch_size);        % ѭ������
f_tau = ifftshift((-Nrg/2:Nrg/2-1)*Fr/Nrg);

[ ~, MOCO_UNIT ] = read_mocodata( moco_file,azimuth_size );
xi = reshape(MOCO_UNIT.forward, [], 1);
yi = reshape(MOCO_UNIT.cross, [], 1);
zi = reshape(MOCO_UNIT.height, [], 1);
% ��ϳ����뺽��
p = polyfit(xi, yi, 1);
yi_ideal = polyval(p, xi);
href = mean(zi(:)); % ����ο��߶�
% ���豣�����δ�������Ϣ
idx = last_pulse_count+1:pulse_count;
yi = yi(idx); zi = zi(idx);
yi_ideal = yi_ideal(idx);

% ����б�����
R0 = near_range + (0:Nrg-1)*delta_r;
cos_alpha = href ./ R0;
sin_alpha = sqrt(1-cos_alpha.^2);
[cos_alpha_mtx, delta_y] = meshgrid(cos_alpha, yi-yi_ideal);
[sin_alpha_mtx, delta_z] = meshgrid(sin_alpha, zi-href);

delta_R = delta_z .* cos_alpha_mtx;
clear delta_z cos_alpha_mtx;
delta_R = delta_R + delta_y .* sin_alpha_mtx;
clear delta_y sin_alpha_mtx

%% �˶�����
current_pulse_count = last_pulse_count+1;
write_data([], out_file, 0);    % ����ջ������ļ�

for k = 1:iter
    bs = (last_pulse_count+pulse_count - current_pulse_count)+1;
    if bs > batch_size
        bs = batch_size;
    end
    
    s0 = read_data(in_file, range_size, 1, current_pulse_count,...
    bs, range_size);
    current_pulse_count = current_pulse_count + bs;
    S0 = fft(s0.').';
    f = repmat(f0+f_tau, bs, 1); 
    hc = exp(1j * 4 * pi * ... % Ҫȡ����ο�����λ�õ�delta_R
        repmat(delta_R((k-1)*batch_size+(1:bs), Nrg/2), 1, Nrg) .* f / c);
    s1 = ifft((S0 .* hc).').';
    disp(['�������ز����У�', num2str(k/iter*100), '%']);
    write_data(s1, out_file, 1);
end

end

