function [  ] = range_resample( moco_file, ...
    in_file, out_file, wave_length, near_range, range_sample_rate,...
    range_size, azimuth_size, pulse_count, last_pulse_count, MAX_MEM_GB)
%RANGE_RESAMPLE 距离向重采样
%   moco_file 运动补偿元数据文件名（即平台轨迹信息）
%   in_file 待处理SAR数据文件名
%   out_file 处理结果输出文件名
%   wave_length 中心波长
%   near_range 最近点最近斜距
%   range_sample_rate 距离向采样率
%   range_size 每个回波信号总的采样点数
%   azimuth_size 用于理想轨迹拟合的总的回波个数
%   pulse_count 本次处理脉冲个数
%   last_pulse_count 上次已经处理过的脉冲个数（用于追加式处理）
%   MAX_MEM_GB 内存大小限制（GB）
%   
%   see fopec_rr() for more efficiently!
c = 299792458;
% parameters convert
MAX_MEM_GB = MAX_MEM_GB * 1024*1024*1024;   % 1次仅处理不超过*GB数据
lambda = wave_length;
f0 = c/lambda;
Fr = range_sample_rate;
Nrg = range_size;
delta_r = c/2/Fr;                           % 距离向采样间距
batch_size = floor(MAX_MEM_GB / 128 / range_size);    % 每次处理的脉冲数量
iter = ceil(pulse_count/batch_size);        % 循环次数
f_tau = ifftshift((-Nrg/2:Nrg/2-1)*Fr/Nrg);

[ ~, MOCO_UNIT ] = read_mocodata( moco_file,azimuth_size );
xi = reshape(MOCO_UNIT.forward, [], 1);
yi = reshape(MOCO_UNIT.cross, [], 1);
zi = reshape(MOCO_UNIT.height, [], 1);
% 拟合出理想航迹
p = polyfit(xi, yi, 1);
yi_ideal = polyval(p, xi);
href = mean(zi(:)); % 理想参考高度
% 仅需保留本次处理部分信息
idx = last_pulse_count+1:pulse_count;
yi = yi(idx); zi = zi(idx);
yi_ideal = yi_ideal(idx);

% 计算斜距误差
R0 = near_range + (0:Nrg-1)*delta_r;
cos_alpha = href ./ R0;
sin_alpha = sqrt(1-cos_alpha.^2);
[cos_alpha_mtx, delta_y] = meshgrid(cos_alpha, yi-yi_ideal);
[sin_alpha_mtx, delta_z] = meshgrid(sin_alpha, zi-href);

delta_R = delta_z .* cos_alpha_mtx;
clear delta_z cos_alpha_mtx;
delta_R = delta_R + delta_y .* sin_alpha_mtx;
clear delta_y sin_alpha_mtx

%% 运动补偿
current_pulse_count = last_pulse_count+1;
write_data([], out_file, 0);    % 先清空或建立该文件

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
    hc = exp(1j * 4 * pi * ... % 要取距离参考中心位置的delta_R
        repmat(delta_R((k-1)*batch_size+(1:bs), Nrg/2), 1, Nrg) .* f / c);
    s1 = ifft((S0 .* hc).').';
    disp(['距离向重采样中：', num2str(k/iter*100), '%']);
    write_data(s1, out_file, 1);
end

end

