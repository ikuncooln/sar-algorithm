function [ delta_R ] = range_space_variant( moco_file,...
    near_range, ref_range, range_sample_rate, range_size, azimuth_size,... 
    width, range_start, pulse_count, last_pulse_count)
%RANGE_SPACE_VARIANT 计算距离向相对于参考位置的空变距离误差
%   moco_file 运动补偿元数据文件名（即平台轨迹信息）
%   near_range 最近点最近斜距
%   ref_range 一阶相位运动补偿时选择的参考距离
%   range_sample_rate 距离向采样率
%   range_size 每个回波信号总的采样点数
%   azimuth_size 用于理想轨迹拟合的总的回波个数
%   width 计算相对距离误差区域的距离向宽度（采样点数）
%   range_start 计算相对距离误差区域的起始距离向采样点
%   pulse_count 本次处理脉冲个数
%   last_pulse_count 上次已经处理过的脉冲个数（用于追加式处理）
c = 299792458;
Fr = range_sample_rate;
Nrg = range_size;
delta_r = c/2/Fr;                           % 距离向采样间距

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

N0 = round((ref_range - near_range) / delta_r) + 1;
assert(N0<=Nrg);
delta_R = delta_R - repmat(delta_R(:, N0), 1, range_size);
delta_R = delta_R(:, range_start+1:range_start+width);

end

