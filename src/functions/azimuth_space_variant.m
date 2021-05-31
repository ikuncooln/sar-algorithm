function [ S ] = azimuth_space_variant( moco_file, S,... 
    lambda, f_etac, Vr, subaperture_num,...
    near_range, range_sample_rate, PRF,...
    azimuth_size, last_pulse_count)
%RANGE_SPACE_VARIANT 方位向空变处理
%   moco_file 运动补偿元数据文件名（即平台轨迹信息）
%   S 待方位空变相位补偿的距离多普勒域信号
%   lambda 波长
%   f_etac 多普勒中心频率
%   Vr 雷达平台运动速度
%   subaperture_num 子孔径数量
%   near_range 最近点最近斜距
%   ref_range 一阶相位运动补偿时选择的参考距离
%   range_sample_rate 距离向采样率
%   PRF 方位采样率
%   azimuth_size 用于理想轨迹拟合的总的回波个数
%   last_pulse_count 上次已经处理过的脉冲个数（用于追加式处理）

c = 299792458;
Fr = range_sample_rate;
Fa = PRF;
delta_r = c/2/Fr;                           % 距离向采样间距
[Naz, Nrg] = size(S);
pulse_count = Naz;

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
R0 = near_range + (0:Nrg-1)*delta_r;
R1 = sqrt((sqrt(R0(:).^2-href^2) - (yi-yi_ideal)).^2 + zi.^2);
R1_grid = repmat(reshape(R1, Naz, 1), 1, Nrg);

f_eta = (0:Fa/Naz:(Naz-1)*Fa/Naz).';                   
f_eta = (round((f_etac-f_eta)/Fa)*Fa+f_eta); 
sub_width = ceil(Naz / subaperture_num);    % 每个子孔径对应的频域宽度
for i = 1:subaperture_num
    id_start = (i-1)*sub_width+1;
    id_end = min(Naz, i*sub_width);
    SUB = S(id_start:id_end, :);               % 子孔径对应频域信号
    sub_fc = f_eta(round((id_start+id_end)/2)); % 子孔径对应中心频率
    sub_theta = asin(sub_fc*lambda/2/Vr);     % 该子孔径在时域对应斜视角
    sub = ifft(SUB, Naz, 1);        %  子孔径对应时域信号
    R0 = near_range + (0:Nrg-1)*delta_r;
    sub_L = R0 * tan(sub_theta);    % 该子孔径在时域所对应方位位置和零多普勒面的距离
    sub_L_grid = repmat(sub_L, Naz, 1);
    delta_Ra = sqrt(R1_grid.^2 + sub_L_grid.^2) -...
        sqrt(repmat(R0, Naz,1).^2+sub_L_grid.^2); % 计算该子孔径所对应方位位置上
                                                  %各个距离门上所对应斜距误差
    delta_Ra = delta_Ra - (R1_grid -repmat(R0, Naz,1)); % 减去参考位置的斜距误差
    sub = sub.*exp(1j*4*pi*delta_Ra/lambda);    % 对该子孔径进行相位补偿
    SUB = fft(sub);
    S(id_start:id_end, :) = SUB(id_start:id_end, :);
end
end

