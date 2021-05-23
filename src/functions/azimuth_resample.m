function [  ] = azimuth_resample( moco_file, in_file, out_file,...
    range_size, azimuth_size, pulse_count, last_pulse_count, MAX_MEM_GB)
%AZIMUTH_RESAMPLE 方位向重采样
%   moco_file 运动补偿元数据文件名（即平台轨迹信息）
%   in_file 待处理SAR数据文件名
%   out_file 处理结果输出文件名
%   range_size 每个回波信号总的采样点数
%   azimuth_size 用于理想轨迹拟合的总的回波个数
%   pulse_count 本次处理脉冲个数
%   last_pulse_count 上次已经处理过的脉冲个数（用于追加式处理）
%   MAX_MEM_GB 内存大小限制（GB）
%   
% parameters convert
MAX_MEM_GB = MAX_MEM_GB * 1024*1024*1024;   % 1次仅处理不超过*GB数据                          % 距离向采样间距
batch_size = floor(MAX_MEM_GB / 128 / range_size);    % 每次处理的脉冲数量
iter = ceil(pulse_count/batch_size);        % 循环次数

[ ~, MOCO_UNIT ] = read_mocodata( moco_file,azimuth_size );
azimuth_pos = MOCO_UNIT.forward;
azimuth_pos_ideal = linspace(azimuth_pos(1), azimuth_pos(end), azimuth_size);

%% 运动补偿
current_pulse_count = last_pulse_count+1;
overlap = 20;    % 为了插值的连续性，交叠一行
write_data([], out_file, 0);    % 先清空或建立该文件
for k = 1:iter
    bs = (last_pulse_count+pulse_count - current_pulse_count)+1;
    if bs > batch_size
        bs = batch_size;
    end

    y_start = current_pulse_count - overlap;
    y_end = current_pulse_count + bs-1 + overlap;
    if y_start <= 0
        y_start = 1;   
    end
    if y_end > last_pulse_count + pulse_count
        y_end = last_pulse_count + pulse_count;
    end
    
    s0 = read_data(in_file, range_size, 1, y_start,...
    y_end-y_start+1, range_size);

    assert(y_end-y_start+1>=bs);
    for col = 1:range_size
        tmp = interp1(azimuth_pos(y_start:y_end), s0(:, col),...
         azimuth_pos_ideal(y_start:y_end), 'linear', 'extrap');
        assert(sum(isnan(tmp)) == 0);
        s0(:, col) = tmp;
        if mod(col, 1000) == 0
            disp(['方位向重采样中：', num2str(col/range_size*100/iter + (k-1)*100/iter), '%']);
        end
    end
    tmp = current_pulse_count-y_start+1;
%     disp(size(s0));
%     disp(['tmp,bs:', num2str(tmp), '-', num2str(bs)]);
    write_data(s0(tmp:tmp+bs-1, :), out_file, 1);
    current_pulse_count = current_pulse_count + bs;
end
end

