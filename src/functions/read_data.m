function [ s0 ] = read_data( file_name, range_size, x0, y0, height, width)
%READ_DATA 从文件中读取SAR原始复数据，文件中数据按快时间优先存储，即先存储完某一
%   方位时刻所发射脉冲所有回波信号后再接着存储下一脉冲回波.每一个复IQ采样点按实虚
%   交替的方式存储，实部虚部均是4位float类型.
%
%   file_name 文件名
%   range_size 每一个脉冲所对应回波信号复采样点数（幅宽）
%   x0 读取数据起始位置的行号坐标（从1开始）
%   y0 读取数据起始位置的列列坐标（从1开始）
%   height 读取数据的高度（复采样点）
%   width 读取数据的宽度（复采样点）
%
%   s0 复信号数据矩阵

%   例：
%       s0 = read_data('data_after_moco.dat', 16384, 1, 1, 4096, 4096)

W = range_size; %  width (pixel) of raw_data(IQ)
s0 = zeros(height, width);  % store the signal

fid = fopen(file_name, 'r');
if fid==-1
    disp(['文件', file_name, '不存在'])
    return;
end
offset = (2*(x0-1) + W*2*(y0-1))*4;
fseek(fid, 0, -1);

for r = 1:height
    fseek(fid, offset, 0);
    offset = (W*2 - width*2)*4;  % update offset
    data_row = fread(fid, width*2,'float32');
    data_row = data_row.';
    s0(r,:) = data_row(1:2:end) + 1j*data_row(2:2:end);
end
fclose(fid);
clear data_row;
end

