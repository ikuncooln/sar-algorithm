function [ ] = write_data(data, file_name, append)
%WRITE_DATA 将SAR复数据写入文件，文件中数据按快时间优先存储，即先存储完某一
%   方位时刻所发射脉冲所有回波信号后再接着存储下一脉冲回波.每一个复IQ采样点按实虚
%   交替的方式存储，实部虚部均是4位float类型.
%
%   包含有SAR复数据的数组，每一行为一个脉冲回波信号
%   file_name 文件名
%   append 1表示以追加方式写入，其他值表示以新建方式写入
if append == 1
    [fid, msg]= fopen(file_name, 'a');
else
    [fid, msg] = fopen(file_name, 'w');
end
if fid == -1
    error(msg);
end

data = reshape(data.', [], 1);
N = length(data);
data_row = zeros(2*N, 1);
data_row(1:2:end) = real(data);
data_row(2:2:end) = imag(data);
clear data
fwrite(fid, data_row, 'float32');
fclose(fid);
clear data_row;
end

