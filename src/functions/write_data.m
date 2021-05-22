function [ ] = write_data(data, file_name, append)
%WRITE_DATA ��SAR������д���ļ����ļ������ݰ���ʱ�����ȴ洢�����ȴ洢��ĳһ
%   ��λʱ���������������лز��źź��ٽ��Ŵ洢��һ����ز�.ÿһ����IQ�����㰴ʵ��
%   ����ķ�ʽ�洢��ʵ���鲿����4λfloat����.
%
%   ������SAR�����ݵ����飬ÿһ��Ϊһ������ز��ź�
%   file_name �ļ���
%   append 1��ʾ��׷�ӷ�ʽд�룬����ֵ��ʾ���½���ʽд��
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

