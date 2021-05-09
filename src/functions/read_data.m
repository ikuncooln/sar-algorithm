function [ s0 ] = read_data( file_name, range_size, x0, y0, height, width)
%READ_DATA ���ļ��ж�ȡSARԭʼ�����ݣ��ļ������ݰ���ʱ�����ȴ洢�����ȴ洢��ĳһ
%   ��λʱ���������������лز��źź��ٽ��Ŵ洢��һ����ز�.ÿһ����IQ�����㰴ʵ��
%   ����ķ�ʽ�洢��ʵ���鲿����4λfloat����.
%
%   file_name �ļ���
%   range_size ÿһ����������Ӧ�ز��źŸ���������������
%   x0 ��ȡ������ʼλ�õ��к����꣨��1��ʼ��
%   y0 ��ȡ������ʼλ�õ��������꣨��1��ʼ��
%   height ��ȡ���ݵĸ߶ȣ��������㣩
%   width ��ȡ���ݵĿ�ȣ��������㣩
%
%   s0 ���ź����ݾ���

%   ����
%       s0 = read_data('data_after_moco.dat', 16384, 1, 1, 4096, 4096)

W = range_size; %  width (pixel) of raw_data(IQ)
s0 = zeros(height, width);  % store the signal

fid = fopen(file_name, 'r');
if fid==-1
    disp(['�ļ�', file_name, '������'])
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

