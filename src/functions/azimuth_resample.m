function [  ] = azimuth_resample( moco_file, in_file, out_file,...
    range_size, azimuth_size, pulse_count, last_pulse_count, MAX_MEM_GB)
%AZIMUTH_RESAMPLE ��λ���ز���
%   moco_file �˶�����Ԫ�����ļ�������ƽ̨�켣��Ϣ��
%   in_file ������SAR�����ļ���
%   out_file ����������ļ���
%   range_size ÿ���ز��ź��ܵĲ�������
%   azimuth_size ��������켣��ϵ��ܵĻز�����
%   pulse_count ���δ����������
%   last_pulse_count �ϴ��Ѿ���������������������׷��ʽ����
%   MAX_MEM_GB �ڴ��С���ƣ�GB��
%   
% parameters convert
MAX_MEM_GB = MAX_MEM_GB * 1024*1024*1024;   % 1�ν���������*GB����                          % ������������
batch_size = floor(MAX_MEM_GB / 128 / range_size);    % ÿ�δ������������
iter = ceil(pulse_count/batch_size);        % ѭ������

[ ~, MOCO_UNIT ] = read_mocodata( moco_file,azimuth_size );
azimuth_pos = MOCO_UNIT.forward;
azimuth_pos_ideal = linspace(azimuth_pos(1), azimuth_pos(end), azimuth_size);

%% �˶�����
current_pulse_count = last_pulse_count+1;
overlap = 20;    % Ϊ�˲�ֵ�������ԣ�����һ��
write_data([], out_file, 0);    % ����ջ������ļ�
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
            disp(['��λ���ز����У�', num2str(col/range_size*100/iter + (k-1)*100/iter), '%']);
        end
    end
    tmp = current_pulse_count-y_start+1;
%     disp(size(s0));
%     disp(['tmp,bs:', num2str(tmp), '-', num2str(bs)]);
    write_data(s0(tmp:tmp+bs-1, :), out_file, 1);
    current_pulse_count = current_pulse_count + bs;
end
end

