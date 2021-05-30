function [ delta_R ] = range_space_variant( moco_file,...
    near_range, ref_range, range_sample_rate, range_size, azimuth_size,... 
    width, range_start, pulse_count, last_pulse_count)
%RANGE_SPACE_VARIANT �������������ڲο�λ�õĿձ�������
%   moco_file �˶�����Ԫ�����ļ�������ƽ̨�켣��Ϣ��
%   near_range ��������б��
%   ref_range һ����λ�˶�����ʱѡ��Ĳο�����
%   range_sample_rate �����������
%   range_size ÿ���ز��ź��ܵĲ�������
%   azimuth_size ��������켣��ϵ��ܵĻز�����
%   width ������Ծ����������ľ������ȣ�����������
%   range_start ������Ծ�������������ʼ�����������
%   pulse_count ���δ����������
%   last_pulse_count �ϴ��Ѿ���������������������׷��ʽ����
c = 299792458;
Fr = range_sample_rate;
Nrg = range_size;
delta_r = c/2/Fr;                           % ������������

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

N0 = round((ref_range - near_range) / delta_r) + 1;
assert(N0<=Nrg);
delta_R = delta_R - repmat(delta_R(:, N0), 1, range_size);
delta_R = delta_R(:, range_start+1:range_start+width);

end

