function [ S ] = azimuth_space_variant( moco_file, S,... 
    lambda, f_etac, Vr, subaperture_num,...
    near_range, range_sample_rate, PRF,...
    azimuth_size, last_pulse_count)
%RANGE_SPACE_VARIANT ��λ��ձ䴦��
%   moco_file �˶�����Ԫ�����ļ�������ƽ̨�켣��Ϣ��
%   S ����λ�ձ���λ�����ľ�����������ź�
%   lambda ����
%   f_etac ����������Ƶ��
%   Vr �״�ƽ̨�˶��ٶ�
%   subaperture_num �ӿ׾�����
%   near_range ��������б��
%   ref_range һ����λ�˶�����ʱѡ��Ĳο�����
%   range_sample_rate �����������
%   PRF ��λ������
%   azimuth_size ��������켣��ϵ��ܵĻز�����
%   last_pulse_count �ϴ��Ѿ���������������������׷��ʽ����

c = 299792458;
Fr = range_sample_rate;
Fa = PRF;
delta_r = c/2/Fr;                           % ������������
[Naz, Nrg] = size(S);
pulse_count = Naz;

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
R0 = near_range + (0:Nrg-1)*delta_r;
R1 = sqrt((sqrt(R0(:).^2-href^2) - (yi-yi_ideal)).^2 + zi.^2);
R1_grid = repmat(reshape(R1, Naz, 1), 1, Nrg);

f_eta = (0:Fa/Naz:(Naz-1)*Fa/Naz).';                   
f_eta = (round((f_etac-f_eta)/Fa)*Fa+f_eta); 
sub_width = ceil(Naz / subaperture_num);    % ÿ���ӿ׾���Ӧ��Ƶ����
for i = 1:subaperture_num
    id_start = (i-1)*sub_width+1;
    id_end = min(Naz, i*sub_width);
    SUB = S(id_start:id_end, :);               % �ӿ׾���ӦƵ���ź�
    sub_fc = f_eta(round((id_start+id_end)/2)); % �ӿ׾���Ӧ����Ƶ��
    sub_theta = asin(sub_fc*lambda/2/Vr);     % ���ӿ׾���ʱ���Ӧб�ӽ�
    sub = ifft(SUB, Naz, 1);        %  �ӿ׾���Ӧʱ���ź�
    R0 = near_range + (0:Nrg-1)*delta_r;
    sub_L = R0 * tan(sub_theta);    % ���ӿ׾���ʱ������Ӧ��λλ�ú����������ľ���
    sub_L_grid = repmat(sub_L, Naz, 1);
    delta_Ra = sqrt(R1_grid.^2 + sub_L_grid.^2) -...
        sqrt(repmat(R0, Naz,1).^2+sub_L_grid.^2); % ������ӿ׾�����Ӧ��λλ����
                                                  %����������������Ӧб�����
    delta_Ra = delta_Ra - (R1_grid -repmat(R0, Naz,1)); % ��ȥ�ο�λ�õ�б�����
    sub = sub.*exp(1j*4*pi*delta_Ra/lambda);    % �Ը��ӿ׾�������λ����
    SUB = fft(sub);
    S(id_start:id_end, :) = SUB(id_start:id_end, :);
end
end

