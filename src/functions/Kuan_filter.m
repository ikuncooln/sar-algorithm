function [ img_Kuan ] = Kuan_filter( I, L )
%KUAN_FILTER ��ͼ�����Kuan�˲�
% I ����Ҷ�ͼ��[0,255]
% L ��Ч����
% img_Lee Lee�˲���ͼ��
% Ŀǰ�˺�����Ϊ����ͼ��߽���������չ�ģ���Ҫ��һ�����ƣ��ο�imfilter������
% ��Lee���ƣ��������� SAR ͼ��ĵ�Ч����ȷ����ͬ��Ȩ�غ��������˲�
% ���ߵ�����ģ��ɼ�������ģʽ
[height, width] = size(I);
I = double(I);          % �������ʱ�������
N = 7;                  % ���򴰿ڴ�СN*N
N_half = floor(N/2);    % �����뾶
N_size = N*N;

% prepare params
Cz2 = 1 / L;

img_Kuan = zeros(height, width);
% filtering
for m = 1:height
    for n = 1:width
        mm = -(N_half-1):N_half + m;
        nn = -(N_half-1):N_half + n;
        mm = mod(mm, height);
        nn = mod(nn, width);
        mm(mm==0) = height;
        nn(nn==0) = width;
        locval = I(mm, nn);
        u = sum(locval(:)) / N_size;            % mean
        v = sum(sum(locval.*locval)) / N_size - u*u;   % variance
        if(u == 0)
            Cy2 = 1;
        else
            Cy2 = v / (u*u);
        end
        k = (1-Cz2/Cy2)/(1+Cz2);
        % caculate current filtered pixel value
        img_Kuan(m,n) = u + k*(I(m,n)-u);
    end
end

img_Kuan = uint8(img_Kuan);
end


