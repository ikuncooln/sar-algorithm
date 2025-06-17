function [ img_Kuan ] = Kuan_filter( I, L )
%KUAN_FILTER 对图像进行Kuan滤波
% I 输入灰度图像[0,255]
% L 等效视数
% img_Lee Lee滤波后图像
% 目前此函数认为输入图像边界是周期拓展的，需要进一步完善（参考imfilter函数）
% 和Lee类似，它依赖于 SAR 图像的等效视数确定不同的权重函数进行滤波
% 将斑点噪声模拟成加性线性模式
[height, width] = size(I);
I = double(I);          % 避免计算时类型溢出
N = 7;                  % 邻域窗口大小N*N
N_half = floor(N/2);    % 滑窗半径
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


