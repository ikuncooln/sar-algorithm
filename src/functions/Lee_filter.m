function [ img_Lee ] = Lee_filter( I )
%LEE_FILTER 对图像进行Lee滤波
% I 输入灰度图像[0,255]
% img_Lee Lee滤波后图像
% 目前此函数认为输入图像边界是周期拓展的，需要进一步完善（参考imfilter函数）

[height, width] = size(I);
I = double(I);          % 避免计算时类型溢出
N = 7;                  % 邻域窗口大小N*N
N_half = floor(N/2);    % 滑窗半径
N_size = N*N;

% prepare params
gu = sum(I(:)) / (height*width);        % 全局均值global mean
gv = sum(sum(I.*I)) / (height*width) - gu*gu;   % 全局方差global variance
if gu == 0
    Cz2 = 1;
else
    Cz2 = gv / (gu*gu);
end

img_Lee = zeros(height, width);
slice = -(N_half-1):N_half;
mbound = height-N_half;
nbound = width-N_half;
% filtering
hwait=waitbar(0,'Lee filtering, please wait>>>>>>>>');
for m = 1:height
    for n = 1:width
        mm = slice + m;
        nn = slice + n;
        if m < N_half || m > mbound
            mm = mod(mm, height);
            mm(mm==0) = height;
        end
        if n < N_half || n > nbound
            nn = mod(nn, width);
            nn(nn==0) = width;
        end
        locval = I(mm, nn);
        u = mean(locval(:));    % mean
        sqrtv = std(locval(:)); 
        v = sqrtv * sqrtv;
        if(v == 0)
            k = 1;
        else
            k = (v - u*u*Cz2) / (v * (1+Cz2));
        end
        
        % caculate current filtered pixel value
        img_Lee(m,n) = u + k*(I(m,n)-u);
    end
    waitbar(m/height,hwait,['Lee filtering: ', num2str(m/height*100), '%']);
end
close(hwait);

img_Lee = uint8(img_Lee);
end

