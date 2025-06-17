

% I 输入灰度图像[0,255]
% img_Lee Lee滤波后图像
% 目前此函数认为输入图像边界是周期拓展的，需要进一步完善（参考imfilter函数）

[height, width] = size(I);
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
        total = sum(locval(:));
        u = total / N_size;             % mean
        v = sum(sum(locval.*locval)) - u*u;   % variance
        if(v == 0)
            k = 0;
        else
            k = (v - u*u*Cz2) / (v * (1+Cz2));
        end
        
        % caculate current filtered pixel value
        img_Lee(m,n) = u + k*(I(m,n)-u);
    end
end
        
        