function [image_upsample,signal_r,quality_r,signal_a,quality_a] = f_point_analyse(target,delta_r,delta_a,ratio,show)
%    输入点目标切片，输出升采样后的结果以及各参数
%输入参数：
%    target是待分析的点目标切片（二维复数矩阵）
%    delta_r、delta_a分别是距离向、方位向采样间距（单位为m）
%    ratio是升采样倍数，缺省值是16
%    show为1时表示画出升采样中间步骤图，其他值则不画出，缺省为0
%返回值：
%    image_upsample是升采样后的结果（二维复数矩阵）
%    signal_r是距离向切片（一维复数向量）
%    quality_r为该切片对应的质量参数，3个值分别为IRW（单位为m）、PSLR、ISLR
%    signal_a是方位向切片（一维复数向量）
%    quality_a为该切片对应的质量参数，3个值分别为IRW（单位为m）、PSLR、ISLR
if ~exist('ratio', 'var'),
    ratio = 16;
end
if ~exist('show', 'var'),
    show = 0;
end
[size_a,size_r] = size(target);
% 原始二维频谱
target_ff = fft2(target);

% 频谱搬移
step = 0.1;
k = 1.1;
while(1)  
    % 找到频谱重心
    bw = repmat(abs(target_ff)>=max(abs(target_ff(:)))/k,3,3);
    [L,num] = bwlabel(bw,4);
    max_area = 0;
    for n = 1:num
        if(length(find(L==n)) > max_area);
            max_n = n;
            max_area = length(find(L==n));
        end
    end
    [az,rg]=find(L==max_n);
    max_a0 = round(mean(az));
    max_r0 = round(mean(rg));     
    % 将频谱搬移至类似P38的形式，希望搬移之后有9个连通区域
    target_ff_shift = circshift(target_ff,-max_a0,1);
    target_ff_shift = circshift(target_ff_shift,-max_r0,2);
    four = repmat(target_ff_shift,2,2); % 正常情况下，搬移之后有9个连通区域   
    test = (abs(four) >= max(abs(four(:)))/k);
    [~,num_block] = bwlabel(test,4);
    if(num_block ==9)
        break;
    else
        k = k+step;
    end
end

pos_a = get_interval(target_ff_shift);
pos_r = get_interval(target_ff_shift.');

% im代表找到的边界点
im = zeros(2*size_a,2*size_r);
for az = 1:2*size_a
    for rg = 1:2*size_r
        im(az,pos_r(ceil(mod(az-0.1,size_a)))) = 1;
        im(pos_a(ceil(mod(rg-0.1,size_r))),rg)=1;
        im(az,size_r+pos_r(ceil(mod(az-0.1,size_a)))) = 1;
        im(size_a+pos_a(ceil(mod(rg-0.1,size_r))),rg)=1;
    end
end
% 将边界点连接起来
[~,num] = bwlabel(im,4);
while(num >1)
    im = bwmorph(im,'dilate');
    [~,num] = bwlabel(im,4);
end
boundary = bwmorph(im,'thin',Inf);

[L_final,~] = bwlabel(1-boundary,4);
L_final = bwmorph(L_final==L_final(size_a+1,size_r+1),'close',Inf);

L_all = L_final*1 + fftshift(L_final,1)*2 + fftshift(L_final,2)*3....
        + fftshift(L_final)*4;
L_all = circshift(L_all,max_a0,1);
L_all = circshift(L_all,max_r0,2);  
chose = [L_all(1,1),L_all(2*size_a,1),L_all(1,2*size_r),L_all(2*size_a,2*size_r)];
chose = chose(chose~=0);
tmp = repmat(target_ff,2,2).*(L_all == mode(chose));

%% 补零
N_fft_a = 2^nextpow2(ratio*size_a);
N_fft_r = 2^nextpow2(ratio*size_r);

target_ff_up = zeros(N_fft_a,N_fft_r);
target_ff_up(1:size_a,1:size_r) = tmp(1:size_a,1:size_r);
target_ff_up(N_fft_a-size_a+1:N_fft_a,1:size_r) = tmp(1+size_a:2*size_a,1:size_r);
target_ff_up(1:size_a,N_fft_r-size_r+1:N_fft_r) = tmp(1:size_a,1+size_r:2*size_r);
target_ff_up(N_fft_a-size_a+1:N_fft_a,N_fft_r-size_r+1:N_fft_r) = tmp(1+size_a:2*size_a,1+size_r:2*size_r);

%% ifft
image_upsample = ifft2(target_ff_up);
if (show == 1)
    figure;subplot(231);imagesc(abs(target_ff));title('(a)原始频谱');
    subplot(232);imagesc(abs(repmat(target_ff,2,2)));title('(b)复制后的频谱');
    subplot(233);imagesc(abs(four));title('(c)搬移、复制后的频谱');
    subplot(234);imagesc(abs(four.*L_final));title('(d)取出的频谱');
    subplot(235);imagesc(abs(tmp));title('(e)补零结果示意');
    subplot(236);imagesc(abs(image_upsample));title('(f)升采样后的时域');
end
%% 距离和方位剖面（直接截取水平和垂直剖面）
% [az_cut,rg_cut] = find(abs(image_upsample) == max(abs(image_upsample(:))),1); % 最大值点
% signal_r = image_upsample(az_cut,:); % 距离剖面
% quality_r = f_IRW_PSLR_ISLR(signal_r);
% quality_r(1) = quality_r(1)*size_r/N_fft_r*delta_r_new;
% signal_a = image_upsample(:,rg_cut); % 方位剖面
% quality_a = f_IRW_PSLR_ISLR(signal_a);
% quality_a(1) = quality_a(1)*size_a/N_fft_a*delta_a_new;
%% 距离和方位剖面（鼠标选取）
figure;imagesc(abs(image_upsample));
[max_az,max_rg] = find(abs(image_upsample) == max(abs(image_upsample(:))),1);
hold on;plot(max_rg,max_az,'r.');
% 第一下选择距离向剖面，第二下选择方位向剖面
[loc_rg,loc_az] = ginput(2);
% 距离向
k_r = (loc_az(1)-max_az)/(loc_rg(1)-max_rg);
y_r = max_az + round(k_r*((1:N_fft_r)-max_rg));
signal_r = image_upsample(sub2ind([N_fft_a,N_fft_r],y_r,1:N_fft_r));
quality_r = f_IRW_PSLR_ISLR(signal_r);
delta_r_new = sqrt(delta_r^2+k_r^2*delta_a^2);
quality_r(1) = quality_r(1)*size_r/N_fft_r*delta_r_new;
% 方位向
k_a = (loc_rg(2)-max_rg)/(loc_az(2)-max_az);
x_a = max_rg + round(k_a*((1:N_fft_a)-max_az));
signal_a = image_upsample(sub2ind([N_fft_a,N_fft_r],1:N_fft_a,x_a));
quality_a = f_IRW_PSLR_ISLR(signal_a);
delta_a_new = sqrt(delta_a^2+k_a^2*delta_r^2);
quality_a(1) = quality_a(1)*size_a/N_fft_a*delta_a_new;
% 画线示意
plot(1:N_fft_r,y_r,'r',x_a,1:N_fft_a,'r');
hold off;
%% 画图
xticks = 0:size_r/N_fft_r:size_r-size_r/N_fft_r;
yticks = 0:size_a/N_fft_a:size_a-size_a/N_fft_a;
figure;
subplot(321);imagesc(xticks,yticks,abs(image_upsample));title('(a)放大后的点目标');xlabel('距离向（采样点）');ylabel('方位向（采样点）');
hold on;plot(((1:N_fft_r)-1)*size_r/N_fft_r,(y_r-1)*size_r/N_fft_r,'r',(x_a-1)*size_a/N_fft_a,((1:N_fft_a)-1)*size_a/N_fft_a,'r');hold off;
subplot(322);contour(xticks,yticks,abs(image_upsample),30);title('(b)放大后的点目标等值线图');set(gca,'YDir','reverse');xlabel('距离向（采样点）');ylabel('方位向（采样点）');
subplot(323);plot(xticks*delta_r_new,20*log10(abs(signal_r)/max(abs(signal_r(:)))));xlim([0,size_r*delta_r_new]);ylim([-35,0]);title('(c)距离剖面图');xlabel('距离向(m)');ylabel('幅度(dB)');
subplot(324);plot(yticks*delta_a_new,20*log10(abs(signal_a)/max(abs(signal_a(:)))));xlim([0,size_a*delta_a_new]);ylim([-35,0]);title('(d)方位剖面图');xlabel('方位向(m)');ylabel('幅度(dB)');
subplot(325);plot(xticks*delta_r_new,angle(signal_r)*180/pi);xlim([0,size_r*delta_r_new]);ylim([-200,200]);title('(e)距离相位');xlabel('距离向(m)');ylabel('相位(°)');
subplot(326);plot(yticks*delta_a_new,angle(signal_a)*180/pi);xlim([0,size_a*delta_a_new]);ylim([-200,200]);title('(f)方位相位');xlabel('方位向(m)');ylabel('相位(°)');
end