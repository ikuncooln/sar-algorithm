function pos = get_interval(target_ff)
% 得到频谱间隔,pos是从1-16的顺序
    [size_a,size_r] = size(target_ff);
    
    % 如何找到一条直线，使得其经过的像素的像素值之和最小?
    % 穷举法（严格来说，没有找到全部直线，但应该足够用了）
    x=1:size_r;
    center = [size_r/2+1,size_a/2+1];
    k=linspace((size_a-center(2))/(1-center(1)),size_a/size_r,100);
    [x,k] = meshgrid(x,k);
    y = round(k.*(x-center(1))+center(2));
    classNo = unique(y,'rows');
    N_all = size(classNo,1);
    
    target_ff_abs = abs(target_ff);
    % 距离向补零位置
    best_sum = Inf;
    for az = 1:size_a
        for rg = 1:size_r
            for line = 1:N_all                
                idx = sub2ind(size(target_ff_abs),ceil(mod(classNo(line,:)+az-center(2)-0.1,size_a)),ceil(mod((1:size_r)+rg-center(1)-0.1,size_r)));
                sum_up = sum(target_ff_abs(idx));
                if(sum_up<best_sum)
                    best_sum = sum_up;
                    [loc_a,loc_r] = ind2sub(size(target_ff_abs),idx);
                end
            end
        end
    end
    
    [~,I] = sort(loc_r);
    pos = loc_a(I);
     
end