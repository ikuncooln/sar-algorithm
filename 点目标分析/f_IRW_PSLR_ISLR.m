function quality = f_IRW_PSLR_ISLR(signal)
%    计算冲激响应宽度、峰值旁瓣比和积分旁瓣比
%    signal是待分析的冲激响应
%    quality的3个值分别代表IRW（采样点）、PSLR、ISLR
quality = zeros(1,3);
signal_dB = 20*log10(abs(signal)/max(abs(signal(:))));
% IRW
quality(1) = sum(signal_dB >= -3);
% PSLR
signal_abs = abs(signal);
[pks,locs]= findpeaks(signal_abs,'MinPeakWidth',quality(1)/4);
main = locs(find(pks == max(pks),1));
second1 = locs(find(pks == max(pks),1)-1);
second2 = locs(find(pks == max(pks),1)+1);
quality(2) = max([signal_dB(second1),signal_dB(second2)]);
% ISLR
zero1 = find(signal_abs(second1:main) == min(signal_abs(second1:main)),1)+second1-1;
zero2 = find(signal_abs(main:second2) == min(signal_abs(main:second2)),1)+main-1;
P_main = sum(signal_abs(zero1:zero2).^2);
P_total = sum(signal_abs.^2);
quality(3) = 10*log10((P_total-P_main)/P_main);
end