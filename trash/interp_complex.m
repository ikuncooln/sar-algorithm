clear all;
close all;
K = 50;
Tr = 1;
Fr = 100;
N = Tr*Fr;

t = (-N/2:N/2-1)/Fr;
t = t(1:5:end);
s0 = exp(1j*pi*K*t.^2);

tq = linspace(min(t(:)), max(t(:)), N*10);
sq = interp1(t, s0, tq);

figure;
plot(imag(s0), 'o');
title('s0');

figure;
plot(imag(sq), 'o');
title('sq');