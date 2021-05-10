clc;
clear;
close all;

load s_rd_20480_16384_win.mat s
[GMG,LS,Dynamic_range,EVA,Mean,Var] = ImageEvaluation(s);