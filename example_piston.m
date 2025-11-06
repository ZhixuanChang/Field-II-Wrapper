clear;
close all;
clc;

% field_init;

c = 1540;
fs = 100e6;

set_field('fs', fs);
set_field('c', c);

tx_aper = xdc_piston(2e-3, 0.1e-3);  % 可以生成一个物理阵元中心点在(0,0,0)的圆形阵元，作平移之后即可得到任意位置上的阵元
% 考虑使用xdc_piston生成圆形换能器组成的稀疏阵列

xdc_display_aper(tx_aper);

xdc_free(tx_aper);
