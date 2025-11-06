% Example of concatenating multiple aperture consisting of rectangular mathematical elements into a single aperture.
%
% It is assumed that field_init is called before this script is executing.
%
% Author: Zhixuan Chang
% Date: Nov. 6, 2025
% Version: 0.1

clear;
close all;
clc;

%%
aper1 = xdc_piston(2e-3, 0.1e-3);
aper2 = xdc_single_rect([-3e-3, -3e-3, 0], 2e-3, 2e-3, 20, 20, [0, 0, 1000]);
aper3 = xdc_single_rect([ 3e-3, -3e-3, 0], 2e-3, 2e-3, 20, 20, [0, 0, 1000]);
aper4 = xdc_single_rect([ 3e-3,  3e-3, 0], 2e-3, 2e-3, 20, 20, [0, 0, 1000]);
aper5 = xdc_single_rect([-3e-3,  3e-3, 0], 2e-3, 2e-3, 20, 20, [0, 0, 1000]);
aper_set = [aper1, aper2, aper3, aper4, aper5];
focus = [0, 0, 1000];

data = xdc_get(aper2);
disp(unique(data(1,:)));

aper = xdc_concatenate(aper_set, focus);

for i = 1 : length(aper_set)
    xdc_free(aper_set(i));
end

xdc_display_aper(aper);

data = xdc_get(aper);
disp(unique(data(1,:)));

xdc_free(aper);
