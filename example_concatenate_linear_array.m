clear;
close all;
clc;

f0              = 5e6;              % Transducer center frequency [Hz]
pitch           = 0.5e-3;           % Pitch [m]
width           = 0.45e-3;          % Width of element
height          = 5e-3;             % Height of element [m]
num_elem        = 64;               % Number of elements in the transducer
num_sub_w       = 5;                % number of sub-elements along width direction
num_sub_h       = 10;               % number of sub-elements along height direction
center          = [0, 0, 0];        % center of the array
angle           = 0;                % angle between the array and the x-axis, [rad]
focus           = [0 0 1000];       % Fixed focal point [m]

aper_ref = xdc_matrix_array(num_elem,pitch,width,1,height,height,num_sub_w,num_sub_h,focus);

aper_set = zeros(1, num_elem);
elem_x = (-(num_elem-1)/2 : (num_elem-1)/2) * pitch;
elem_y = zeros(1, num_elem);
elem_z = zeros(1, num_elem);
for i = 1 : num_elem
    aper_set(i) = xdc_single_rect([elem_x(i),elem_y(i),elem_z(i)],width,height,num_sub_w,num_sub_h,focus);
end

aper_test = xdc_concatenate(aper_set, focus);

for i = 1 : num_elem
    xdc_free(aper_set(i));
end

data_ref = xdc_get(aper_ref);
data_test = xdc_get(aper_test);

xdc_display_aper(aper_ref);
title('reference aperture');
xdc_display_aper(aper_test);
title('test aperture');

xdc_free(aper_ref);
xdc_free(aper_test);
