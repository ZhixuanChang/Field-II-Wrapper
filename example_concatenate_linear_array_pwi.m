% Verify the xdc_concatenate() method by generating a linear array and implementing PWI.
%
% It is assumed that field_init has been called before.
%
% Author: Zhixuan Chang
% Date: Nov. 6, 2025
% Version: 0.1

clear;
close all;
clc;

%% generate linear array aperture
% generate each element aperture
num_elem = 64;
pitch = 0.5e-3;
width = 0.45e-3;
height = 5e-3;
focus = [0, 0, 1000];
fc = 5e6;

elem_x = (-(num_elem-1)/2 : (num_elem-1)/2) * pitch;
elem_y = zeros(1, num_elem);
elem_z = zeros(1, num_elem);

% dx = 0.1e-3;
% n_sub_x = ceil(width / dx);
% n_sub_y = ceil(height / dx);
n_sub_x = 5;
n_sub_y = 10;

aper_set = zeros(1, num_elem);

for i = 1 : num_elem
    aper_set(i) = xdc_single_rect([elem_x(i), elem_y(i), elem_z(i)], width, height, n_sub_x, n_sub_y, focus);
end

tx_aperture = xdc_concatenate(aper_set, focus);
rx_aperture = xdc_concatenate(aper_set, focus);

for i = 1 : num_elem
    xdc_free(aper_set(i));
end

% xdc_display_aper(tx_aperture);

%% set global parameters
fs = 100e6;
c = 1540;

set_field('fs', fs);
set_field('c', c);

%% set aperture impulse response and excitation signal
impulse_response = sin(2*pi*fc*(0:1/fs:2/fc));
impulse_response = impulse_response .* hanning(length(impulse_response))';
excitation = sin(2*pi*fc*(0:1/fs:2/fc));

xdc_impulse(tx_aperture, impulse_response);
xdc_impulse(rx_aperture, impulse_response);
xdc_excitation(tx_aperture, excitation);
xdc_apodization(tx_aperture, 0, ones(1, num_elem));
xdc_apodization(rx_aperture, 0, ones(1, num_elem));

%% calculate scattering signals
num_emit = 21;
steered_angles = linspace(-10, 10, 21) / 180 * pi;
tx_delay = sin(steered_angles') * elem_x / c;
tx_delay = tx_delay - min(tx_delay, [], 'all');

scat_pos = [
    0, 0, 20e-3;
    0, 0, 30e-3;
    0, 0, 40e-3;
    0, 0, 50e-3;
    -10e-3, 0, 40e-3;
    10e-3, 0, 40e-3;
];
scat_amp = ones(size(scat_pos, 1), 1);

rf_data_set = cell(1, num_emit);
t_start_set = zeros(1, num_emit);

% xdc_times_focus(rx_aperture, 0, zeros(1, num_elem));
for i = 1 : num_emit
    xdc_times_focus(tx_aperture, 0, tx_delay(i, :));
    [rf_cur, t_start] = calc_scat_multi(tx_aperture, rx_aperture, scat_pos, scat_amp);
    rf_data_set{i} = rf_cur;
    t_start_set(i) = t_start;
end

xdc_free(tx_aperture);
xdc_free(rx_aperture);

rfdata = f_rf_comb(rf_data_set, t_start_set, fs);

%% reconstructing images
roi_width = 30e-3;
roi_depth = 40e-3;
roi_pos = [0, 35e-3];
roi_dx = 0.1e-3;
roi_dz = 0.1e-3;
roi_nx = round(roi_width/roi_dx);
roi_nz = round(roi_depth/roi_dz);
roi_x_vec = (-(roi_nx-1)/2 : (roi_nx-1)/2) * roi_dx + roi_pos(1);
roi_z_vec = (-(roi_nz-1)/2 : (roi_nz-1)/2) * roi_dz + roi_pos(2);

[roi_x, roi_z] = meshgrid(roi_x_vec, roi_z_vec);

tx_tof = zeros(roi_nz, roi_nx, num_emit);
rx_tof = zeros(roi_nz, roi_nx, num_elem);

f = waitbar(0, 'Calculating incident ToF');
for i = 1 : num_emit
    waitbar(i/num_emit, f, 'Calculating incident ToF');
    tx_tof(:, :, i) = ((roi_x - elem_x(1)) * sin(steered_angles(i)) + (roi_z - elem_z(1)) * cos(steered_angles(i))) / c + tx_delay(i, 1);
end
close(f);

f = waitbar(0, 'Calculating scattering ToF');
for i = 1 : num_elem
    waitbar(i/num_elem, f, 'Calculating scattering ToF');
    rx_tof(:, :, i) = sqrt((roi_x - elem_x(i)).^2 + (roi_z - elem_z(i)).^2) / c;
end
close(f);

iqdata = hilbert(rfdata);
img = zeros(roi_nz, roi_nx);
t_array = (0 : size(rfdata,1)-1) / fs;

tic
f = waitbar(0, 'Reconstructing image');
for i = 1 : num_emit
    waitbar(i/num_emit, f, 'Reconstructing image');
    for j = 1 : num_elem
        img = img + interp1(t_array, iqdata(:, j, i), tx_tof(:, :, i) + rx_tof(:, :, j), "linear", 0);
    end
end
close(f);
toc

img = abs(img);

figure();
imagesc(roi_x_vec*1e3, roi_z_vec*1e3, img);
axis equal tight;
xlabel('x (mm)');
ylabel('z (mm)');
