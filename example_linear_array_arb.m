% Example of use of the new Field II program running under Matlab
%
% This example shows how a linear array B-mode system scans an image
%
% This script assumes that the field_init procedure has been called
%
% Example by Joergen Arendt Jensen, Version 2.0, March 22, 2011.
clear
clc

% Generate the transducer apertures for send and receive
f0              = 5e6;              % Transducer center frequency [Hz]
fs              = 100e6;            % Sampling frequency [Hz]
c               = 5900;             % Speed of sound [m/s]
lambda          = c/f0;             % Wave length [m]
pitch           = 0.6e-3;           % Pitch [m]
width           = 0.5e-3;           % Width of element
height          = 5e-3;             % Height of element [m]
num_elem        = 64;               % Number of elements in the transducer
num_sub_w       = 5;                % number of sub-elements along width direction
num_sub_h       = 10;               % number of sub-elements along height direction
center          = [0, 0, 0];        % center of the array
angle           = 0;                % angle between the array and the x-axis, [rad]
focus           = [0 0 1000];       % Fixed focal point [m]

% Set the sampling frequency
set_sampling(fs);
set_field('c', c);
set_field('fs', fs);

% Generate aperture for emission
emit_aperture = xdc_linear_array_arb(num_elem, pitch, width, height, num_sub_w, num_sub_h, center, angle, focus);
% emit_aperture = xdc_linear_array(num_elem, width, height, pitch-width, 5, 50, focus);

% Set the impulse response and excitation of the emit aperture
impulse_response = sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response = impulse_response .* hanning(max(size(impulse_response)))';
xdc_impulse(emit_aperture, impulse_response);
excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (emit_aperture, excitation);

% Generate aperture for reception
receive_aperture = xdc_linear_array_arb(num_elem, pitch, width, height, num_sub_w, num_sub_h, center, angle, focus);
% receive_aperture = xdc_linear_array(num_elem, width, height, pitch-width, 5, 50, focus);

% Set the impulse response for the receive aperture
xdc_impulse(receive_aperture, impulse_response);

% Load the computer phantom
scatter_pos = zeros(400, 3);
scatter_pos(:, 1) = 0;
scatter_pos(:, 2) = (-(400-1)/2:(400-1)/2) * 0.1e-3;
scatter_pos(:, 3) = 50e-3;
% scatter_pos = [0, 0, 50e-3];
scatter_amp = ones(400, 1);

figure();
scatter3(scatter_pos(:, 1), scatter_pos(:, 2), scatter_pos(:, 3), 'o');

%%

% Do FMC

% Pre-allocate some storage
rf_data_set = cell(1, num_elem);
times = zeros(1, num_elem);

for i = 1 : num_elem
    tx_apo = zeros(1, num_elem);
    tx_apo(i) = 1;
    xdc_apodization(emit_aperture, 0, tx_apo);
    rx_apo = ones(1, num_elem);
    xdc_apodization(receive_aperture, 0, rx_apo);
    [v, t1] = calc_scat_multi(emit_aperture, receive_aperture, scatter_pos, scatter_amp);
    rf_data_set{i} = v;
    times(i) = t1;
end

% Free space for apertures
xdc_free(emit_aperture)
xdc_free(receive_aperture)

rf_env = abs(hilbert(v));
[~, tof] = max(rf_env, [], 1);
plot(tof)
tof = tof / fs + t1;
plot(tof)
elem_x = (-(num_elem-1)/2 : (num_elem-1)/2) * pitch;
ref = (sqrt((elem_x(1)-scatter_pos(1,1)).^2 + scatter_pos(1,3).^2) + sqrt((elem_x-scatter_pos(1,1)).^2 + scatter_pos(1,3).^2)) / c;
figure(1);
plot(tof);
hold on;
plot(ref);
hold off;
err = tof - ref;
fprintf("The maximum error is %e us\n", 1e6 * max(err, [], "all"));

num_samp = round(times * fs);
for i = 1 : num_elem
    num_samp(i) = num_samp(i) + size(rf_data_set{i}, 1);
end
num_samp = max(num_samp, [], "all");

rf_data = zeros(num_samp, num_elem, num_elem);
for i = 1 : num_elem
    t_start = round(times(i)*fs);
    rf_data(t_start : t_start + size(rf_data_set{i}, 1) - 1, :, i) = rf_data_set{i};
end

%% TFM imaging
roi_width = 50e-3;
roi_depth = 40e-3;
roi_pos = [0, 50e-3];
dx = 0.1e-3;
dz = 0.1e-3;
nx = round(roi_width / dx);
nz = round(roi_depth / dz);
x_vec = (-(nx-1)/2 : (nx-1)/2) * dx + roi_pos(1);
z_vec = (-(nz-1)/2 : (nz-1)/2) * dz + roi_pos(2);
[x, z] = meshgrid(x_vec, z_vec);

focal_tof = zeros(nz, nx, num_elem);
for i = 1 : num_elem
    focal_tof(:, :, i) = sqrt((x - elem_x(i)).^2 + z.^2) / c;
end

%%
img = zeros(nz, nx);
iq_data = hilbert(rf_data);
t_array = (0:size(rf_data, 1)-1) / fs;
for i = 1 : num_elem
    for j = 1 : num_elem
        img = img + interp1(t_array, iq_data(:, j, i), focal_tof(:, :, i) + focal_tof(:, :, j), 'linear', 0);
    end
end
img = abs(img);

figure();
imagesc(x_vec, z_vec, img);
axis equal tight;
hold on;
plot(scatter_pos(1,1), scatter_pos(1,3), 'ro');
hold off;
