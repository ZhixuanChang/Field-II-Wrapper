% Example of use of the new Field II program running under Matlab
%
% This example shows how a linear array B-mode system scans an image
%
% This script assumes that the field_init procedure has been called
%
% Example by Joergen Arendt Jensen, Version 2.0, March 22, 2011.
clear
close all
clc

%% Set global parameters
% Generate the transducer apertures for send and receive
fs              = 100e6;            % Sampling frequency [Hz]
c               = 1540;             % Speed of sound [m/s]

% sampling frequency
set_field('fs', fs);
set_field('c', c);

% Set scatterers
scatter_pos = [
    0, 0, 30e-3;
    0, 0, 40e-3;
    0, 0, 50e-3;
    0, 0, 60e-3;
    -10e-3, 0, 50e-3;
    10e-3,  0, 50e-3;
    ];
scatter_amp = ones(size(scatter_pos, 1), 1);

figure();
scatter3(scatter_pos(:, 1), scatter_pos(:, 2), scatter_pos(:, 3), 'o');

%% set aperture
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

% generate impulse response
impulse_response = sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response = impulse_response .* hanning(max(size(impulse_response)))';

% generate excitation signal
excitation=sin(2*pi*f0*(0:1/fs:2/f0));

% Generate aperture for emission
% emit_aperture = xdc_linear_array_arb(num_elem, pitch, width, height, num_sub_w, num_sub_h, center, angle, focus);
emit_aperture = xdc_matrix_array(num_elem,pitch,width,1,height,height,num_sub_w,num_sub_h,focus);
xdc_impulse(emit_aperture, impulse_response);
xdc_excitation(emit_aperture, excitation);
xdc_apodization(emit_aperture, 0, ones(1, num_elem));

% Generate aperture for reception
% receive_aperture = xdc_linear_array_arb(num_elem, pitch, width, height, num_sub_w, num_sub_h, center, angle, focus);
receive_aperture = xdc_matrix_array(num_elem,pitch,width,1,height,height,num_sub_w,num_sub_h,focus);
xdc_impulse(receive_aperture, impulse_response);
xdc_apodization(receive_aperture, 0, ones(1, num_elem));

%% Calculate scattering echoes

% calculate emission delays
elem_x = (-(num_elem-1)/2 : (num_elem-1)/2) * pitch;
tx_angle = linspace(-10, 10, 21) / 180 * pi;
tx_delay = sin(tx_angle') * elem_x / c;
tx_delay = tx_delay - min(tx_delay, [], "all");

% Pre-allocate some storage
rf_data_set = cell(1, num_elem);
t_start_set = zeros(1, num_elem);

f = waitbar(0, "Calculating ...");
for i = 1 : length(tx_angle)
    waitbar(i / num_elem, f, "Calculating ...");
    xdc_times_focus(emit_aperture, 0, tx_delay(i, :));
    [v, t1] = calc_scat_multi(emit_aperture, receive_aperture, scatter_pos, scatter_amp);
    rf_data_set{i} = v;
    t_start_set(i) = t1;
end
close(f);

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

% num_samp = round(times * fs);
% for i = 1 : num_elem
%     num_samp(i) = num_samp(i) + size(rf_data_set{i}, 1);
% end
% num_samp = max(num_samp, [], "all");
% 
% rf_data = zeros(num_samp, num_elem, num_elem);
% for i = 1 : num_elem
%     t_start = round(times(i)*fs);
%     rf_data(t_start : t_start + size(rf_data_set{i}, 1) - 1, :, i) = rf_data_set{i};
% end

[rf_data, t_start] = f_rf_comb(rf_data_set, t_start_set, fs);

%% TFM imaging
roi_width = 40e-3;
roi_depth = 40e-3;
roi_pos = [0, 45e-3];
dx = 0.1e-3;
dz = 0.1e-3;
nx = round(roi_width / dx);
nz = round(roi_depth / dz);
x_vec = (-(nx-1)/2 : (nx-1)/2) * dx + roi_pos(1);
z_vec = (-(nz-1)/2 : (nz-1)/2) * dz + roi_pos(2);
[x, z] = meshgrid(x_vec, z_vec);

tx_tof = zeros(nz, nx, num_elem);
for i = 1 : length(tx_angle)
    tx_tof(:, :, i) = ((x - elem_x(1)) * sin(tx_angle(i)) + z * cos(tx_angle(i))) / c + tx_delay(i, 1);
end
rx_tof = zeros(nz, nx, num_elem);
f = waitbar(0, "Calculating focal delay");
for i = 1 : num_elem
    waitbar(i/num_elem, f, "Calculating focal delay");
    rx_tof(:, :, i) = sqrt((x - elem_x(i)).^2 + z.^2) / c - t_start;
end
close(f);

%%
img = zeros(nz, nx);
iq_data = hilbert(rf_data);
t_array = (0:size(rf_data, 1)-1) / fs;
f = waitbar(0, "Reconstructing image");
for i = 1 : length(tx_angle)
    waitbar(i/num_elem, f, "Reconstructing image");
    for j = 1 : num_elem
        img = img + interp1(t_array, iq_data(:, j, i), tx_tof(:, :, i) + rx_tof(:, :, j), 'linear', 0);
    end
end
close(f);
img = abs(img);

figure();
imagesc(x_vec*1e3, z_vec*1e3, img);
axis equal tight;
hold on;
plot(scatter_pos(:,1)*1e3, scatter_pos(:,3)*1e3, 'ro');
hold off;
xlabel('x (mm)');
ylabel('z (mm)');
