%%
% Example Name: example_cross_array
% Author: Z. Chang
% Date: January 06, 2026
% Comment: example for underwater imaging with cross array and multiple beams imaging.
%%
clear;

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

%% set aperture
f0              = 2.25e6;              % Transducer center frequency [Hz]
pitch           = 0.4e-3;           % Pitch [m]
width           = 0.35e-3;          % Width of element
height          = 5e-3;             % Height of element [m]
num_elem        = 64;               % Number of elements in the transducer
num_sub_w       = 5;                % number of sub-elements along width direction
num_sub_h       = 10;               % number of sub-elements along height direction
center          = [0, 0, 0];        % center of the array
angle           = 0;                % angle between the array and the x-axis, [rad]
focus           = [0, 0, 1000];

% generate impulse response
impulse_response = sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response = impulse_response .* hanning(max(size(impulse_response)))';

% generate excitation signal
excitation=sin(2*pi*f0*(0:1/fs:2/f0));

% Generate aperture for emission
center_h = [0, 0, 0];
h_aper = xdc_matrix_array(num_elem,pitch,width,1,height,height,num_sub_w,num_sub_h,focus,center_h);
xdc_impulse(h_aper, impulse_response);
xdc_excitation(h_aper, excitation);
xdc_apodization(h_aper, 0, ones(1, num_elem));

% Generate aperture for reception
center_v = [0, num_elem*pitch/2+3e-3, 0];
v_aper = xdc_matrix_array(1,height,height,num_elem,pitch,width,num_sub_h,num_sub_w,focus,center_v);
xdc_impulse(v_aper, impulse_response);
xdc_apodization(v_aper, 0, ones(1, num_elem));

xdc_display_aper(h_aper);
xdc_display_aper(v_aper);

return;

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
    xdc_times_focus(h_aper, 0, tx_delay(i, :));
    [v, t1] = calc_scat_multi(h_aper, v_aper, scatter_pos, scatter_amp);
    rf_data_set{i} = v;
    t_start_set(i) = t1;
end
close(f);

% Free space for apertures
xdc_free(h_aper)
xdc_free(v_aper)

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
