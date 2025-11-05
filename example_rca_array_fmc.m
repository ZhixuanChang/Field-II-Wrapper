% Example of use of the new Field II program running under MATLAB.
%
% This example shows how a RCA transducer acquire RF data in FMC mode. The row elements are arranged along the y-axis,
% while the column elements are arranged along x-axis. The RF data is an array with a shape of Nt x Nv x Nu, where Nt
% is the sampling length, Nv is the number of reception elements, Nu is the number of emission cycles. For a RCA array
% with Nr row elements and Nc column elements, Nv = max(Nr, Nc), Nu = Nr + Nc.
%
% This script assumes that the field_init procedure has been called.
%
% Example by Joergen Arendt Jensen, Version 2.0, March 22, 2011.
% Revised by Zhixuan Chang, Version 2.1, November 5, 2025.

clear
close all
clc

%% Parameters input
% parameters listed in this block need to be specified by user
fs              = 100e6;                         % Sampling frequency [Hz]
c               = 1540;                          % Speed of sound [m/s]

f0              = 1e6;                           % Transducer center frequency [Hz]
pitch           = 0.5e-3;                        % Pitch [m]
width           = 0.45e-3;                       % Width of element
num_elem        = 64;                            % Number of elements in the transducer

center          = [0, 0, 0];                     % center of the array in world coordinates
focus           = [0, 0, 1000];                  % Fixed focal point [m], set focus(3) to a large value to approximate a flat transducer

% scatter_pos = [                                  % scatterers position, each row for a scatterer, the three columns are x, y, z-coordinate, respectively, unit: [m]
%     -10e-3, -10e-3, 50e-3;
%     -10e-3, 0e-3,   50e-3;
%     -10e-3, 10e-3,  50e-3;
%     0e-3,   -10e-3, 50e-3;
%     0e-3,   0e-3,   50e-3;
%     0e-3,   10e-3,  50e-3;
%     10e-3,  -10e-3, 50e-3;
%     10e-3,  0e-3,   50e-3;
%     10e-3,  10e-3,  50e-3;
%     ];
% scatter_amp = ones(size(scatter_pos, 1), 1);     % scattering amplitude, each row for a scatterer, the number of rows must match that of scatter_pos
scatter_pos = [10e-3, 10e-3, 50e-3];
scatter_amp = [1];

%% Set global parameters
set_field('fs', fs);
set_field('c', c);

%% Set aperture
delta_x         = c/f0/10;                        % maximum size of a single patch of an element
height          = num_elem*pitch-(pitch-width);  % Height of element [m]
num_sub_w       = ceil(width/delta_x);           % number of sub-elements along width direction
num_sub_h       = ceil(height/delta_x);          % number of sub-elements along height direction
% num_sub_w       = 1;           % number of sub-elements along width direction
% num_sub_h       = 5;          % number of sub-elements along height direction

impulse_response = sin(2*pi*f0*(0:1/fs:1/f0));
impulse_response = impulse_response .* hanning(max(size(impulse_response)))';
excitation = sin(2*pi*f0*(0:1/fs:1/f0));
% excitation = excitation .* hanning(max(size(excitation)))';

%% Set scatterers

% figure();
% scatter3(scatter_pos(:, 1), scatter_pos(:, 2), scatter_pos(:, 3), 'o');
% title("Scatterer distribution");
% xlabel x
% ylabel y
% zlabel z

%% Calculate scattering echoes
% Pre-allocate some storage
rf_data_set = cell(1, num_elem * 2);
times = zeros(1, num_elem * 2);

% row emission, column reception
% emit_aperture = xdc_linear_array_arb(num_elem, pitch, width, height, num_sub_w, num_sub_h, center, pi/2, focus);
emit_aperture = xdc_matrix_array(1,height,height,num_elem,pitch,width,num_sub_h,num_sub_w,focus);
xdc_impulse(emit_aperture, impulse_response);
xdc_excitation(emit_aperture, excitation);
% xdc_display_aper(emit_aperture);
% title("Emit aperture");

% receive_aperture = xdc_linear_array_arb(num_elem, pitch, width, height, num_sub_w, num_sub_h, center, 0, focus);
receive_aperture = xdc_matrix_array(num_elem,pitch,width,1,height,height,num_sub_w,num_sub_h,focus);
xdc_impulse(receive_aperture, impulse_response);
xdc_apodization(receive_aperture, 0, ones(1, num_elem));
% xdc_display_aper(receive_aperture);
% title("Receive aperture");
% return

f = waitbar(0, "Calculating Row Emission Column Reception ...");
for i = 1 : num_elem
    try
        waitbar(i / num_elem, f, "Calculating Row Emission Column Reception ...");
        tx_apod = zeros(1, num_elem);
        tx_apod(i) = 1;
        xdc_apodization(emit_aperture, 0, tx_apod);
        [rf_cur, t_start] = calc_scat_multi(emit_aperture, receive_aperture, scatter_pos, scatter_amp);
        rf_data_set{i} = rf_cur;
        times(i) = t_start;
    catch e
        fprintf("[ERROR] RTCR simulation is broken, Message: %s\n", e.message);
    end
end
close(f);

xdc_free(emit_aperture);
xdc_free(receive_aperture);


% column emission, row reception
% emit_aperture = xdc_linear_array_arb(num_elem, pitch, width, height, num_sub_w, num_sub_h, center, 0, focus);
emit_aperture = xdc_matrix_array(num_elem,pitch,width,1,height,height,num_sub_w,num_sub_h,focus);
xdc_impulse(emit_aperture, impulse_response);
xdc_excitation(emit_aperture, excitation);

% receive_aperture = xdc_linear_array_arb(num_elem, pitch, width, height, num_sub_w, num_sub_h, center, pi/2, focus);
receive_aperture = xdc_matrix_array(1,height,height,num_elem,pitch,width,num_sub_h,num_sub_w,focus);
xdc_impulse(receive_aperture, impulse_response);
xdc_apodization(receive_aperture, 0, ones(1, num_elem));

f = waitbar(0, 'Calculating Column Emission Row Reception ...');
for i = 1 : num_elem
    try
        waitbar(i / num_elem, f, 'Calculating Column Emission Row Reception ...');
        tx_apod = zeros(1, num_elem);
        tx_apod(i) = 1;
        xdc_apodization(emit_aperture, 0, tx_apod);
        [rf_cur, t_start] = calc_scat_multi(emit_aperture, receive_aperture, scatter_pos, scatter_amp);
        rf_data_set{i+num_elem} = rf_cur;
        times(i+num_elem) = t_start;
    catch e
        fprintf("[ERROR] CTRR simulation is broken, Message: %s\n", e.message);
    end
end
close(f);

xdc_free(emit_aperture);
xdc_free(receive_aperture);

% return
%% Check A Scan ToF
% ind = 1;
err_max = 0;
for ind = 1 : num_elem

rf_env = abs(hilbert(rf_data_set{ind}));
% figure();
% tiledlayout(2, 2, "TileSpacing", "compact", "Padding", "compact");
% nexttile();
% imagesc(rf_env);

[~, tof] = max(rf_env, [], 1);
tof = tof / fs + times(ind);

tx_elem_x = zeros(1, num_elem);
tx_elem_y = (-(num_elem-1)/2 : (num_elem-1)/2) * pitch;
tx_elem_z = zeros(1, num_elem);

rx_elem_x = (-(num_elem-1)/2 : (num_elem-1)/2) * pitch;
rx_elem_y = zeros(1, num_elem);
rx_elem_z = zeros(1, num_elem);

ref = (sqrt((tx_elem_y(ind) - scatter_pos(1,2)).^2 + (tx_elem_z(ind) - scatter_pos(1,3)).^2) ...
    + sqrt((rx_elem_x - scatter_pos(1,1)).^2 + (rx_elem_z - scatter_pos(1,3)).^2)) / c;
% nexttile();
% % figure();
% plot(tof);
% hold on;
% plot(ref);
% hold off;
% legend("ToF", "Reference", 'Location', 'southwest');

err_vec = abs(tof - ref);
err_max_cur = max(err_vec, [], 'all');
err_max = max(err_max_cur, err_max);
% fprintf("[INFO] The maximum error is %e us\n", 1e6 * max(err, [], "all"));


rf_env = abs(hilbert(rf_data_set{ind+num_elem}));
% % figure();
% nexttile();
% imagesc(rf_env);

[~, tof] = max(rf_env, [], 1);
tof = tof / fs + times(ind + num_elem);

tx_elem_x = (-(num_elem-1)/2 : (num_elem-1)/2) * pitch;
tx_elem_y = zeros(1, num_elem);
tx_elem_z = zeros(1, num_elem);
rx_elem_x = zeros(1, num_elem);
rx_elem_y = (-(num_elem-1)/2 : (num_elem-1)/2) * pitch;
rx_elem_z = zeros(1, num_elem);
ref = (sqrt((tx_elem_x(ind) - scatter_pos(1,1)).^2 + (tx_elem_z(ind) - scatter_pos(1,3)).^2) ...
    + sqrt((rx_elem_y - scatter_pos(1,2)).^2 + (rx_elem_z - scatter_pos(1,3)).^2)) / c;

% % figure();
% nexttile();
% plot(tof);
% hold on;
% plot(ref);
% hold off;
% legend("ToF", "Reference", 'Location', 'southwest');

err_vec = abs(tof - ref);
err_max_cur = max(err_vec, [], 'all');
err_max = max(err_max_cur, err_max);
% fprintf("[INFO] The maximum error is %e us\n", 1e6 * max(err, [], "all"));

end

fprintf("[INFO] The maximum error is %e us\n", 1e6 * err_max);

% return

num_samp = round(times * fs);
for i = 1 : length(times)
    num_samp(i) = num_samp(i) + size(rf_data_set{i}, 1);
end
num_samp = max(num_samp, [], "all");

rf_data = zeros(num_samp, num_elem, num_elem*2);
for i = 1 : length(times)
    t_start = round(times(i)*fs);
    rf_data(t_start : t_start + size(rf_data_set{i}, 1) - 1, :, i) = rf_data_set{i};
end
