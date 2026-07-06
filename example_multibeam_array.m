clear;
close all;

%% Global configurations
c = 1500; % speed of sound in sea
fc = 1e6; % center frequency
fs = 100e6; % sampling rate

l0 = c / fc;
pitch = l0 / 2; % array pitch
kerf = 0.05e-3;
width = pitch - kerf;
height = 10e-3;
num_elem = 128;
focus = [0, 0, 1000];

aper = num_elem * pitch;
fprintf("Aperture size: %f\n", aper);
% return;

% scatter distribution
% dx = 1e-3;
% h = 5;
% mean_depth = 3.5;
% seabed_terrain = generate_seabed_terrain(dx, h, mean_depth);

% figure(1);
% plot(seabed_terrain(:, 1), seabed_terrain(:, 2));
% axis tight;
% return;

% scat_pos = [seabed_terrain(:,1), zeros(size(seabed_terrain,1),1), seabed_terrain(:,2)];
% scat_amp = ones(size(seabed_terrain,1),1);

scat_x_vec = -300e-3 : 100e-3 : 300e-3;
scat_y_vec = 0;
scat_z_vec = 200e-3 : 100e-3 : 400e-3;

[scat_x_mesh, scat_y_mesh, scat_z_mesh] = meshgrid(scat_x_vec, scat_y_vec, scat_z_vec);
scat_pos = [scat_x_mesh(:), scat_y_mesh(:), scat_z_mesh(:)];
scat_amp = ones(size(scat_pos,1), 1);

% beamforming parameters
steered_angles = (-60 : 0.2 : 60) / 180 * pi;

fprintf("Scatter number: %d\n", length(scat_amp));
fprintf("Emission cycles: %d\n", length(steered_angles));

% return

%% Execute simulation
field_init;

set_field('c', c);
set_field('fs', fs);

math_size = l0 / 10;
num_sub_w = round(width / math_size);
num_sub_h = round(height / math_size);

tx_aper = xdc_linear_array(num_elem, width, height, kerf, num_sub_w,num_sub_h, focus);
rx_aper = xdc_linear_array(num_elem, width, height, kerf, num_sub_w,num_sub_h, focus);

impulse_response = sin(2*pi*fc*(0:1/fs:2/fc));
impulse_response = impulse_response .* hanning(max(size(impulse_response)))';

excitation = sin(2*pi*fc*(0:1/fs:2/fc));

xdc_impulse(tx_aper, impulse_response);
xdc_impulse(rx_aper, impulse_response);
xdc_excitation(tx_aper, excitation);
xdc_apodization(tx_aper, 0, ones(1, num_elem));
xdc_apodization(rx_aper, 0, ones(1, num_elem));

elem_x = (-(num_elem-1)/2 : (num_elem-1)/2) * pitch;
tx_delay = sin(steered_angles') * elem_x ./ c;
tx_delay = tx_delay - min(tx_delay(:));

rf_data_set = cell(1, length(steered_angles));
t_start_set = zeros(1, num_elem);

bar_length = 40;    % 进度条的视觉长度 (字符数)
reverse_str = '';   % 初始化退格字符串 (关键：用于清除上一行)

fprintf("Simulation start ...\n");

for i = 1 : length(steered_angles)
    xdc_times_focus(tx_aper, 0, tx_delay(i, :));
    [v, t1] = calc_scat_multi(tx_aper, rx_aper, scat_pos, scat_amp);
    rf_data_set{i} = v;
    t_start_set(i) = t1;
    
    percent = i / length(steered_angles);
    num_equals = round(percent * bar_length);
    if num_equals == 0
        bar_str = ['[', repmat(' ', 1, bar_length), ']'];
    elseif num_equals == bar_length
        bar_str = ['[', repmat('=', 1, bar_length), ']'];
    else
        bar_str = ['[', repmat('=', 1, num_equals-1), '>', ...
                   repmat(' ', 1, bar_length - num_equals), ']'];
    end
    msg = sprintf('%s %.1f%% (%d/%d)', bar_str, percent * 100, i, length(steered_angles));
    fprintf([reverse_str, msg]);
    reverse_str = repmat(sprintf('\b'), 1, length(msg));
end

fprintf("Simulation finished!\n");

xdc_free(tx_aper);
xdc_free(rx_aper);

field_end;

%% Saving data

[rfdataRaw, t_start] = f_rf_comb(rf_data_set, t_start_set, fs);
disp(size(rfdataRaw));

rfdata = rfdataRaw(1:5:end, :, :);
fs_new = fs / 5;
save("multibeam_sim01.mat", "fs_new", "c", "fc", "steered_angles", "t_start");

h5create("multibeam_sim01_rfdata.h5", "/rfdata", size(rfdata), "Datatype", "single");
h5write("multibeam_sim01_rfdata.h5", "/rfdata", single(rfdata));
