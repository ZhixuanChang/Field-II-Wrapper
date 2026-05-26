% This script is used to simulate ultrasonic waves

clear;
close all;

%% Global parameters
fc = 1e6;
fs = 100e6;
c = 1500;

nx = 16;
ny = 16;
pitch = 0.75e-3;
width = 0.7e-3;
kerf = pitch - width;
math_size = 0.1e-3;
focus = [0, 0, 1000];

scat_x_vec = -4e-3 : 4e-3 : 4e-3;
scat_y_vec = -4e-3 : 4e-3 : 4e-3;
scat_z_vec = 15e-3;

[scat_x, scat_y, scat_z] = ndgrid(scat_x_vec, scat_y_vec, scat_z_vec);
scat_pos = [reshape(scat_x, [], 1), reshape(scat_y, [], 1), reshape(scat_z, [], 1)];
scat_amp = ones(size(scat_pos,1), 1);

% figure(1);
% scatter3(scat_pos(:,1), scat_pos(:,2), scat_pos(:,3));

% return;

%% 

field_init();

num_sub_x = round(width / math_size);
num_sub_y = round(width / math_size);

tx_aper = xdc_linear_multirow(nx, width, ny, width * ones(ny, 1), kerf, kerf, num_sub_x, num_sub_y, focus);
rx_aper = xdc_linear_multirow(nx, width, ny, width * ones(ny, 1), kerf, kerf, num_sub_x, num_sub_y, focus);

impulse_response = sin(2 * pi * fc * (0 : 1/fs : 2/fc));
impulse_response = impulse_response .* hanning(length(impulse_response))';
excitation = sin(2 * pi * fc * (0 : 1/fs : 2/fc));
excitation = excitation .* hanning(length(excitation))';

num_elem = nx * ny;

xdc_impulse(tx_aper, impulse_response);
xdc_impulse(rx_aper, impulse_response);
xdc_excitation(tx_aper, impulse_response);
xdc_apodization(rx_aper, 0, ones(1, num_elem));
xdc_times_focus(tx_aper, 0, zeros(1, num_elem));
xdc_times_focus(rx_aper, 0, zeros(1, num_elem));

rf_data_set = cell(1, num_elem);
t_start_set = zeros(1, num_elem);

bar_length = 40;    % 进度条的视觉长度 (字符数)
reverse_str = '';   % 初始化退格字符串 (关键：用于清除上一行)

tic

fprintf("Simulation start ...\n");
for i = 1 : num_elem
    tx_apod = zeros(1, num_elem);
    tx_apod(i) = 1;
    xdc_apodization(tx_aper, 0, tx_apod);
    [rf_data_set{i}, t_start_set(i)] = calc_scat_multi(tx_aper, rx_aper, scat_pos, scat_amp);

    percent = i / num_elem;
    num_equals = round(percent * bar_length);
    if num_equals == 0
        bar_str = ['[', repmat(' ', 1, bar_length), ']'];
    elseif num_equals == bar_length
        bar_str = ['[', repmat('=', 1, bar_length), ']'];
    else
        bar_str = ['[', repmat('=', 1, num_equals-1), '>', ...
                   repmat(' ', 1, bar_length - num_equals), ']'];
    end
    msg = sprintf('%s %.1f%% (%d/%d)', bar_str, percent * 100, i, num_elem);
    fprintf([reverse_str, msg]);
    reverse_str = repmat(sprintf('\b'), 1, length(msg));
end
fprintf("\nSimulation finished!\n");

toc

field_end();

[rfdata, rx_delay] = f_rf_comb(rf_data_set, t_start_set, fs);

disp(size(rfdata));

% return

%% save data

% save configuration
config.fc = fc;
config.fs = fs;
config.bandwidth = 0.8;
elem_x_vec = (-(nx - 1) / 2 : (nx - 1) / 2) * pitch;
elem_y_vec = (-(ny - 1) / 2 : (ny - 1) / 2) * pitch;
[elem_y, elem_x] = meshgrid(elem_y_vec, elem_x_vec);
config.elemPos = [reshape(elem_x, [], 1), reshape(elem_y, [], 1), zeros(num_elem, 1)];
config.nx = nx;
config.ny = ny;
config.nt = size(rfdata,1);
config.c = c;
config.rxDelay = rx_delay;
config.roiPos.x = 0e-3;
config.roiPos.y = 0e-3;
config.roiPos.z = 15e-3;
config.roiPixelSize.x = 0.5e-3;
config.roiPixelSize.y = 0.5e-3;
config.roiPixelSize.z = 0.5e-3;
config.roiSize.x = 15e-3;
config.roiSize.y = 15e-3;
config.roiSize.z = 10e-3;

fid = fopen("matlab_ex01.json", "w");
fwrite(fid, jsonencode(config, 'PrettyPrint', true));
fclose(fid);

file_name = "matlab_ex01.h5";
h5create(file_name, "/rfdata", size(rfdata), "Datatype", "single");
h5write(file_name, "/rfdata", rfdata);
h5create(file_name, "/t_start", size(rx_delay), "Datatype", "single");
h5write(file_name, "/t_start", rx_delay);
h5create(file_name, "/fs", size(fs), "Datatype", "single");
h5write(file_name, "/fs", fs);
h5disp(file_name);
