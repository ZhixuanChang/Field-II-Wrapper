% This script is used to simulate ultrasonic waves

clear;
close all;

%% Global parameters
fc = 1e6;
fs = 100e6;
c = 1500;

nx = 16;
ny = 16;
pitch = 10.0e-3;
width = 8.5e-3;
kerf = pitch - width;
math_size = 0.25e-3;
focus = [0, 0, 1000];

scat_x_vec = -50e-3 : 50e-3 : 50e-3;
scat_y_vec = -50e-3 : 50e-3 : 50e-3;
scat_z_vec = 120e-3;

[scat_x, scat_y, scat_z] = ndgrid(scat_x_vec, scat_y_vec, scat_z_vec);
scat_center_pos = [reshape(scat_x, [], 1), reshape(scat_y, [], 1), reshape(scat_z, [], 1)];

scat_pos = [];
for i = 1 : size(scat_center_pos, 1)
    scat_tmp = gen_sphere_solid(scat_center_pos(i, :), 1e-3, 0.5e-3);
    scat_pos = [scat_pos; scat_tmp];
end
scat_amp = ones(size(scat_pos,1), 1);

disp("The number of scatterers:");
disp(size(scat_amp,1));

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
rfdata = rfdata ./ max(abs(rfdata(:)));

disp(size(rfdata));

div = 10;
rfdata = rfdata(1:div:end, :, :);
fs = fs / div;

% return

%% save configuration

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
config.roiPos.z = 120e-3;
config.roiPixelSize.x = 1e-3;
config.roiPixelSize.y = 1e-3;
config.roiPixelSize.z = 1e-3;
config.roiSize.x = 180e-3;
config.roiSize.y = 160e-3;
config.roiSize.z = 40e-3;

fid = fopen("data/prj_uw_target_pose_identify.json", "w");
fwrite(fid, jsonencode(config, 'PrettyPrint', true));
fclose(fid);

%% save RF data

file_name = "data/prj_uw_target_pose_identify.h5";
if exist(file_name, 'file')
    delete(file_name);
end
h5create(file_name, "/rfdata", size(rfdata), "Datatype", "single");
h5write(file_name, "/rfdata", single(rfdata));
h5create(file_name, "/t_start", size(rx_delay), "Datatype", "single");
h5write(file_name, "/t_start", single(rx_delay));
h5create(file_name, "/fs", size(fs), "Datatype", "single");
h5write(file_name, "/fs", single(fs));
h5disp(file_name);

return;

%% Postprocessing

close all;

% Volume viewing
img = h5read("data/prj_uw_target_pose_identify_img.h5", "/img");
volumeViewer(img);

% MIP curves along each axis
img_nx = round(config.roiSize.x / config.roiPixelSize.x);
img_ny = round(config.roiSize.y / config.roiPixelSize.y);
img_nz = round(config.roiSize.z / config.roiPixelSize.z);

img_dx = config.roiPixelSize.x;
img_dy = config.roiPixelSize.y;
img_dz = config.roiPixelSize.z;

img_x_vec = (-img_nx/2 : (img_nx/2 - 1)) * img_dx + config.roiPos.x;
img_y_vec = (-img_ny/2 : (img_ny/2 - 1)) * img_dy + config.roiPos.y;
img_z_vec = (-img_nz/2 : (img_nz/2 - 1)) * img_dz + config.roiPos.z;

[~, x_slice_ind] = min(abs(img_x_vec - 0));
img_y_slice = squeeze(img(x_slice_ind));
img_y_slice = img_y_slice';

figure();
tiledlayout(2, 1, "TileSpacing", "compact", "Padding", "compact");
nexttile();
imagesc(img_y_vec * 1e3, img_z_vec * 1e3, img_y_slice);
xlabel("y (mm)");
ylabel("z (mm)");
nexttile();
plot(img_y_vec, max(img_y_slice));
xlabel("y (mm)");
ylabel("Amplitude");
% f_show_image(img(x_slice_ind, :, :));
