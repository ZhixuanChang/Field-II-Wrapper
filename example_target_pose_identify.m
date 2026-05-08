% Use a group of acoustic reflectors to identify the pose of the target.
%
% 2026-04-14
% Zhixuan Chang

clear;
close all;

%% Global parameters
fc = 1e6;
fs = 100e6;
c = 1500;

% acoustic reflectors position in body frame
% Details:
%   - Column 1~3: the center position of each reflector pair.
%   - Column 4~6: the direction of each reflector pair.
%   - Column 7: the distance between the two reflectors in each pair.

% reflector_pose_body_frame = [
%      0.0e-3,   0.0e-3,   0.0e-3,          1,          0,          0,   8.0e-3;
%     70.0e-3,   0.0e-3,   0.0e-3,          0,          1,          0,  11.0e-3;
%      0.0e-3,  60.0e-3,   0.0e-3,          0,          0,          1,  14.5e-3;
%     25.0e-3,  25.0e-3,  50.0e-3,  1/sqrt(2),  1/sqrt(2),          0,  18.5e-3;
%     50.0e-3,  40.0e-3,  30.0e-3,          0,  1/sqrt(2),  1/sqrt(2),  23.0e-3;
% ];

% scat_pos = zeros(10, 3);

% for i = 1 : size(reflector_pose_body_frame, 1)
%     scat_pos((i-1)*2+1, 1) = reflector_pose_body_frame(i,1) ...
%         - 0.5 * reflector_pose_body_frame(i,4) * reflector_pose_body_frame(i,7);
%     scat_pos((i-1)*2+1, 2) = reflector_pose_body_frame(i,2) ...
%         - 0.5 * reflector_pose_body_frame(i,5) * reflector_pose_body_frame(i,7);
%     scat_pos((i-1)*2+1, 3) = reflector_pose_body_frame(i,3) ...
%         - 0.5 * reflector_pose_body_frame(i,6) * reflector_pose_body_frame(i,7);
%     scat_pos((i-1)*2+2, 1) = reflector_pose_body_frame(i,1) ...
%         + 0.5 * reflector_pose_body_frame(i,4) * reflector_pose_body_frame(i,7);
%     scat_pos((i-1)*2+2, 2) = reflector_pose_body_frame(i,2) ...
%         + 0.5 * reflector_pose_body_frame(i,5) * reflector_pose_body_frame(i,7);
%     scat_pos((i-1)*2+2, 3) = reflector_pose_body_frame(i,3) ...
%         + 0.5 * reflector_pose_body_frame(i,6) * reflector_pose_body_frame(i,7);
% end

scat_x_vec = -10e-3 : 5e-3 : 10e-3;
scat_y_vec = -10e-3 : 5e-3 : 10e-3;
scat_z_vec = 15e-3;

[scat_x, scat_y, scat_z] = ndgrid(scat_x_vec, scat_y_vec, scat_z_vec);
scat_pos = [reshape(scat_x, [], 1), reshape(scat_y, [], 1), reshape(scat_z, [], 1)];

% figure(1);
% scatter3(scat_pos(:,1), scat_pos(:,2), scat_pos(:,3));
% axis tight;

scat_amp = ones(size(scat_pos,1), 1);

pitch_x = 1.5e-3;
pitch_y = 1.5e-3;
num_elem_x = 16;
num_elem_y = 16;
kerf_x = 0.1e-3;
kerf_y = 0.1e-3;
math_size = 0.2e-3;
focus = [0, 0, 1000];

%% Execute simulation

field_init();

width = pitch_x - kerf_x;
height = pitch_y - kerf_y;
num_sub_x = round(width / math_size);
num_sub_y = round(height / math_size);
height_vec = ones(num_elem_y, 1) * height;
tx_aper = xdc_linear_multirow(num_elem_x, width, num_elem_y, height_vec, kerf_x, kerf_y, num_sub_x, num_sub_y, focus);
rx_aper = xdc_linear_multirow(num_elem_x, width, num_elem_y, height_vec, kerf_x, kerf_y, num_sub_x, num_sub_y, focus);

impulse_response = sin(2*pi*fc*(0:1/fs:2/fc));
impulse_response = impulse_response .* hanning(length(impulse_response))';
excitation = sin(2*pi*fc*(0:1/fs:2/fc));
excitation = excitation .* hanning(length(impulse_response))';

num_elem = num_elem_x * num_elem_y;

xdc_impulse(tx_aper, impulse_response);
xdc_impulse(rx_aper, impulse_response);
xdc_excitation(tx_aper, impulse_response);
xdc_apodization(rx_aper, 0, ones(1, num_elem));
xdc_times_focus(tx_aper, 0, zeros(1, num_elem));
xdc_times_focus(rx_aper, 0, zeros(1, num_elem));

bar_length = 40;    % 进度条的视觉长度 (字符数)
reverse_str = '';   % 初始化退格字符串 (关键：用于清除上一行)

rf_data_set = cell(1, num_elem);
t_start_set = zeros(1, num_elem);

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

field_end();

%% Save data
fs_div = 4;

[rfdata, rx_delay] = f_rf_comb(rf_data_set, t_start_set, fs);
rfdata = rfdata(1:4:end, :, :);
fs_new = fs / fs_div;

rfdata = rfdata / max(abs(rfdata(:)));
rfdata = round(rfdata * 32767); % 定点化

fid = fopen("target_pose_identify_sim02.bin", "wb");
fwrite(fid, rfdata, 'int16');
fclose(fid);

config.fc = fc;
config.fs = fs_new;
config.bandwidth = 0.8;
elem_x_vec = (-(num_elem_x - 1) / 2 : (num_elem_x - 1) / 2) * pitch_x;
elem_y_vec = (-(num_elem_y - 1) / 2 : (num_elem_y - 1) / 2) * pitch_y;
[elem_y, elem_x] = meshgrid(elem_y_vec, elem_x_vec);
config.elemPos = [reshape(elem_x, [], 1), reshape(elem_y, [], 1), zeros(num_elem, 1)];
config.nx = num_elem_x;
config.ny = num_elem_y;
config.nt = size(rfdata,1);
config.c = c;
config.rxDelay = rx_delay;
% config.roiPos.x = mean(reflector_pose_body_frame(:,1));
% config.roiPos.y = mean(reflector_pose_body_frame(:,2));
% config.roiPos.z = mean(reflector_pose_body_frame(:,3));
% config.roiPixelSize.x = 1e-3;
% config.roiPixelSize.y = 1e-3;
% config.roiPixelSize.z = 1e-3;
% config.roiSize.x = 160e-3;
% config.roiSize.y = 160e-3;
% config.roiSize.z = max(abs(reflector_pose_body_frame(:,3)-config.roiPos(3)))*2 + 20e-3;
config.roiPos.x = 0e-3;
config.roiPos.y = 0e-3;
config.roiPos.z = 15e-3;
config.roiPixelSize.x = 0.5e-3;
config.roiPixelSize.y = 0.5e-3;
config.roiPixelSize.z = 0.5e-3;
config.roiSize.x = 25e-3;
config.roiSize.y = 25e-3;
config.roiSize.z = 10e-3;

fid = fopen("target_pose_identify_sim02.json", "w");
fwrite(fid, jsonencode(config, 'PrettyPrint', true));
fclose(fid);

%% 