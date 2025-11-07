% example of sparse array imaging.

clear;
close all;
clc;

%%
elem_wx = 2e-3;
elem_wy = 2e-3;
pitch_x = 2e-3;
pitch_y = 2e-3;

% 三维仿真，从一个32x32的面阵换能器中随机选择256个阵元
nx_full = 32;
ny_full = 32;
elem_x_full_vec = (-(nx_full-1)/2 : (nx_full-1)/2) * pitch_x;
elem_y_full_vec = (-(ny_full-1)/2 : (ny_full-1)/2) * pitch_y;

[elem_x_full, elem_y_full] = meshgrid(elem_x_full_vec, elem_y_full_vec);
elem_x_full_seq = reshape(elem_x_full, numel(elem_x_full), 1);
elem_y_full_seq = reshape(elem_y_full, numel(elem_y_full), 1);

random_order = randperm(nx_full * ny_full);
num_sparse = 256;
elem_ind = random_order(1:num_sparse);
elem_ind = sort(elem_ind);
elem_x_seq = elem_x_full_seq(elem_ind);
elem_y_seq = elem_y_full_seq(elem_ind);

elem_pos = [elem_x_seq, elem_y_seq, zeros(num_sparse, 1)];

math_size = 0.5e-3;
focus = [0, 0, 1000];

tx_aperture = xdc_sparse_array_rect(elem_pos, elem_wx, elem_wy, math_size, focus);
rx_aperture = xdc_sparse_array_rect(elem_pos, elem_wx, elem_wy, math_size, focus);

% xdc_display_aper(tx_aperture);

%
fs = 100e6;
c = 1500;

set_field('fs', fs);
set_field('c', c);

%
fc = 1e6;

impulse_response = sin(2*pi*fc*(0 : 1/fs : 2/fc));
impulse_response = impulse_response .* hanning(numel(impulse_response))';

excitation = sin(2*pi*fc*(0 : 1/fs : 2/fc));

xdc_impulse(tx_aperture, impulse_response);
xdc_impulse(rx_aperture, impulse_response);

xdc_excitation(tx_aperture, excitation);

xdc_apodization(rx_aperture, 0, ones(1, num_sparse));

%
scat_x_vec = -25e-3 : 25e-3 : 25e-3;
scat_y_vec = -25e-3 : 25e-3 : 25e-3;
scat_z_vec = 50e-3;
[scat_x, scat_y, scat_z] = meshgrid(scat_x_vec, scat_y_vec, scat_z_vec);
scat_x = reshape(scat_x, numel(scat_x), 1);
scat_y = reshape(scat_y, numel(scat_y), 1);
scat_z = reshape(scat_z, numel(scat_z), 1);
scat_pos = [scat_x, scat_y, scat_z];
scat_amp = ones(size(scat_pos, 1), 1);

figure()
scatter3(scat_x, scat_y, scat_z);

f_div = 4;

rf_data_set = cell(1, num_sparse);
t_start_set = zeros(1, num_sparse);
f = waitbar(0, 'Calculating RF data');
for i = 1 : num_sparse
    waitbar(i/num_sparse, f, 'Calculating RF data');
    tx_apod = zeros(1, num_sparse);
    tx_apod(i) = 1;
    xdc_apodization(tx_aperture, 0, tx_apod);
    [rf_cur, t_start] = calc_scat_multi(tx_aperture, rx_aperture, scat_pos, scat_amp);
    rf_data_set{i} = rf_cur(1:f_div:end, :);
    t_start_set(i) = t_start;
end
close(f);

fs = fs / f_div;

[rfdata, t_start] = f_rf_comb(rf_data_set, t_start_set, fs);
rfdata = single(rfdata);

%
xdc_free(tx_aperture);
xdc_free(rx_aperture);

%%
roi_size = [64e-3, 64e-3, 20e-3];
roi_pos = [0, 0, 50e-3];
roi_dx = [1e-3, 1e-3, 1e-3];

roi_nx = round(roi_size ./ roi_dx);
roi_x_vec = (-(roi_nx(1)-1)/2 : (roi_nx(1)-1)/2) * roi_dx(1) + roi_pos(1);
roi_y_vec = (-(roi_nx(2)-1)/2 : (roi_nx(2)-1)/2) * roi_dx(2) + roi_pos(2);
roi_z_vec = (-(roi_nx(3)-1)/2 : (roi_nx(3)-1)/2) * roi_dx(3) + roi_pos(3);

[x, y, z] = meshgrid(roi_x_vec, roi_y_vec, roi_z_vec);
tx_tof = zeros(roi_nx(2), roi_nx(1), roi_nx(3), num_sparse);

x = single(x);
y = single(y);
z = single(z);
elem_pos = single(elem_pos);

f = waitbar(0, 'Calculating focal delay');
for i = 1 : num_sparse
    waitbar(i/num_sparse, f, 'Calculating focal delay');
    tx_tof(:, :, :, i) = sqrt((x-elem_pos(i,1)).^2 + (y-elem_pos(i,2)).^2 + (z-elem_pos(i,3)).^2) / c;
end
close(f);


num_samp = size(rfdata, 1);
t_array = (0:(num_samp-1)) / fs + t_start;

img = single(zeros(roi_nx(2), roi_nx(1), roi_nx(3), num_sparse));

tic
iqdata = hilbert(rfdata);
% f = waitbar(0, 'Reconstructing image');
parfor i = 1 : num_sparse
    % waitbar(i/num_sparse, f, 'Reconstructing image');
    for j = 1 : num_sparse
        img = img + interp1(t_array, iqdata(:,j,i), tx_tof(:,:,:,i)+tx_tof(:,:,:,j), 'linear', 0);
    end
end
% close(f);
img = sum(img, 4);
img = abs(img);
toc

%%
% figure(101);
% [~, z_ind] = min(abs(roi_z_vec - 50e-3));
% imagesc(roi_x_vec, roi_y_vec, img(:, :, z_ind));
% axis equal tight;
% colorbar;
