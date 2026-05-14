clear;
close all;

%% Global Parameters
fs_sim          = 100e6;            % Sampling frequency [Hz]
c               = 1500;             % Speed of sound [m/s]
f0              = 1.0e6;            % Transducer center frequency [Hz]
pitch           = 0.75e-3;          % Pitch [m]
width           = 0.7e-3;           % Width of element
height          = 5e-3;             % Height of element [m]
num_elem        = 128;              % Number of elements in the transducer
num_sub_w       = 5;                % number of sub-elements along width direction
num_sub_h       = 10;               % number of sub-elements along height direction
focus           = [0 0 1000];       % Fixed focal point [m]
pulse_width     = 2;                % The excitation pulse width in count of period
div             = 10;               % Sampling frequency division number

aper            = pitch * num_elem; % The physical aperture size of the array
scat_x_range    = [-0.75 * aper, 0.75 * aper]; % The x range of the scatterers
scat_z_range    = [0.2 * aper, 1.5 * aper];    % The z range of the scatterers

num_samples     = 200;              % The number of samples in a single dataset file

%% Execute Simulation
field_init();

set_field('fs', fs_sim);
set_field('c',  c);

impulse_response = sin(2 * pi * f0 * (0 : 1/fs_sim : pulse_width/f0));
impulse_response = impulse_response .* hanning(length(impulse_response))';

excitation = sin(2 * pi * f0 * (0 : 1/fs_sim : pulse_width/f0));

tx_aper = xdc_linear_array(num_elem, width, height, pitch-width, num_sub_w, num_sub_h, focus);
rx_aper = xdc_linear_array(num_elem, width, height, pitch-width, num_sub_w, num_sub_h, focus);

xdc_impulse(tx_aper, impulse_response);
xdc_excitation(tx_aper, excitation);
xdc_times_focus(tx_aper, 0, zeros(1,num_elem));

xdc_impulse(rx_aper, impulse_response);
xdc_apodization(rx_aper, 0, ones(1,num_elem));

fs = fs_sim / div;

dataset = cell(1, num_samples);
tstartset = zeros(1,num_samples);
nt_max = 0;

fprintf('Simulating %d samples...\n', num_samples);
tic_total = tic;
for j = 1 : num_samples
    elapsed = toc(tic_total);
    avg_time = elapsed / j;
    eta = avg_time * (num_samples - j);
    fprintf('\r  %d/%d (%.1f%%) | elapsed: %.1fs | avg: %.2fs/iter | ETA: %.1fs ', ...
        j, num_samples, j/num_samples*100, elapsed, avg_time, eta);

    num_scat = randi([4, 12]);
    scat_x_vec = rand(num_scat, 1) * diff(scat_x_range) + scat_x_range(1);
    scat_z_vec = rand(num_scat, 1) * diff(scat_z_range) + scat_z_range(1);
    
    scatter_pos = [scat_x_vec, zeros(num_scat, 1), scat_z_vec];
    scatter_amp = ones(num_scat, 1);
    
    rf_data_set = cell(1, num_elem);
    t_start_set = zeros(1, num_elem);
    
    for i = 1 : num_elem
        tx_apod = zeros(1,num_elem);
        tx_apod(i) = 1;
        xdc_apodization(tx_aper, 0, tx_apod);
        [v, t1] = calc_scat_multi(tx_aper, rx_aper, scatter_pos, scatter_amp);
        rf_data_set{i} = v;
        t_start_set(i) = t1;
    end
    
    [rf_data_sim, t_start] = f_rf_comb(rf_data_set, t_start_set, fs_sim);
    rf_data_sim = rf_data_sim / max(abs(rf_data_sim(:)));
    
    
    dataset{j} = single(rf_data_sim(1:div:end,:,:));
    nt_max = max(nt_max, size(dataset{j},1));
    tstartset(j) = t_start;
end
fprintf('\nSimulation complete. Total time: %.1fs\n', toc(tic_total));

field_end();

% rf_data_tmp = zeros(nt_max, num_elem, num_elem, 100, "single");

rf_data_set2 = zeros(nt_max, num_elem, num_elem, num_samples, "single");
for i = 1 : num_samples
    rf_data_set2(1 : size(dataset{1,i},1), :, :, i) = dataset{1,i};
end

rf_data_set2 = single(rf_data_set2);
tstartset = single(tstartset);

file_name = "sparse_fmc_dataset_002.h5";
h5create(file_name, "/rf_data", size(rf_data_set2), 'Datatype','single');
h5write(file_name, "/rf_data", rf_data_set2);
h5create(file_name, "/t_start", size(tstartset),'Datatype','single');
h5write(file_name, "/t_start", tstartset);
h5disp(file_name);

return;

%% Postprocessing

div = 5;
fs_div = fs_sim / div;
rf_data_div = rf_data_sim(1:div:end,:,:);
nt = size(rf_data_div, 1);

roi_width = diff(scat_x_range);
roi_depth = diff(scat_z_range);
roi_pos = [mean(scat_x_range), mean(scat_z_range)];
l0 = c / f0;
dx = 0.5 * l0;
dz = 0.5 * l0;
nx = round(roi_width / dx);
nz = round(roi_depth / dz);
x_vec = (-(nx-1)/2 : (nx-1)/2) * dx + roi_pos(1);
z_vec = (-(nz-1)/2 : (nz-1)/2) * dz + roi_pos(2);
[x, z] = meshgrid(x_vec, z_vec);

elem_x = (-(num_elem-1)/2 : (num_elem-1)/2) * pitch;

rx_tof = zeros(nz, nx, num_elem);
f = waitbar(0, "Calculating focal delay");
for i = 1 : num_elem
    waitbar(i/num_elem, f, "Calculating focal delay");
    rx_tof(:, :, i) = sqrt((x - elem_x(i)).^2 + z.^2) / c;
end
close(f);

img = zeros(nz, nx);
iq_data = hilbert(rf_data_div);
t_array = (0:nt-1) / fs_div;
f = waitbar(0, "Reconstructing image");
for i = 1 : num_elem
    waitbar(i/num_elem, f, "Reconstructing image");
    for j = 1 : num_elem
        img = img + interp1(t_array, iq_data(:, j, i), rx_tof(:, :, i) + rx_tof(:, :, j) - t_start, 'linear', 0);
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
