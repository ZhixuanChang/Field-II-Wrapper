function [rf, t_start] = calc_scan_rf(fc, fs, tx_aper, rx_aper, tx_delay, tx_apod, scat_pos, scat_amp, c, scanning_path)
% CALC_SCAN_RF calculates the scattering waves with emission and reception apertures being placed at a group of
%   positions.
%
% INPUT
%   fc          - the center frequency of the emission and the reception aperture, [Hz].
%   fs          - the sampling rate of the simulation, [Hz].
%   tx_aper     - the created emission aperture.
%   rx_aper     - the created reception aperture.
%   tx_delay    - the emission delays for each physical elements in each emission cycle, in size of [nu x nv] where nu
%                 is the number of emission cycles and nv is the number of physical elements in the emission aperture.
%   tx_apod     - the emission apodization for each physical elements in each emission cycle. The size and meaning of 
%                 the elements are the same as tx_delay.
%   scat_pos    - the position of the scatterers in unit of [m]. The size is [ns x 3] where ns is the number of
%                 scatterers. The first, second, and last columns are the x-, y-, z-coordinates, respectively.
%   scat_amp    - the scattering amplitude of the scatterers. The size is [ns x 1].
%   c           - the sound of speed in the medium, [m/s].
%
%
% OUTPUT
%   rf          - the RF data in size of [nt, nv, nu, ns], where nt is the length of the temporal sequence, nv is the
%                 number of reception channels, nu is the number of emission cycles in a frame of FMC, ns is the number
%                 of scanning points specified by the scanning_path. The scattering waves should be aligned with same
%                 t_start.
%   t_start     - the reception start time.
%
% NOTES
%   - To avoid freeing and re-initializing apertures at each scanning position, this function equivalently moves
%     the scatterers in the opposite direction. If the aperture moves to position [x, y, z], the scatterers are
%     shifted by [-x, -y, -z].
%   - All RF data from different scanning positions are aligned to the same t_start.
%
% EXAMPLE
%   % Create scanning path
%   scanning_path = scan_zigzag(10, 10, 0.5e-3, 0.5e-3);
%
%   % Run simulation
%   [rf, t0] = calc_scan_rf(fc, fs, tx_aper, rx_aper, tx_delay, tx_apod, ...
%                           scat_pos, scat_amp, c, scanning_path);
%   % rf size: [nt x nv x nu x 100]
%
% See also
%   calc_rf, scan_zigzag

    %% Input validation
    if nargin < 10
        error('calc_scan_rf requires 10 input arguments');
    end

    % Check scanning_path dimensions
    if size(scanning_path, 2) ~= 3
        error('scanning_path must be [ns x 3] matrix');
    end

    num_scan_pos = size(scanning_path, 1);  % Number of scanning positions

    %% Loop through each scanning position
    fprintf('Starting scanning simulation with %d positions...\n', num_scan_pos);

    rf_cells = cell(num_scan_pos, 1);     % Store RF data for each scan position
    t_starts = zeros(num_scan_pos, 1);    % Store t_start for each scan position

    for i_scan = 1:num_scan_pos
        % Get current scanning position
        scan_pos = scanning_path(i_scan, :);  % [x, y, z]

        % Move scatterers in opposite direction (equivalent to moving apertures)
        % If aperture moves to [x, y, z], scatterers move by [-x, -y, -z]
        scat_pos_shifted = scat_pos - scan_pos;

        % Calculate RF data for this scanning position
        [rf_temp, t_start_temp] = calc_rf(fc, fs, tx_aper, rx_aper, ...
                                          tx_delay, tx_apod, ...
                                          scat_pos_shifted, scat_amp, c);

        % Store results
        rf_cells{i_scan} = rf_temp;
        t_starts(i_scan) = t_start_temp;

        % Display progress
        if mod(i_scan, max(1, floor(num_scan_pos/10))) == 0 || i_scan == num_scan_pos
            fprintf('Scanning progress: %d/%d (%.1f%%)\n', i_scan, num_scan_pos, 100*i_scan/num_scan_pos);
        end
    end

    %% Align all RF data to the same t_start
    % Find the minimum start time across all scanning positions
    t_start = min(t_starts);

    % Calculate delay in samples for each scanning position
    delay_samples = round((t_starts - t_start) * fs);

    % Find maximum length needed
    max_len = 0;
    for i_scan = 1:num_scan_pos
        current_len = delay_samples(i_scan) + size(rf_cells{i_scan}, 1);
        if current_len > max_len
            max_len = current_len;
        end
    end

    % Get dimensions from first RF data
    % rf_cells{1} has size [nt x nv x nu]
    nv = size(rf_cells{1}, 2);  % Number of receive channels
    nu = size(rf_cells{1}, 3);  % Number of emissions

    % Preallocate output matrix [nt x nv x nu x ns]
    rf = zeros(max_len, nv, nu, num_scan_pos);

    % Fill aligned data
    for i_scan = 1:num_scan_pos
        start_idx = delay_samples(i_scan) + 1;
        nt_i = size(rf_cells{i_scan}, 1);
        end_idx = start_idx + nt_i - 1;

        % rf_cells{i_scan} has size [nt_i x nv x nu]
        rf(start_idx:end_idx, :, :, i_scan) = rf_cells{i_scan};
    end

    %% Display summary
    fprintf('\nScanning simulation completed:\n');
    fprintf('  Number of scanning positions: %d\n', num_scan_pos);
    fprintf('  Number of emissions per position: %d\n', nu);
    fprintf('  Number of receive channels: %d\n', nv);
    fprintf('  RF data size: [%d x %d x %d x %d]\n', size(rf, 1), size(rf, 2), size(rf, 3), size(rf, 4));
    fprintf('  Global start time: %.6f us\n', t_start*1e6);
    fprintf('  Time span: %.6f us\n', max_len/fs*1e6);
end
