function [rfdata, t_start] = f_rf_comb(rf_data_set, t_start_set, fs)
% Align RF data of multiple emission cycles.
%
% The calculated RF signals of Field is unaligned. To obtain the spatial spectrum, the signal must be aligned first.
%
% Parameters
% rf_data_set - a cell contains all scattering signals.
% t_start_set - the start time vector.
% fs - the sampling rate.
%
% Return
% rfdata - the aligned RF data.
% t_start - the unified reception delay for all channels.
%
% Example
% rfdata = h_rf_comb(rf_data_set, t_start_set, fs);
%
% See also
% calc_scat_all, calc_scat_multi

num_emit = length(t_start_set);
num_elem = size(rf_data_set{1}, 2);

t_start = min(t_start_set);
t_start_set = t_start_set - t_start;

num_padd = round(t_start_set * fs);
num_samp = zeros(1, num_emit);
for i = 1 : num_emit
    num_samp(i) = num_samp(i) + size(rf_data_set{i}, 1);
end

num_samp = num_samp + num_padd;
num_samp = max(num_samp, [], 'all');

rfdata = zeros(num_samp, num_elem, num_emit);
for i = 1 : num_emit
    rfdata(num_padd(i)+1 : num_padd(i)+size(rf_data_set{i},1), :, i) = rf_data_set{i};
end

end
