function img = f_fmc_das(rf_data, t_array, delay)

if size(rf_data,3) ~= size(rf_data,2)
    error("The input rf_data is not a valid FMC dataset.");
end
if ~isvector(t_array)
    error("The input t_array is not a valid temporal label.");
end
if size(rf_data,1) ~= length(t_array)
    error("The size of rf_data does not match t_array.");
end
if size(rf_data,2) ~= size(delay,3)
    error("The size of rf_data does not match delay.");
end

img = zeros(size(delay,1), size(delay,2));

iq_data = hilbert(rf_data);
for i = 1 : size(rf_data,3)
    for j = 1 : size(rf_data,2)
        img = img + interp1(t_array, iq_data(:,j,i), delay(:,:,i) + delay(:,:,j), "linear", 0.0);
    end
end
img = abs(img);
