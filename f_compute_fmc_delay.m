function delay = f_compute_fmc_delay(roi_size, roi_pixel_size, roi_pos, elem_pos, c)

img_size = round(roi_size ./ roi_pixel_size);

x_vec = (-(img_size(1)-1)/2 : (img_size(1)-1)/2) * roi_pixel_size(1) + roi_pos(1);
y_vec = (-(img_size(2)-1)/2 : (img_size(2)-1)/2) * roi_pixel_size(2) + roi_pos(2);
[x, y] = meshgrid(x_vec, y_vec);

delay = zeros(img_size(2), img_size(1), size(elem_pos,1));

for i = 1 : size(elem_pos,1)
    delay(:,:,i) = sqrt((x - elem_pos(i,1)).^2 + (y - elem_pos(i,2)).^2) / c;
end