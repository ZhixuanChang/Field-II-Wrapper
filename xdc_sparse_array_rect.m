function aper = xdc_sparse_array_rect(phy_elem_pos, phy_elem_wx, phy_elem_wy, math_elem_size, focus)
% generate sparse array.
%
% Parameters
% phy_elem_pos - the center position of each physical element. Each row for an element. The first, second, and third
%                column are x, y, and z coordinates, respectively.
% phy_elem_wx - the size of the physical element along x-axis, [m].
% phy_elem_wy - the size of the physical element along y-axis, [m].
% math_elem_size - the size of a mathematical element, [m].
% focus - the fixed focus point of the array.
%
% Example
% aper = xdc_sparse_array_rect(phy_elem_pos, phy_elem_wx, phy_elem_wy, math_elem_size, focus);
%
% See also
% xdc_concatenate, xdc_rectangles, xdc_single_rect

num_elem = size(phy_elem_pos, 1);
n_sub_x = ceil(phy_elem_wx / math_elem_size);
n_sub_y = ceil(phy_elem_wy / math_elem_size);

aper_set = zeros(1, num_elem);
for i = 1 : num_elem
    aper_set(i) = xdc_single_rect(phy_elem_pos(i, :), phy_elem_wx, phy_elem_wy, n_sub_x, n_sub_y, focus);
end

aper = xdc_concatenate(aper_set, focus);

for i = 1 : num_elem
    xdc_free(aper_set(i));
end

end
