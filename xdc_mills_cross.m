function [horiz_aper, vert_aper] = xdc_mills_cross(num_elem, pitch, width, height, n_sub_x, n_sub_y, focus)
% Generate Mills cross array with two perpendicular linear arrays.
%
% The Mills cross array consists of a horizontal array parallel to the x-axis and a vertical array parallel to the
% y-axis. Both arrays share the same element width, height, and pitch parameters, but are oriented perpendicular to
% each other. The horizontal array elements have width in x-direction and height in y-direction, while the vertical
% array elements have height in x-direction and width in y-direction (rotated 90 degrees).
%
% Mills cross arrays are typically used for multi-beam imaging. The limited element height creates larger beam
% divergence compared to RCA arrays.
%
% Parameters
% num_elem - the number of elements in each array (both arrays have the same number of elements).
% pitch - the center distance between adjacent elements, unit: [m]. Same for both arrays.
% width - the element width, unit: [m]. For horizontal array, this is the size along x-axis.
%         For vertical array, this is the size along y-axis.
% height - the element height, unit: [m]. For horizontal array, this is the size along y-axis.
%          For vertical array, this is the size along x-axis.
% n_sub_x - the subdivision number along x-axis for horizontal array elements (along y-axis for vertical array).
% n_sub_y - the subdivision number along y-axis for horizontal array elements (along x-axis for vertical array).
% focus - the fixed focus point of both apertures, should be given as [x, y, z], unit: [m].
%
% Return
% horiz_aper - the horizontal array aperture (parallel to x-axis).
% vert_aper - the vertical array aperture (parallel to y-axis).
%
% Example
% [horiz_aper, vert_aper] = xdc_mills_cross(64, 0.3e-3, 0.27e-3, 5e-3, 3, 10, [0, 0, 1000]);
%
% See also
% xdc_matrix_array, xdc_single_rect, xdc_concatenate
%
% Author: Zhixuan Chang (https://zhixuanchang.github.io/)
% Date: Feb. 8, 2026
% Version: 0.1

if nargin < 7
    error('xdc_mills_cross requires 7 input arguments.');
end

if numel(focus) ~= 3
    error('focus must be a vector with 3 elements [x, y, z].');
end

% Generate horizontal array (parallel to x-axis)
% Elements arranged along x-axis, centered at y=0, z=0
% Element dimensions: width (in x) × height (in y)
% Using xdc_matrix_array: (nx, pitch_x, width_x, ny, pitch_y, width_y, n_sub_x, n_sub_y, focus, center)
horiz_aper = xdc_matrix_array(num_elem, pitch, width, 1, height, height, n_sub_x, n_sub_y, focus, [0, 0, 0]);

% Generate vertical array (parallel to y-axis)
% Elements arranged along y-axis, centered at x=0, z=0
% Element dimensions: height (in x) × width (in y) - rotated 90 degrees
% For vertical array: nx=1 (single column), ny=num_elem (elements along y)
vert_aper = xdc_matrix_array(1, height, height, num_elem, pitch, width, n_sub_y, n_sub_x, focus, [0, 0, 0]);

end
