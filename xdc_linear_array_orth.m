function aper = xdc_linear_array_orth(varargin)
% Generate linear array aperture with elements arranged along y-axis.
%
% The linear array is arranged along y-axis, and the normal direction of the probe plane must be parallel to
% z-axis. Each element is a rectangle.
%
% Usage
% aper = xdc_linear_array_orth(n_elements, width, height, pitch, n_sub_y, n_sub_x, focus, center);
%   n_elements - the number of elements along y-axis.
%   width - the width of the element along y-axis (pitch direction), unit: [m].
%   height - the height of the element along x-axis (elevation direction), unit: [m].
%   pitch - the center distance between adjacent elements along y-axis, unit: [m].
%   n_sub_y - the subdivision number of each physical element along y-axis (width direction).
%   n_sub_x - the subdivision number of each physical element along x-axis (height direction).
%   focus - the fixed focus point of the aperture, should be given as a matrix with a size of [1, 3], unit: [m].
%   center - the center of the array, unit: [m].
%
%   focus and center are optional, the default value of focus is [0, 0, 1000], while the default value of center
%     is [0, 0, 0].
%
% See also
% xdc_single_rect, xdc_rectangles, xdc_matrix_array
%
% Author: Zhixuan Chang (chang_zhixuan@outlook.com)
% Date: June 04, 2026
% Version: 0.1
% Comment: Create file.

if nargin == 6
    % xdc_linear_array_orth(n_elements, width, height, pitch, n_sub_y, n_sub_x);
    n_elements = varargin{1};
    width = varargin{2};
    height = varargin{3};
    pitch = varargin{4};
    n_sub_y = varargin{5};
    n_sub_x = varargin{6};
    focus = [0, 0, 1000];
    center = [0, 0, 0];
elseif nargin == 7
    % xdc_linear_array_orth(n_elements, width, height, pitch, n_sub_y, n_sub_x, focus);
    n_elements = varargin{1};
    width = varargin{2};
    height = varargin{3};
    pitch = varargin{4};
    n_sub_y = varargin{5};
    n_sub_x = varargin{6};
    focus = varargin{7};
    center = [0, 0, 0];
elseif nargin == 8
    % xdc_linear_array_orth(n_elements, width, height, pitch, n_sub_y, n_sub_x, focus, center);
    n_elements = varargin{1};
    width = varargin{2};
    height = varargin{3};
    pitch = varargin{4};
    n_sub_y = varargin{5};
    n_sub_x = varargin{6};
    focus = varargin{7};
    center = varargin{8};
else
    error("Invalid arguments number. Use help xdc_linear_array_orth to get the supported calling statements.");
end

% Only one element along x-axis, elements are arranged along y-axis
nx = 1;
ny = n_elements;

x_vec = center(1);
y_vec = (-(ny-1)/2 : (ny-1)/2) * pitch + center(2);

sub_width_x = height / n_sub_x;
sub_width_y = width / n_sub_y;

n_sub = n_sub_x * n_sub_y;

rect = zeros(nx * ny, 19);
for i = 1 : ny
    ind = i;
    sub_x_vec = (-(n_sub_x-1)/2 : (n_sub_x-1)/2) * sub_width_x + x_vec;
    sub_y_vec = (-(n_sub_y-1)/2 : (n_sub_y-1)/2) * sub_width_y + y_vec(i);
    [sub_x, sub_y] = meshgrid(sub_x_vec, sub_y_vec);
    sub_x = reshape(sub_x, numel(sub_x), 1);
    sub_y = reshape(sub_y, numel(sub_y), 1);
    rect((ind-1)*n_sub_x*n_sub_y+1 : ind*n_sub_x*n_sub_y, :) = [ind * ones(n_sub,1), ...
        sub_x - sub_width_x/2, sub_y - sub_width_y/2, zeros(n_sub,1) + center(3), ...
        sub_x + sub_width_x/2, sub_y - sub_width_y/2, zeros(n_sub,1) + center(3), ...
        sub_x + sub_width_x/2, sub_y + sub_width_y/2, zeros(n_sub,1) + center(3), ...
        sub_x - sub_width_x/2, sub_y + sub_width_y/2, zeros(n_sub,1) + center(3), ...
        ones(n_sub,1), sub_width_x*ones(n_sub,1), sub_width_y*ones(n_sub,1), ...
        sub_x, sub_y, zeros(n_sub,1) ...
        ];
end

x = x_vec * ones(ny, 1);
y = y_vec(:);

x_center = [x, y, zeros(ny, 1)];

aper = xdc_rectangles(rect, x_center, focus);

end
