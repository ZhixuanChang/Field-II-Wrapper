function aper = xdc_matrix_array(varargin)
% Generate 2-D matrix array aperture.
%
% The matrix must be parallel to x-axis and y-axis, and the normal direction of the probe plane must be parallel to
% z-axis.
%
% Usage
% aper = xdc_matrix_array(nx, pitch_x, width_x, ny, pitch_y, width_y, n_sub_x, n_sub_y, focus, center);
%   nx - the number of elements along x-axis.
%   pitch_x - the center distance between adjacent elements along x-axis, unit: [m].
%   width_x - the width of the element along x-axis, unit: [m].
%   ny - the number of elements along y-axis.
%   pitch_y - the center distance between adjacent elements along y-axis, unit: [m].
%   width_y - the width of the element along y-axis, unit: [m].
%   n_sub_x - the subdivision number of each physical element along x-axis.
%   n_sub_y - the subdivision number of each physical element along y-axis.
%   focus - the fixed focus point of the aperture, should be given as a matrix with a size of [1, 3], unit: [m].
%   center - the center of the array, unit: [m].
%
%   focus and center are optional, the default value of focus is [0, 0, 1000], while the default value of center
%     is [0, 0, 0].
%
% See also
% xdc_single_rect, xdc_rectangles
%
% Author: Zhixuan Chang (chang_zhixuan@outlook.com)
% Date: January 06, 2026
% Version: 0.2
% Comment: Add argument `center` which allow the user to set the center of the array.
%
% Author: Zhixuan Chang (chang_zhixuan@outlook.com)
% Date: January 06, 2026
% Version: 0.1
% Comment: Create file.
% 
if nargin == 8
    % xdc_matrix_array(nx, pitch_x, width_x, ny, pitch_y, width_y, n_sub_x, n_sub_y);
    nx = varargin{1};
    pitch_x = varargin{2};
    width_x = varargin{3};
    ny = varargin{4};
    pitch_y = varargin{5};
    width_y = varargin{6};
    n_sub_x = varargin{7};
    n_sub_y = varargin{8};
    focus = [0, 0, 1000];
    center = [0, 0, 0];
elseif nargin == 10
    % xdc_matrix_array(nx, pitch_x, width_x, ny, pitch_y, width_y, n_sub_x, n_sub_y, focus);
    nx = varargin{1};
    pitch_x = varargin{2};
    width_x = varargin{3};
    ny = varargin{4};
    pitch_y = varargin{5};
    width_y = varargin{6};
    n_sub_x = varargin{7};
    n_sub_y = varargin{8};
    focus = varargin{9};
    center = [0, 0, 0];
elseif nargin == 11
    % xdc_matrix_array(nx, pitch_x, width_x, ny, pitch_y, width_y, n_sub_x, n_sub_y, focus, center);
    nx = varargin{1};
    pitch_x = varargin{2};
    width_x = varargin{3};
    ny = varargin{4};
    pitch_y = varargin{5};
    width_y = varargin{6};
    n_sub_x = varargin{7};
    n_sub_y = varargin{8};
    focus = varargin{9};
    center = varargin{10};
else
    error("Invalid arguments number. Use help xdc_matrix_array to get the supported calling statements.");
end

x_vec = (-(nx-1)/2 : (nx-1)/2) * pitch_x + center(1);
y_vec = (-(ny-1)/2 : (ny-1)/2) * pitch_y + center(2);

sub_width_x = width_x / n_sub_x;
sub_width_y = width_y / n_sub_y;

n_sub = n_sub_x * n_sub_y;

rect = zeros(nx * ny, 19);
for i = 1 : ny
    for j = 1 : nx
        ind = (i-1) * nx + j;
        sub_x_vec = (-(n_sub_x-1)/2 : (n_sub_x-1)/2) * sub_width_x + x_vec(j);
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
end

[x, y] = meshgrid(x_vec, y_vec);
x = reshape(x, numel(x), 1);
y = reshape(y, numel(y), 1);

x_center = [x, y, zeros(nx*ny,1)];

aper = xdc_rectangles(rect, x_center, focus);

end
