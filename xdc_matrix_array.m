function aper = xdc_matrix_array(nx, pitch_x, width_x, ny, pitch_y, width_y, n_sub_x, n_sub_y, focus)
% Generate 2-D matrix array aperture.
%
% The matrix must be parallel to x-axis and y-axis, and the normal direction of the probe plane must be parallel to
% z-axis.
%
% Parameters
% nx - the number of elements along x-axis.
% pitch_x - the center distance between adjacent elements along x-axis.
% width_x - the width of the element along x-axis.
% ny - the number of elements along y-axis.
% pitch_y - the center distance between adjacent elements along y-axis.
% width_y - the width of the element along y-axis.
% n_sub_x - the subdivision number of each physical element along x-axis.
% n_sub_y - the subdivision number of each physical element along y-axis.
% focus - the fixed focus point of the aperture, should be in the shape of [x, y, z].
%
% Example
% aper = xdc_matrix_array(nx, pitch_x, width_x, ny, pitch_y, width_y, n_sub_x, n_sub_y, focus);
%
% See also
% xdc_single_rect, xdc_rectangles

x_vec = (-(nx-1)/2 : (nx-1)/2) * pitch_x;
y_vec = (-(ny-1)/2 : (ny-1)/2) * pitch_y;

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
            sub_x - sub_width_x/2, sub_y - sub_width_y/2, zeros(n_sub,1), ...
            sub_x + sub_width_x/2, sub_y - sub_width_y/2, zeros(n_sub,1), ...
            sub_x + sub_width_x/2, sub_y + sub_width_y/2, zeros(n_sub,1), ...
            sub_x - sub_width_x/2, sub_y + sub_width_y/2, zeros(n_sub,1), ...
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
