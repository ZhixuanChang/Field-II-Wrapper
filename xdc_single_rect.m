function aper = xdc_single_rect(center, len_x, len_y, n_sub_x, n_sub_y, focus)
% Generate an aperture contains a single rectangular physical element.
%
% The aperture must be parallel to x-axis and y-axis, and its normal direction must be parallel to z-axis.
% The center of the element can be in arbitrary position.
%
% Parameters
% center - the center of the element, [x,y,z].
% len_x - the length of the element along x-axis, [m].
% len_y - the length of the element along y-axis, [m].
% n_sub_x - subdivision number along x-axis.
% n_sub_y - subdivision number along y-axis.
% focus - the fixed focus point of the aperture.
%
% Example
% aper = xdc_single_rect(center, len_x, len_y, n_sub_x, n_sub_y, focus);
% 
% See also
% xdc_rectangles

n_sub = n_sub_x * n_sub_y;

sub_len_x = len_x / n_sub_x;
sub_len_y = len_y / n_sub_y;

sub_x_vec = (-(n_sub_x-1)/2 : (n_sub_x-1)/2) * sub_len_x + center(1);
sub_y_vec = (-(n_sub_y-1)/2 : (n_sub_y-1)/2) * sub_len_y + center(2);
[sub_x, sub_y] = meshgrid(sub_x_vec, sub_y_vec);
sub_x_seq = reshape(sub_x, numel(sub_x), 1);
sub_y_seq = reshape(sub_y, numel(sub_y), 1);

z_seq = center(3)*ones(n_sub,1);

rects = [ones(n_sub,1), ...
    sub_x_seq-sub_len_x/2, sub_y_seq-sub_len_y/2, z_seq, ...
    sub_x_seq+sub_len_x/2, sub_y_seq-sub_len_y/2, z_seq, ...
    sub_x_seq+sub_len_x/2, sub_y_seq+sub_len_y/2, z_seq, ...
    sub_x_seq-sub_len_x/2, sub_y_seq+sub_len_y/2, z_seq, ...
    ones(n_sub,1), sub_len_x*ones(n_sub,1), sub_len_y*ones(n_sub,1), ...
    sub_x_seq, sub_y_seq, z_seq];

aper = xdc_rectangles(rects,center,focus);

end
