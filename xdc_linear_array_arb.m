function [aper] = xdc_linear_array_arb(num_elem, pitch, width, height, num_sub_w, num_sub_h, center, angle, focus)
% Generate aperture for Field II linear array simulation.
%
% Parameters
% num_elem - the number of elements in the linear array.
% pitch - the center distance between the adjacent elements.
% width - the width of an element, also the size of the element along x-axis.
% height - the height of an element, also the size of the element along y-axis.
% num_sub_w - the divided number of sub-elements along width direction.
% num_sub_h - the divided number of sub-elements along height direction.
% center - the center of the array, [x, y, z]. The value of z must be 0.
% angle - the angle between the array line and the x-axis.
% focus - the fixed focus of the aperture.
%
% Examples
% aper = xdc_linear_array_arb(num_elem, pitch, width, height, center, angle);
%
% See also
% xdc_rectangles

num_elem_rects = num_sub_w * num_sub_h;

elem_x = (-(num_elem-1)/2 : (num_elem-1)/2) * pitch * cos(angle) + center(1);
elem_y = (-(num_elem-1)/2 : (num_elem-1)/2) * pitch * sin(angle) + center(2);
elem_z = zeros(1, num_elem);

sub_width = width / num_sub_w;
sub_height = height / num_sub_h;

rect = zeros(num_elem_rects * num_elem, 19);
for i = 1 : num_elem
    sub_elem_x = (-(num_sub_w-1)/2 : (num_sub_w-1)/2) * sub_width + elem_x(i);
    sub_elem_y = (-(num_sub_h-1)/2 : (num_sub_h-1)/2) * sub_height + elem_y(i);
    [sub_elem_x, sub_elem_y] = meshgrid(sub_elem_x, sub_elem_y);
    sub_elem_x = reshape(sub_elem_x, numel(sub_elem_x), 1);
    sub_elem_y = reshape(sub_elem_y, numel(sub_elem_y), 1);
    sub_elem_z = zeros(num_elem_rects, 1);
    % rect((i-1)*num_elem_rects+1 : i*num_elem_rects, :) = [i * ones(num_elem_rects,1), ...
    %     sub_elem_x - sub_width/2, sub_elem_y - sub_height/2, sub_elem_z, ...
    %     sub_elem_x + sub_width/2, sub_elem_y - sub_height/2, sub_elem_z, ...
    %     sub_elem_x + sub_width/2, sub_elem_y + sub_height/2, sub_elem_z, ...
    %     sub_elem_x - sub_width/2, sub_elem_y + sub_height/2, sub_elem_z, ...
    %     ones(num_elem_rects,1), sub_width*ones(num_elem_rects,1), sub_height*ones(num_elem_rects,1), ...
    %     sub_elem_x, sub_elem_y, sub_elem_z];
    rect((i-1)*num_elem_rects+1 : i*num_elem_rects, :) = [i * ones(num_elem_rects,1), ...
        (sub_elem_x - sub_width/2 - elem_x(i)) * cos(angle) - (sub_elem_y - sub_height/2 - elem_y(i)) * sin(angle) + elem_x(i), (sub_elem_x - sub_width/2 - elem_x(i)) * sin(angle) + (sub_elem_y - sub_height/2 - elem_y(i)) * cos(angle) + elem_y(i), sub_elem_z, ...
        (sub_elem_x + sub_width/2 - elem_x(i)) * cos(angle) - (sub_elem_y - sub_height/2 - elem_y(i)) * sin(angle) + elem_x(i), (sub_elem_x + sub_width/2 - elem_x(i)) * sin(angle) + (sub_elem_y - sub_height/2 - elem_y(i)) * cos(angle) + elem_y(i), sub_elem_z, ...
        (sub_elem_x + sub_width/2 - elem_x(i)) * cos(angle) - (sub_elem_y + sub_height/2 - elem_y(i)) * sin(angle) + elem_x(i), (sub_elem_x + sub_width/2 - elem_x(i)) * sin(angle) + (sub_elem_y + sub_height/2 - elem_y(i)) * cos(angle) + elem_y(i), sub_elem_z, ...
        (sub_elem_x - sub_width/2 - elem_x(i)) * cos(angle) - (sub_elem_y + sub_height/2 - elem_y(i)) * sin(angle) + elem_x(i), (sub_elem_x - sub_width/2 - elem_x(i)) * sin(angle) + (sub_elem_y + sub_height/2 - elem_y(i)) * cos(angle) + elem_y(i), sub_elem_z, ...
        ones(num_elem_rects,1), sub_width*ones(num_elem_rects,1), sub_height*ones(num_elem_rects,1), ...
        sub_elem_x, sub_elem_y, sub_elem_z];
end

x_center = [elem_x', elem_y', elem_z'];

aper = xdc_rectangles(rect, x_center, focus);

end
