function aper = xdc_concatenate(aper_set, focus)
% Concatenate multiple apertures.
%
% The concatenated aperture contains all physical elements of all apertures in aper_set. The physical number of each
% element will increment by the order of aper_set and internal physical number.
% 
% Remember to manually free apertures in aper_set.
%
% Parameters
% aper_set - a set of apertures to be concatenated. Each aperture in this set must consist of rectangles mathematical
% elements only.
% focus - the fixed focus point of the concatenated aperture, should be with a shape of [x, y, z].
%
% Examples
% aper = xdc_concatenate(aper_set);
%
% See also
% xdc_retangles, xdc_piston, xdc_single_rect, xdc_matrix_array
%
% Author: Zhixuan Chang (https://zhixuanchang.github.io/)
% Date: Nov. 6, 2025
% Version: 0.1

aper_set = reshape(aper_set, numel(aper_set), 1);
rect = [];
center = [];
phy_ind = 1;
for i = 1 : length(aper_set)
    aper_elem = xdc_get(aper_set(i));
    center = [center; aper_elem(24:26,1)'];
    % the physical number of elements returned by xdc_get is start from 0, so the first column in rect should be
    % incremented by 1.
    rect = [rect;
            [aper_elem(1,:)'+phy_ind, aper_elem(11:22,:)', aper_elem(5,:)', aper_elem(3:4,:)', aper_elem(8:10,:)']];
    phy_ind = max(rect(:, 1), [], 'all') + 1;
end

aper = xdc_rectangles(rect, center, focus);

end
