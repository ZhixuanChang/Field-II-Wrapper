function xdc_display_aper(aper)
% Display the aperture.
%
% Parameters
% aper - the aperture to be displayed.
%
% See also
% xdc_show, xdc_get

aper_data = xdc_get(aper);
corner_data = aper_data(11:22, :);
corner_data = reshape(corner_data, 3, 4, []);

figure();
for i = 1 : size(corner_data, 3)
    hold on;
    if aper_data(1, i) == 0
        patch(corner_data(1, :, i), corner_data(2, :, i), 'magenta', 'EdgeColor', 'black', 'LineWidth', 1);
    else
        patch(corner_data(1, :, i), corner_data(2, :, i), 'cyan', 'EdgeColor', 'black', 'LineWidth', 1);
    end
end
hold off;
axis equal tight;

end
