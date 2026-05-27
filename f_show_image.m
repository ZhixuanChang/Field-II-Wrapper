function f_show_image(img)
% Plot 1D, 2D or 3D image.

img_sq = squeeze(img);
nd = ndims(img_sq);

if isvector(img_sq)
    figure;
    plot(img_sq);
elseif nd == 2
    figure;
    imagesc(img_sq);
elseif nd == 3
    volumeViewer(img_sq);
else
    warning('f_show_image: unsupported dimension %d.', nd);
end