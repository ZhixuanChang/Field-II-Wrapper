function aper = xdc_matrix_array(nx, pitch_x, width_x, ny, pitch_y, width_y, n_sub_x, n_sub_y, focus)

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
