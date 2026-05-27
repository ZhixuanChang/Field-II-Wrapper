function scat_pos = gen_sphere_solid(center, radius, math_size)
    pos_range = [center' - radius * 2, center' + radius * 2];
    nx = round(diff(pos_range(1,:)) / math_size);
    ny = round(diff(pos_range(2,:)) / math_size);
    nz = round(diff(pos_range(3,:)) / math_size);

    x_vec = (-(nx-1)/2 : (nx-1)/2) * math_size + center(1);
    y_vec = (-(ny-1)/2 : (ny-1)/2) * math_size + center(2);
    z_vec = (-(nz-1)/2 : (nz-1)/2) * math_size + center(3);
    [x_mesh, y_mesh, z_mesh] = ndgrid(x_vec, y_vec, z_vec);

    center_dis = sqrt((x_mesh - center(1)).^2 + (y_mesh - center(2)).^2 + (z_mesh - center(3)).^2);
    mask = center_dis < radius;
    scat_x = x_mesh(mask);
    scat_y = y_mesh(mask);
    scat_z = z_mesh(mask);

    scat_pos = [scat_x, scat_y, scat_z];
end
