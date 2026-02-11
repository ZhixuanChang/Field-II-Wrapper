function scanning_path = scan_zigzag(nx, ny, sx, sy)
% SCAN_ZIGZAG generate a zigzag scanning path.
%
% INPUT
%   nx - the number of scanning points along x-axis.
%   ny - the number of scanning points along y-axis.
%   sx - the scanning step size along x-axis in unit of [m].
%   sy - the scanning step size along y-axis in unit of [m].
%
% OUTPUT
%   scanning_path - the coordinates of the scanning points in size of [N, 3] where N is the number of scanning points in
%                   total. N is always equal to nx * ny. The first column is the x-coordinate, the second column is the
%                   y-coordinate, and the last column is the z-coordinate (fixed to 0 since the normal zigzag scanning
%                   is implemented in a plane).
%
% USAGE
%   nx = 100;
%   ny = 100;
%   sx = 0.5e-3;
%   sy = 0.5e-3;
%   scanning_path = scan_zigzag(nx, ny, sx, sy);
%
% See also

% Initialize output matrix
N = nx * ny;
scanning_path = zeros(N, 3);

% Generate x coordinates with proper spacing
% Center the grid at origin
x_start = -(nx - 1) / 2 * sx;
x_coords = x_start + (0:nx-1) * sx;

% Generate y coordinates with proper spacing
% Center the grid at origin
y_start = -(ny - 1) / 2 * sy;
y_coords = y_start + (0:ny-1) * sy;

% Generate zigzag pattern
idx = 1;
for iy = 1:ny
    if mod(iy, 2) == 1
        % Odd rows: left to right
        for ix = 1:nx
            scanning_path(idx, 1) = x_coords(ix);
            scanning_path(idx, 2) = y_coords(iy);
            scanning_path(idx, 3) = 0;  % z-coordinate is always 0
            idx = idx + 1;
        end
    else
        % Even rows: right to left (zigzag)
        for ix = nx:-1:1
            scanning_path(idx, 1) = x_coords(ix);
            scanning_path(idx, 2) = y_coords(iy);
            scanning_path(idx, 3) = 0;  % z-coordinate is always 0
            idx = idx + 1;
        end
    end
end

end
