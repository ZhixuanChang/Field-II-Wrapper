function scan_pos = scan_sun(dt,T)
%SCAN_SUN calculates the scanning position of the planet scanning mechanism.
%
% INPUT
%   dt - step time. The angular speed is fixed to 6e4 [deg/s].
%   T - scanning elapse.
%
% OUTPUT
%   scan_pos - the scanning position of the probe, which is a matrix with a size of [N, 2]. N is the number of scanning
%              points. Each row is a coordinate. The first column is the x-coordinate, while the second column is the
%              y-coordinate.

r = 13; % mm
d = 3.0;  % mm
D = 2.5;  % mm

omega_input = 10000; % rpm 逆时针旋转
omega = omega_input * 2 * pi / 60;

omega_d = (omega*d/2) * ((r-d)/(r-d/2))/d;
omega_D = -omega_d * d / (r-d);

% t = 0:0.0001:1;
t = 0:dt:T;
x = d*cos(omega_d*t) + D*cos(omega_D*t);
y = d*sin(omega_d*t) + D*sin(omega_D*t);

% figure;
% plot(x,y);
% hold on;

% x_sample = x(1:100:end);
% y_sample = y(1:100:end);
% plot(x_sample,y_sample,LineStyle="none",Marker="*",MarkerSize=3);
% axis equal;
% xlabel("x (mm)");
% ylabel("y (mm)");
% xlim([-7,7]);
% ylim([-7,7]);

scan_pos = [x', y'];
