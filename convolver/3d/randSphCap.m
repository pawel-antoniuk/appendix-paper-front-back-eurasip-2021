% Paweł Antoniuk 2021
% Bialystok University of Technology

function [azel] = randSphCap(N, direction, capAngle)

capAngleRad = deg2rad(capAngle);
directionRad = deg2rad(direction);

z = (1 - cos(capAngleRad)).* rand(1, N) + cos(capAngleRad);
az = rand(1, N) * 2 * pi;
el = atan2(z, sqrt(1-z.^2));

[x,y,z] = sph2cart(az, el, 1);

rotx = @(t) [1 0 0; 0 cos(t) -sin(t); 0 sin(t) cos(t)];
rotz = @(t) [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1];

xyz = rotz(-directionRad(1)-pi/2) ...
    * rotx(directionRad(2)-pi/2) ...
    * [x; y; z];

[az,el] = cart2sph(xyz(1, :), xyz(2, :), xyz(3, :));
azel = rad2deg([az; el]');
azel(1, :) = wrapTo360(azel(1, :));
azel(2, :) = wrapTo180(azel(2, :));

% y = rotx(directionRad(2)) * [x; y; z];

% y=cos(el).*sin(az);
% x=cos(el).*cos(az);
% z=sin(el);

% a = 2 * pi * rand(1, N);
% r = rand(1, N) + rand(1, N);
% r(r > 1) = 2 - r(r > 1);
% az = sin(a) .* r .* capAngle;
% el = cos(a) .* r .* capAngle;

% z = rand(1, N) * (1 - cos(capAngle)) + cos(capAngle);
% phi = rand(1, N) * 2 * pi;
% x = sqrt(1-z.^2).*cos(phi);
% y = sqrt(1-z.^2).*sin(phi);
% az = atan2(y, x);
% el = atan2(z, sqrt(x.^2 + y.^2));
% azElPairs = [az;el]';

% az = capAngle * (2 * rand(1, N) - 1);
% el = (capAngle.^2 - az.^2) .^ (1/2) .* (2 * rand(1, N) - 1);
% azElPairs = [az;el]';