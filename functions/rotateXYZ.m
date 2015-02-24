% ax, ay, az are rotation angles in degrees

function [X2r Y2r Z2r Rtot] = rotateXYZ(X1r, Y1r, Z1r, ax, ay, az)

Rx = [1 0 0; 0 cosd(ax) -1*sind(ax); 0 sind(ax) cosd(ax)];
Ry = [cosd(ay) 0 sind(ay); 0 1 0; -1*sind(ay) 0 cosd(ay)];
Rz = [cosd(az) -1*sind(az) 0; sind(az) cosd(az) 0; 0 0 1];

Rtot = Rz*Ry*Rx

coordsOut = Rtot*[X1r; Y1r; Z1r];

X2r = coordsOut(1,:);
Y2r = coordsOut(2,:);
Z2r = coordsOut(3,:);



