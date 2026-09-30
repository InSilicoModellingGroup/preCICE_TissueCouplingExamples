SetFactory("OpenCASCADE");

R = 0.5;
H = 1.0;
lc = 0.04;

Cylinder(1) = {0, 0, 0, 0, 0, H, R};

// Sync and classify surfaces by bounding-box centre z / type
s() = Surface{:};
For i In {0:#s()-1}
  bb() = BoundingBox Surface{s[i]};
  zc = 0.5*(bb[2]+bb[5]);
  dx = bb[3]-bb[0]; dy = bb[4]-bb[1]; dz = bb[5]-bb[2];
  If (dz < 0.1*H && zc > 0.9*H)
    topSurfaces[i] = s[i];
  ElseIf (dz < 0.1*H && zc < 0.1*H)
    bottomSurfaces[i] = s[i];
  Else
    wallSurfaces[i] = s[i];
  EndIf
EndFor

Physical Volume("internal") = {1};
Physical Surface("top") = topSurfaces[];
Physical Surface("bottom") = bottomSurfaces[];
Physical Surface("wall") = wallSurfaces[];

Characteristic Length{ PointsOf{ Volume{1}; } } = lc;
Mesh.Algorithm = 6;
Mesh.Algorithm3D = 1;
