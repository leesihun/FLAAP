// Cantilever aluminium plate wing for FLAAP (gmsh script, metres)
// 0.5 m chord (x) x 1.0 m span (y), structured 6 x 12 MITC4 quads.
// Same plate as examples/ansys and examples/calculix.
//   structure gmsh / structure_file plate_wing.geo / clamp_group Root / thickness 0.004
c = 0.5;  s = 1.0;  nc = 6;  ns = 12;
Point(1) = {0, 0, 0};  Point(2) = {c, 0, 0};
Point(3) = {c, s, 0};  Point(4) = {0, s, 0};
Line(1) = {1, 2};  Line(2) = {2, 3};  Line(3) = {3, 4};  Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Curve{1, 3} = nc + 1;
Transfinite Curve{2, 4} = ns + 1;
Transfinite Surface{1};
Recombine Surface{1};
Physical Curve("Root") = {1};      // y = 0 edge: clamped
Physical Surface("Skin") = {1};
