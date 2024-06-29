Merge '../common/profile_blade_se1050_4000p.geo';
Merge '../common/calc.geo';

s2 = 1;

BSpline(1) = {2145:2070};
BSpline(2) = {2070:9};
Line(3) = {9, 3992};
BSpline(4) = {3992:2145};

  pfirst = newp;
Point(pfirst) = {xmin - pitch, 0.01, 0, s2};
Point(newp) = {xmin, 0.025, 0, s2};
Point(newp) = {0.5*(xmin+xmax), 0.03, 0, s2};
Point(newp) = {xmax, -0.035, 0, s2};
  plast = newp;
Point(plast) = {xmax + pitch, -0.08, 0, s2};

Translate {0, -pitch, 0} {
  Duplicata { Point{pfirst:plast}; }
}

//HRANICE
Line(5) = {4006, 4007};
BSpline(6) = {4007, 4008, 4009};
Line(7) = {4009, 4010};
Line(8) = {4010, 4005};
Line(9) = {4005, 4004};
BSpline(10) = {4004, 4003, 4002};
Line(11) = {4002, 4001};
Line(12) = {4001, 4006};

Line Loop(1) = {5, 6, 7, 8, 9, 10, 11, 12};
Line Loop(2) = {2, 3, 4, 1};
Plane Surface(1) = {1, 2};

dz = -0.1*pitch;
Extrude {0, 0, dz} {
  Surface{1};
  Layers{1};
  Recombine;
}

Physical Surface("inlet") = {57};
Physical Surface("outlet") = {41};
Physical Surface("wall") = {61, 73, 69, 65};
Physical Surface("periodbeg") = {29, 33, 37};
Physical Surface("periodend") = {53, 49, 45};
Physical Surface("wall2") = {1, 74};
Physical Volume("interior") = {1};

//INLET, OUTLET
Transfinite Line {12} = 20 Using Progression 1;
Transfinite Line {8} = 40 Using Progression 1;

//BLADE
Transfinite Line {1} = 12 Using Bump 2;
Transfinite Line {2} = 120 Using Bump 0.2;
Transfinite Line {3} = 4 Using Bump 1;
Transfinite Line {4} = 120 Using Bump 0.2;

//PERIOD
Transfinite Line {5} = 25 Using Progression 1;
Transfinite Line {6} = 80 Using Progression 0.985;
Transfinite Line {7} = 60 Using Progression 1.01;

Periodic Line {11} = {5} Translate {0, pitch, 0};
Periodic Line {10} = {6} Translate {0, pitch, 0};
Periodic Line {9} = {7} Translate {0, pitch, 0};

//Merge './settings.geo';
