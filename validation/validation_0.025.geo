sz = 0.1;

//+
Point(1) = {0, 0, 0, sz};
//+
Point(2) = {0, 1, 0, sz};
//+
Point(3) = {2, 1, 0, sz};
//+
Point(4) = {2, 0, 0, sz};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {2, 3, 4, 1};
//+
Plane Surface(1) = {1};
//+
Physical Curve("LEFT", 1) = {1};
//+
Physical Curve("UPPER", 2) = {2};
//+
Physical Curve("RIGHT", 3) = {3};
//+
Physical Curve("LOWER", 4) = {4};
//+
