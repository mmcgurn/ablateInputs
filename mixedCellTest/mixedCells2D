//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, .01, 0, 1.0};
//+
Point(3) = {.01, .01, 0, 1.0};
//+
Point(4) = {.01, 0, 0, 1.0};
//+
Point(5) = {.02, 0, 0, 1.0};
//+
Point(6) = {.03, 0, 0, 1.0};
//+
Point(7) = {.03, 0.01, 0, 1.0};
//+
Point(8) = {.02, 0.01, 0, 1.0};
//+
Line(1) = {2, 1};
//+
Line(2) = {1, 4};
//+
Line(3) = {4, 5};
//+
Line(4) = {5, 6};
//+
Line(5) = {6, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 3};
//+
Line(8) = {3, 2};
//+
Line(9) = {3, 4};
//+
Line(10) = {8, 5};
//+
Curve Loop(1) = {1, 2, -9, 8};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, 3, -10, 7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {10, 4, 5, 6};
//+
Plane Surface(3) = {3};
//+
Physical Curve("inlet", 11) = {1};
//+
Physical Curve("outlet", 12) = {5};
//+
Physical Curve("wall", 13) = {8, 7, 6, 4, 3, 2};
//+
Physical Surface("flow", 14) = {1, 2, 3};
//+
Recombine Surface {1, 3};
