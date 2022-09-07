// Gmsh project created on Tue Aug  2 10:42:23 2022
SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0, 0, 0, .1, 0, 0, .01, 2*Pi};
//+
Physical Surface("sides", 4) = {1};
//+
Physical Surface("outlet", 5) = {2};
//+
Physical Surface("inlet", 6) = {3};
