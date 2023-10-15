// Gmsh project created on Thu Aug 11 11:07:54 2022

L = 5; 
h = 1;
//mesh size dx
dx = h;
//+
Point(1) = {0, 0, 0, dx};
//+
Point(2) = {5, 0, 0, dx};
//+
Point(3) = {-5, 0, 0, dx};
//+
Point(4) = {10, 0, 0, dx};
//+
Point(5) = {-10, 0, 0, dx};
//+
Circle(1) = {4, 1, 5};
//+
Circle(2) = {5, 1, 4};
//+
Circle(3) = {2, 1, 3};
//+
Circle(4) = {3, 1, 2};
//+
Physical Curve(100) = {2, 1, 3, 4};
Curve Loop(1) = {1,2};
Curve Loop(2) = {3,4};
Plane Surface(1) = {1,2};
//+
Physical Surface(1000) = {1};