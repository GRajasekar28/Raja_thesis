// Gmsh project created on Sat Feb 26 12:21:56 2022
//+

//Radius and length between supports (mm)
R = 75;
s = 120;

// Notch length and width (mm)
a = 15;   
w = 0;

// l_c
lc = 4;  

// mesh size at selected points
m1 = 8;
m2 = 8;

//loading length
le = 2;
theta = le/R;
//pi = 3.14159265358979323846;
//angle = theta*pi/180;

Point(1) = {0, 0, 0, m2};
//+
Point(2) = {R-.5*s-.5*le, 0, 0, le/8};
//+
Point(3) = {R-.5*s+.5*le, 0, 0, le/8};
//
Point(4) = {R, 0, 0, m1};
//+
Point(5) = {R, a, 0, 2*w};
//+
Point(6) = {R, R, 0, le/8};
//+
Point(7) = {R-R*Sin(theta), R*Cos(theta), 0, le/8};
//



//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Circle(7) = {7, 4, 1};

//+
Curve Loop(1) = {1,2,3,4,5,6,7};
//+
Plane Surface(1) = {1};
//+
Physical Curve(101) = {2};
//+
Physical Curve(102) = {6};
//+
Physical Curve(103) = {5};
//+
Physical Curve(1) = {1,3,4,7};
//+
Physical Surface(28) = {1};
//+


Field[1] = Box;
//+
Field[1].Thickness = lc;
//+
Field[1].VIn = lc/10;
//+
Field[1].VOut = 10;
//+
Field[1].XMax = R;
//+
Field[1].XMin = R-lc;
//+
Field[1].YMin = a/4;
Field[1].YMax = R;
//+
Background Field = 1;
