// Gmsh project created on Sat Feb 26 12:21:56 2022
//+

//Radius and length between supports
R = 75;
s = 120;

// Notch length and width (mm)
a = 10;   
w = .35;

// 2*l_c
lc = 20;  

// mesh size at selected points
m1 = 8;
m2 = 8;

Point(1) = {0, 0, 0, m2};
//+
Point(2) = {R-.5*s, 0, 0, m1};
//+
Point(3) = {R-w/2, 0, 0, m1};
//+
Point(4) = {R-w/2, a, 0, 2*w};
//+
Point(5) = {R+w/2, a, 0, 2*w};
//+
Point(6) = {R+w/2, 0, 0, m1};
//+
Point(7) = {R+s/2, 0, 0, m1};
//+
Point(8) = {2*R, 0, 0, m2};
//+

Point(9) = {R, R, 0, m1};
//+  
Point(10) = {R, 0, 0, 1.0};
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
Line(7) = {7, 8};
//+
Circle(8) = {8, 10, 9};
//+
Circle(9) = {9, 10, 1};
//+
Curve Loop(1) = {1, 2,3,4,5,6,7,8,9};
//+
Plane Surface(1) = {1};
//+
Physical Point(12) = {2};
//+
Physical Point(13) = {7};
//+
Physical Point(26) = {9};
//+
Physical Curve(27) = {1,2,3,4,5,6,7,8,9};
//+
Physical Surface(28) = {1};
//+


Field[1] = Box;
//+
Field[1].Thickness = lc;
//+
Field[1].VIn = 3*w;
//+
Field[1].VOut = 8;
//+
Field[1].XMax = R+lc/2;
//+
Field[1].XMin = R-lc/2;
//+
Field[1].YMin = a/3;
Field[1].YMax = .95*R;
//+
Background Field = 1;
