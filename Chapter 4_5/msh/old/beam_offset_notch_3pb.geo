// Gmsh project created on Sat Jul 16 22:44:54 2022

// length of beam
L = 375;
// span of beam (between supports)
S = 330;
//width of beam
b = 100;
// offset length
of = 65;

//crack width 2w
w = .5;
//crack  length a
a = 19;

// mesh size at points
h = 2*w; 
h1 = 10;

// lenght of loading line
le = 2;

//size of elements at critical zone
dx = 5/5;

//+
Point(1) = {0, 0, 0, 1.3*h1};
//+
Point(2) = {.5*L-.5*S-.5*le, 0, 0, .3*le};
//+
Point(3) = {.5*L-.5*S+.5*le, 0, 0, .3*le};
//+
Point(4) = {.5*L-of-w, 0, 0, .5*h1};
//+
Point(5) = {.5*L-of-w,a,0,h};
//+
Point(6) = {.5*L-of+w,a,0,h};
//+
Point(7) = {.5*L-of+w,0,0,.5*h1};
//+
Point(8) = {.5*L+.5*S-.5*le,0,0,.3*le};
//+
Point(9) = {.5*L+.5*S+.5*le,0,0,.3*le};
//+
Point(10) = {L,0,0,1.3*h1};
//+
Point(11) = {L,b,0,1.5*h1};
//+
Point(12) = {.5*L+le,b,0,.35*le};
//+
Point(13) = {.5*L-le,b,0,.35*le};
//+
Point(14) = {0,b,0,1.5*h1};

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
Line(8) = {8, 9};
//+
Line(9) = {9, 10};
//+
Line(10) = {10, 11};
//+
Line(11) = {11,12};
//+
Line(12) = {12, 13};
//+
Line(13) = {13, 14};
//+
Line(14) = {14,1};
//+



Physical Curve(101) = {2};
//+
Physical Curve(102) = {8};
//+ 
Physical Curve(103) = {12};
//+
Physical Curve(1) = {1,  3, 4, 5, 6, 7, 9, 10, 11, 13, 14};
//+


Curve Loop(1) = {1,2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,13,14};
//+
Plane Surface(1) = {1};
//+
Physical Surface(1000) = {1};
//+

//thickness of critical zone
t = 60;
Field[1] = Box;
//+
Field[1].Thickness = t;
//+
Field[1].VIn = dx;
//+
Field[1].VOut = 35;
//+
Field[1].XMax = 180;
//+
Field[1].XMin = 94;
//+
Field[1].YMax = 100;
//+
Field[1].YMin = a/4;
//+
Background Field = 1;
