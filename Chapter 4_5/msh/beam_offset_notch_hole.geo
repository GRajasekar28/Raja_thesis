// Gmsh project created on Sat Jul 16 22:44:54 2022

//crack width 2w
w = .03;
//crack  length a
a = 1;

// mesh size at points
h = 2*w; 
h1 = 1;

// lenght of loading line
le = .1;

//+
Point(1) = {0, 0, 0, 1.3*h1};
//+
Point(2) = {1-.5*le, 0, 0, .3*le};
//+
Point(3) = {1+.5*le, 0, 0, .3*le};
//+
Point(4) = {4-w, 0, 0, .5*h1};
//+
Point(5) = {4-w,a,0,h};
//+
Point(6) = {4+w,a,0,h};
//+
Point(7) = {4+w,0,0,.5*h1};
//+
Point(8) = {19-.5*le,0,0,.3*le};
//+
Point(9) = {19+.5*le,0,0,.3*le};
//+
Point(10) = {20,0,0,1.3*h1};
//+
Point(11) = {20,8,0,1.5*h1};
//+
Point(12) = {10+le,8.,0,.35*le};
//+
Point(13) = {10-le,8.,0,.35*le};
//+
Point(14) = {0,8,0,1.5*h1};

r1 = .25;
Point(15) = {6,8-1.25,0,.35*h1};
Point(16) = {6+r1,8-1.25,0,.35*h1};
Point(17) = {6-r1,8-1.25,0,.35*h1};
Point(18) = {6,8-3.25,0,.35*h1};
Point(19) = {6+r1,8-3.25,0,.35*h1};
Point(20) = {6-r1,8-3.25,0,.35*h1};
Point(21) = {6,8-5.25,0,.35*h1};
Point(22) = {6+r1,8-5.25,0,.35*h1};
Point(23) = {6-r1,8-5.25,0,.35*h1};





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
//+
Circle(15) = {16, 15, 17};
//+
Circle(16) = {17, 15, 16};
//+
Circle(17) = {19, 18, 20};
//+
Circle(18) = {20, 18, 19};
//+
Circle(19) = {22, 21, 23};
//+
Circle(20) = {23, 21, 22};



Physical Curve(101) = {2};
//+
Physical Curve(102) = {8};
//+ 
Physical Curve(103) = {12};
//+
Physical Curve(1) = {1,  3, 4, 5, 6, 7, 9, 10, 11, 13, 14,15,16,17,18,19,20};
//+

Curve Loop(1) = {15,16};
Curve Loop(2) = {17,18};
Curve Loop(3) = {19,20};


Curve Loop(4) = {1,2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,13,14};
//+
Plane Surface(1) = {4,-1,-2,-3};
//+
Physical Surface(1000) = {1};
//+
Field[1] = Box;
//+
Field[1].Thickness = 5;
//+
Field[1].VIn = w;
//+
Field[1].VOut = h1;
//+
Field[1].XMax = 4+4;
//+
Field[1].XMin = 4-1;
//+
Field[1].YMax = 8;
//+
Field[1].YMin = 0;
//+
Background Field = 1;

