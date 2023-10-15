// Gmsh project created on Mon Aug 22 11:24:32 2022
L = 1;
l =1;
a = .3*L;
b = .8*L;
c1 = L -2*a;
c2 = L - b;
t = b/100;
r1 = .15;
r2 = .38;
Point(1) = {0, 0, 0, l};
//+
Point(2) = {L, 0, 0, l};
//+
Point(3) = {L, L , 0, l};
//+
Point(4) = {0, L, 0, l};


Point(6) = {.5,.5,.0,l};
Point(7) = {.5-r1, 0.5, 0, l};
//+
Point(8) = {.5+r1, 0.5, 0, l};
//+
Point(9) = {.5-r2, 0.5, 0, l};
//+
Point(10) = {.5+r2, 0.5, 0, l};




//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};

Line(4) = {4,1};


Circle(5) = {8, 6, 7};
//+
Circle(6) = {7, 6, 8};
//+
Circle(7) = {10, 6, 9};
//+
Circle(8) = {9, 6, 10};





//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Curve Loop(2) = {5,6};
//+
Curve Loop(3) = {7,8};
//+

Plane Surface(1) = {1,3};

Plane Surface(2) = {2};
//+
Plane Surface(3) = {3,2};
Physical Curve(1000) ={1,2,3,4};
Physical Surface(100) = {2};
//+
Physical Surface(101) = {3};
Physical Surface(1000) = {1};
//+
//+
//+

