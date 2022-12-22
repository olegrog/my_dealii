H = 0.8;  // height of callus
L = 1.3;  // lenght of callus
D = 0.2;  // depth of fracture
W = 0.2;  // width of fracture
C = 1;    // mesh length scale

Point(1) = {0, 0, 0};
Point(2) = {L, 0, 0};
Point(3) = {0, H, 0};
Point(4) = {0, -D, 0};
Point(5) = {W, -D, 0};
Point(6) = {W, 0, 0};

Ellipse(1) = {3, 1, 2, 2};
Line(2) = {2, 6};
Line(3) = {6, 5};
Line(4) = {5, 4};
Line(5) = {4, 3};

Line Loop(1) = {1:5};
Plane Surface(1) = {1};

Mesh.Algorithm = 8;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = C;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;

Show "*";
