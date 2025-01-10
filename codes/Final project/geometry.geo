SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 4.00, 4.00};
Disk(2) = {0, 0, 0, 0.50};
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }
Physical Surface("Plate") = {1};
Physical Curve("Boundary") = {4, 3, 2, 1};
