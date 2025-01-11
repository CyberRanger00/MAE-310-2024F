SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 2.00, 2.00};
Disk(2) = {0, 0, 0, 0.50};
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }
Mesh.CharacteristicLengthMax = 1.54;
Physical Surface("Plate") = {1};
Physical Curve("Boundary") = {3, 2, 1};
