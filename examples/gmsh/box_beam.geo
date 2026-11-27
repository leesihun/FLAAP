// Solid aluminium cantilever bar for FLAAP (gmsh script, OpenCASCADE, metres)
// 1.0 m (y) x 0.10 m (x) x 0.02 m (z); unstructured tetrahedra.
// Use mesh_order 2 (tet10) and clamp_group Root.
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 0.10, 1.0, 0.02};
Physical Surface("Root") = {3};     // the y = 0 face (OCC box face order: x-, x+, y-, y+, z-, z+)
Physical Volume("Bar") = {1};
Mesh.MeshSizeMax = 0.025;
