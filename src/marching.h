
namespace CGL{
typedef struct {
   Vector3D p[3];
} TRIANGLE;

typedef struct {
   Vector3D p[8];
   double val[8];
} GRIDCELL;

int Polygonise(GRIDCELL grid,double isolevel,TRIANGLE *triangles);
Vector3D VertexInterp(double isolevel,Vector3D p1,Vector3D p2, double valp1, double valp2);

}// namespace CGL