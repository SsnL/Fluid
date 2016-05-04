#include "CGL/CGL.h"
#include "CGL/Vector3D.h"
#include <vector>

namespace CGL{

struct TRIANGLE {
   Vector3D p[3];
};

struct GRIDCELL {
   Vector3D p[8];
   double val[8];
};

std::vector<TRIANGLE *> polygonise(GRIDCELL grid,double isolevel);
Vector3D VertexInterp(double isolevel,Vector3D p1,Vector3D p2, double valp1, double valp2);

} // namespace CGL
