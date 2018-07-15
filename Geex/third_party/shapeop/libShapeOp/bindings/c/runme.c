#include "API.h"
#include <stdio.h>

//print array
void printArray( ShapeOpScalar *a, int n) {
  int r, c;
  for (c = 0; c < 3; c++) {
    for (r = 0; r < n; r++)
      printf("%f ", a[r * 3 + c]);
    printf("\n");
  }
}

int main() {
  ShapeOpScalar p[12]; //column major
  p[0] = 0.; p[3] = 0.5; p[6] = 0.5; p[9] = 0.;
  p[1] = 0.; p[4] = 0.; p[7] = 1.; p[10] = 1.;
  p[2] = 0.; p[5] = 1.; p[8] = 0.; p[11] = 1.;
  printf("Input points:\n");
  printArray(p, 4);
  struct ShapeOpSolver *s = shapeop_create();
  shapeop_setPoints(s, p, 4);
  ShapeOpScalar weight = 1.0;
  //add a plane constraint to all the vertices.
  {
    int ids[4];
    ids[0] = 0; ids[1] = 1; ids[2] = 2; ids[3] = 3;
    shapeop_addPlaneConstraint(s, ids, 4, weight);
  }
  //add a closeness constraint to the 1st vertex.
  {
    shapeop_addClosenessConstraint(s, 0, weight);
  }
  //add a closeness constraint to the 4th vertex.
  {
    shapeop_addClosenessConstraint(s, 3, weight);
  }
  shapeop_init(s);
  shapeop_solve(s, 10);
  shapeop_getPoints(s, p, 4);
  shapeop_delete(s);
  printf("Output points:\n");
  printArray(p, 4);
  return 0;
}
