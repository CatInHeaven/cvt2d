#include "Solver.h"
#include "Constraint.h"
#include <iostream>
int main() {
  ShapeOp::Matrix34 p; //column major
  p << 0., 0.5, 0.5, 0.,
       0., 0., 1., 1.,
       0., 1., 0., 1.;
  std::cout <<  "Input points:" << std::endl;
  std::cout << p << std::endl;
  ShapeOp::Solver s;
  s.setPoints(p);
  ShapeOp::Scalar weight = 1.0;
  //add a plane constraint to all the vertices.
  {
    std::vector<int> id_vector;
    id_vector.push_back(0); id_vector.push_back(1); id_vector.push_back(2); id_vector.push_back(3);
    auto c = std::make_shared<ShapeOp::PlaneConstraint>(id_vector, weight, s.getPoints());
    s.addConstraint(c);
  }
  //add a closeness constraint to the 1st vertex.
  {
    std::vector<int> id_vector;
    id_vector.push_back(0);
    auto c = std::make_shared<ShapeOp::ClosenessConstraint>(id_vector, weight, s.getPoints());
    s.addConstraint(c);
  }
  //add a closeness constraint to the 4th vertex.
  {
    std::vector<int> id_vector;
    id_vector.push_back(3);
    auto c = std::make_shared<ShapeOp::ClosenessConstraint>(id_vector, weight, s.getPoints());
    s.addConstraint(c);
  }
  s.initialize();
  s.solve(10);
  p = s.getPoints();
  std::cout << "Output points:" << std::endl;
  std::cout << p << std::endl;
  return 0;
}
