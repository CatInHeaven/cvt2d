///////////////////////////////////////////////////////////////////////////////
// This file is part of ShapeOp, a lightweight C++ library
// for static and dynamic geometry processing.
//
// Copyright (C) 2014 Sofien Bouaziz <sofien.bouaziz@gmail.com>
// Copyright (C) 2014 LGG EPFL
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
///////////////////////////////////////////////////////////////////////////////
#include "window.h"
#include "renderer.h"
///////////////////////////////////////////////////////////////////////////////
int main() {
  float width = 1024.0f;
  float height = 768.0f;
  GLWindow w("ShapeOp Demo", width, height);
  GLCHECK //The first error is due to glewInit!
  GLRenderer r(width, height);
  GLCHECK
  for (int i = -1; i <= 1; ++i)
    for (int j = -1; j <= 1; ++j) {
      //for (int k = -3; k < 3; ++k) {
        std::cout << i << " " << j << std::endl;
        auto m = std::make_shared<GLSurfaceMeshObject>("neutral.obj", Eigen::Vector3f((i+1.0)/2.0, (j+1.0)/2.0, 1.0f));
        Eigen::Matrix4f t = Eigen::Matrix4f::Identity();
        t(0, 3) = 0.7*i;
        t(1, 3) = 0.0;
        t(2, 3) = 0.7*j;
        m->setModel(t);
        r.addMesh(m);
        GLCHECK
      }
//  {
//    auto m = std::make_shared<GLSurfaceMeshObject>("neutral.obj", Eigen::Vector3f(1.0f, 1.0f, 1.0f));
//    Eigen::Matrix4f t = Eigen::Matrix4f::Identity();
//    t(0, 3) = 0.0f;
//    m->setModel(t);
//    r.addMesh(m);
//    GLCHECK
//  }
//  {
//    auto m = std::make_shared<GLSurfaceMeshObject>("neutral.obj", Eigen::Vector3f(0.0f, 1.0f, 0.0f));
//    Eigen::Matrix4f t = Eigen::Matrix4f::Identity();
//    t(0, 3) = -0.7f;
//    m->setModel(t);
//    r.addMesh(m);
//    GLCHECK
//  }
//  {
//    auto m = std::make_shared<GLSurfaceMeshObject>("neutral.obj", Eigen::Vector3f(0.0f, 0.0f, 1.0f));
//    Eigen::Matrix4f t = Eigen::Matrix4f::Identity();
//    t(0, 3) = 0.7f;
//    m->setModel(t);
//    r.addMesh(m);
//    GLCHECK
//  }
  {
    auto m = std::make_shared<GLSurfaceMeshObject>("cube.obj", Eigen::Vector3f(1.0f, 1.0f, 1.0f));
    Eigen::Matrix4f t = 100.0 * Eigen::Matrix4f::Identity();
    t(1, 3) = -50.2;
    t(3, 3) = 1.0;
    m->setModel(t);
    r.addMesh(m);
    GLCHECK
  }
//  {
//    auto m = std::make_shared<GLSurfaceMeshObject>("cube.obj", Eigen::Vector3f(1.0f, 1.0f, 1.0f));
//    Eigen::Matrix4f t = 100.0 * Eigen::Matrix4f::Identity();
//    t(3, 3) = 1.0;
//    m->setModel(t);
//    r.addMesh(m);
//    GLCHECK
//  }
  auto c = std::make_shared<Camera>(Eigen::Vector3f(0.0f, 1.0f, 2.0f), Eigen::Vector3f::Zero(), width / height);
  GLCHECK
  r.setCamera(c);
  auto l = std::make_shared<Camera>(Eigen::Vector3f(0.0f, 5.0f, 2.0f), Eigen::Vector3f::Zero(), width / height);
  GLCHECK
  r.setLight(l);
  auto d = std::bind(&GLRenderer::display, r);
  w.loop(d, c);
  return 0;
}
///////////////////////////////////////////////////////////////////////////////
