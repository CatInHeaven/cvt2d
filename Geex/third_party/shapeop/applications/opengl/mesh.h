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
#ifndef MESH_H
#define MESH_H
///////////////////////////////////////////////////////////////////////////////
#include <Types.h>
#include <GL/glew.h>
#include <GL/glfw.h>
#include <OpenGP/Surface_mesh.h>
#include <OpenGP/surface_mesh/IO.h>
///////////////////////////////////////////////////////////////////////////////
class GLSurfaceMeshObject {
 public:
  GLSurfaceMeshObject(const std::string &file, const Eigen::Vector3f &color = Eigen::Vector3f(1.0f, 0.0f, 0.0f)) : color_(color) {
    opengp::read_mesh(mesh_, file);
    model_ = Eigen::Matrix4f::Identity();
    init();
  }
  void setModel(const Eigen::Matrix4f &model) {
    model_ = model;
  }
  const Eigen::Matrix4f &getModel() const {
    return model_;
  }
  const Eigen::Vector3f &getColor() const { return color_; }
  void display() {
    //update();
    glBindVertexArray(vertexArrayID_);
    glDrawElements(GL_TRIANGLES, numIndices_, GL_UNSIGNED_INT, 0);
  }
 private:
  void update() {
    mesh_.update_vertex_normals();
    glBindVertexArray(vertexArrayID_);
    updateVertexBuffer();
    updateNormalBuffer();
  }
  void updateVertexBuffer() {
    auto vpoints = mesh_.get_vertex_property<opengp::Vec3>("v:point");
    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer_);
    glBufferData(GL_ARRAY_BUFFER, mesh_.n_vertices() * sizeof(opengp::Vec3), vpoints.data(), GL_STATIC_DRAW);
  }
  void updateNormalBuffer() {
    auto vnormals = mesh_.get_vertex_property<opengp::Vec3>("v:normal");
    glBindBuffer(GL_ARRAY_BUFFER, normalBuffer_);
    glBufferData(GL_ARRAY_BUFFER, mesh_.n_vertices() * sizeof(opengp::Vec3), vnormals.data(), GL_STATIC_DRAW);
  }
  void init() {
    mesh_.update_vertex_normals();
    /// Vertex Array
    glGenVertexArrays(1, &vertexArrayID_);
    glBindVertexArray(vertexArrayID_);

    /// Vertex Buffer
    glGenBuffers(1, &vertexBuffer_);
    updateVertexBuffer();
    /// Vertex Attribute ID
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer_); //TODO TO REMOVER
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

    /// Normal Buffer
    glGenBuffers(1, &normalBuffer_);
    updateNormalBuffer();
    /// Vertex Attribute ID
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, normalBuffer_); //TODO TO REMOVE
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

    /// Index Buffer
    std::vector<unsigned int> indices;
    for (auto fit = mesh_.faces_begin(); fit != mesh_.faces_end(); ++fit) {
      unsigned int n = mesh_.valence(*fit);
      auto vit = mesh_.vertices(*fit);
      for (unsigned int v = 0; v < n; ++v) {
        indices.push_back((*vit).idx());
        ++vit;
      }
    }
    numIndices_ = static_cast<int>(indices.size());
    GLuint triangleBuffer;
    glGenBuffers(1, &triangleBuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, triangleBuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);
  }
  opengp::Surface_mesh mesh_;
  Eigen::Matrix4f model_;
  Eigen::Vector3f color_;
  GLuint vertexArrayID_;
  GLuint vertexBuffer_;
  GLuint normalBuffer_;
  int numIndices_;
};
///////////////////////////////////////////////////////////////////////////////
class FullScreenQuad {
 public:
  FullScreenQuad() {
    init();
  }
  void display() {
    glBindVertexArray(vertexArrayID_);
    glDrawArrays(GL_TRIANGLES, 0, 6);
  }
 private:
  void init() {
    /// Vertex Array
    glGenVertexArrays(1, &vertexArrayID_);
    glBindVertexArray(vertexArrayID_);
    /// Vertex Buffer
    GLfloat quad[] = {
      -1.0f, -1.0f, 0.0f,
      1.0f, -1.0f, 0.0f,
      -1.0f,  1.0f, 0.0f,
      -1.0f,  1.0f, 0.0f,
      1.0f, -1.0f, 0.0f,
      1.0f,  1.0f, 0.0f,
    };
    GLuint vertexBuffer;
    glGenBuffers(1, &vertexBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quad) , quad, GL_STATIC_DRAW);
    /// Vertex Attribute ID
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
  }
  GLuint vertexArrayID_;
};
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
