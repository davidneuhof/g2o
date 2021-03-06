// g2o - General Graph Optimization
// Copyright (C) 2011 R. Kuemmerle, G. Grisetti, W. Burgard
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "vertex_qplane.h"

#include "g2o/stuff/opengl_wrapper.h"

namespace g2o
{

  VertexQPlane::VertexQPlane(){
    color << .2, .2, .2;
    //color << .9, .1, .1;
  }
  
  bool VertexQPlane::read(std::istream& is) {
    Vector4D lv;
    for (int i=0; i<4; i++)
      is >> lv[i];
    setEstimate(QPlane3D(lv));
    is >> color(0) >> color(1) >> color(2);
    return true;
  }

  bool VertexQPlane::write(std::ostream& os) const {
    Vector4D lv=_estimate.toVector();
    for (int i=0; i<4; i++){
      os << lv[i] << " ";
    }
    os << color(0) << " " << color(1) << " " << color(2) << " ";
    return os.good();
  }

#ifdef G2O_HAVE_OPENGL

  VertexQPlaneDrawAction::VertexQPlaneDrawAction(): DrawAction(typeid(VertexQPlane).name())
  {
  }

  bool VertexQPlaneDrawAction::refreshPropertyPtrs(HyperGraphElementAction::Parameters* params_)
  {
    if (!DrawAction::refreshPropertyPtrs(params_))
      return false;
    if (_previousParams){
      _planeWidth = _previousParams->makeProperty<FloatProperty>(_typeName + "::PLANE_WIDTH", 3);
      _planeHeight = _previousParams->makeProperty<FloatProperty>(_typeName + "::PLANE_HEIGHT", 3);
    } else {
      _planeWidth = 0;
      _planeHeight = 0;
    }
    return true;
  }

  HyperGraphElementAction* VertexQPlaneDrawAction::operator()(HyperGraph::HyperGraphElement* element, 
                 HyperGraphElementAction::Parameters* params_)
  {
    if (typeid(*element).name()!=_typeName)
      return 0;

    refreshPropertyPtrs(params_);
    if (! _previousParams)
      return this;
    
    if (_show && !_show->value())
      return this;

    VertexQPlane* that = static_cast<VertexQPlane*>(element);
    double d = that->estimate().distance();
    double azimuth = QPlane3D::azimuth(that->estimate().normal());
    double elevation = QPlane3D::elevation(that->estimate().normal());
    glColor3f(float(that->color(0)), float(that->color(1)), float(that->color(2)));
    glPushMatrix();
    glRotatef(float(RAD2DEG(azimuth)), 0.f, 0.f, 1.f);
    glRotatef(float(RAD2DEG(elevation)), 0.f, -1.f, 0.f);
    glTranslatef(float(d), 0.f ,0.f);
    
    if (_planeWidth && _planeHeight){
      glBegin(GL_QUADS);
      glNormal3f(-1.f, 0.f, 0.f);
      glVertex3f(0.f, -_planeWidth->value(), -_planeHeight->value());
      glVertex3f(0.f,  _planeWidth->value(), -_planeHeight->value());
      glVertex3f(0.f,  _planeWidth->value(),  _planeHeight->value());
      glVertex3f(0.f, -_planeWidth->value(),  _planeHeight->value());
      glEnd();
    }

    glPopMatrix();
    return this;
  }
#endif

}
