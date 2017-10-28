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

#ifndef G2O_QPLANE3D_H_
#define G2O_QPLANE3D_H_

#include "g2o_types_slam3d_addons_api.h"
#include "g2o/stuff/misc.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace g2o {
  
  class G2O_TYPES_SLAM3D_ADDONS_API QPlane3D {
  public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
      friend QPlane3D operator*(const Eigen::Isometry3d& t, const QPlane3D& plane);

      QPlane3D() {
          Eigen::Vector4d v;
          v << 1., 0., 0., -1.;
          fromVector(v);
      }

      QPlane3D(const Eigen::Vector4d& v) {
          fromVector(v);
      }


      inline Eigen::Vector4d toVector() const {
          return _coeffs;
      }

      inline const Eigen::Vector4d& coeffs() const { return _coeffs; }

      inline void fromVector(const Eigen::Vector4d& coeffs_) {
          _coeffs = coeffs_;
          normalize(_coeffs);
      }

      static double azimuth(const Eigen::Vector3d& v) {
          return atan2(v(1), v(0));
      }

      static  double elevation(const Eigen::Vector3d& v) {
          return atan2(v(2), v.head<2>().norm());
      }

      double distance() const {
          double n = _coeffs.head<3>().norm();
          return -_coeffs(3) * (1. / n);
      }

      Eigen::Vector3d normal() const {
          double n = _coeffs.head<3>().norm();
          return _coeffs.head<3>() * (1. / n);
      }


      static Eigen::Matrix3d rotation(const Eigen::Vector3d& v) {
          double _azimuth = azimuth(v);
          double _elevation = elevation(v);
          return (Eigen::AngleAxisd(_azimuth, Eigen::Vector3d::UnitZ())*Eigen::AngleAxisd(-_elevation, Eigen::Vector3d::UnitY())).toRotationMatrix();
      }

      inline void oplus(const Eigen::Vector3d& v) {
          double theta = v.norm();
          // from on kaess' plane operations: Simultaneous localization and mapping with infinite planes
          double S;
          if (theta < 0.0001) {
              S = 0.5 + theta*theta / 48.;
          }
          else {
              S = sin(0.5*theta) / theta;
          }
          double C = cos(0.5*theta);
          Eigen::Quaterniond q(_coeffs);
          Eigen::Quaterniond dq(C, S*v(0), S*v(1), S*v(2));
          Eigen::Quaterniond newq = dq * q;
          _coeffs = newq.coeffs();
          normalize(_coeffs);
      }

      inline Eigen::Vector3d ominus(const QPlane3D& plane) {
          // DN: equivalent to log map formula:
          Eigen::Quaterniond qEstimate(_coeffs);
          Eigen::Quaterniond qMeasured(plane.coeffs());
          Eigen::Quaterniond deltaQ = qEstimate.conjugate() * qMeasured;
          Eigen::AngleAxisd deltaAngleAxis(deltaQ);
          if (deltaAngleAxis.angle() > M_PI) deltaAngleAxis.angle() -= 2.*M_PI;
          if (deltaAngleAxis.angle() < -M_PI) deltaAngleAxis.angle() += 2.*M_PI;
          Eigen::Vector3d error = deltaAngleAxis.axis() * deltaAngleAxis.angle();
          return error;
      }

      //protected:

      static inline void normalize(Eigen::Vector4d& coeffs) {
          double norm = coeffs.norm();
          coeffs = coeffs * (1. / norm);
          if (coeffs(3) < 0) coeffs *= -1;
      }

      Eigen::Vector4d _coeffs;
  };

  inline QPlane3D operator*(const Eigen::Isometry3d& t, const QPlane3D& plane) {
      Eigen::Vector4d v = plane._coeffs;
      Eigen::Vector4d v2;
      Eigen::Matrix3d R = t.rotation();
      v2.head<3>() = R*v.head<3>();
      v2(3) = v(3) - t.translation().dot(v2.head<3>());
      return QPlane3D(v2);
  };

}

#endif
