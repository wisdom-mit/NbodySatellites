/*

Copyright (C) 2019-2022 Jack Wisdom

This file is part of NbodySatellites.  NbodySatellites is software supporting the work
of Jack Wisdom .

NbodySatellites is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

NbodySatellites is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with nbody-satellites.  If not, see <https://www.gnu.org/licenses/>.

*/

Quaternion quaternion_multiply(Quaternion q1, Quaternion q2)
{
  Quaternion q;

  q.r = q1.r*q2.r - (q1.x*q2.x + q1.y*q2.y + q1.z*q2.z);
  q.x = q1.r*q2.x + q2.r*q1.x + q1.y*q2.z - q2.y*q1.z;
  q.y = q1.r*q2.y + q2.r*q1.y + q1.z*q2.x - q2.z*q1.x;
  q.z = q1.r*q2.z + q2.r*q1.z + q1.x*q2.y - q2.x*q1.y;

  return(q);
}

Quaternion quaternion_conjugate(Quaternion q)
{
  Quaternion r;

  r.r = q.r;
  r.x = -q.x;
  r.y = -q.y;
  r.z = -q.z;

  return(r);
}

Quaternion quaternion_normalize(Quaternion q)
{
  Quaternion r;
  double m;

  m = sqrt(q.r*q.r + q.x*q.x + q.y*q.y + q.z*q.z);

  r.r = q.r/m;
  r.x = q.x/m;
  r.y = q.y/m;
  r.z = q.z/m;

  return(r);
}

Quaternion rotation_matrix_to_quaternion(Matrix M)
{
  Quaternion q;
  double q0_2, q1_2, q2_2, q3_2;
  double q0q1, q0q2, q0q3;
  double q1q2, q1q3, q2q3;
  
  q0_2 = 0.25*(1.0 + M.xx + M.yy + M.zz);
  q1_2 = 0.25*(1.0 + M.xx - M.yy - M.zz);
  q2_2 = 0.25*(1.0 - M.xx + M.yy - M.zz);
  q3_2 = 0.25*(1.0 - M.xx - M.yy + M.zz);

  q0q1 = 0.25*(M.zy - M.yz);
  q0q2 = 0.25*(M.xz - M.zx);
  q0q3 = 0.25*(M.yx - M.xy);
  q1q2 = 0.25*(M.xy + M.yx);
  q1q3 = 0.25*(M.xz + M.zx);
  q2q3 = 0.25*(M.yz + M.zy);

  if((q0_2 > q1_2) && (q0_2 > q2_2) && (q0_2 > q3_2)) {
    q.r = sqrt(q0_2);
    q.x = q0q1/q.r;
    q.y = q0q2/q.r;
    q.z = q0q3/q.r;
  } else if ((q1_2 > q0_2) && (q1_2 > q2_2) && (q1_2 > q3_2)) {
    q.x = sqrt(q1_2);
    q.r = q0q1/q.x;
    q.y = q1q2/q.x;
    q.z = q1q3/q.x;
  } else if ((q2_2 > q0_2) && (q2_2 > q1_2) && (q2_2 > q3_2)) {
    q.y = sqrt(q2_2);
    q.r = q0q2/q.y;
    q.x = q1q2/q.y;
    q.z = q2q3/q.y;
  } else {
    q.z = sqrt(q3_2);
    q.r = q0q3/q.z;
    q.x = q1q3/q.z;
    q.y = q2q3/q.z;
  }
  return(q);
}

Matrix quaternion_to_rotation_matrix(Quaternion q) 
{
  Matrix M;
  double m2;

  m2 = q.r*q.r + q.x*q.x + q.y*q.y + q.z*q.z;

  M.xx = (q.r*q.r + q.x*q.x - q.y*q.y - q.z*q.z)/m2;
  M.xy = 2.0*(q.x*q.y - q.r*q.z)/m2;
  M.xz = 2.0*(q.x*q.z + q.r*q.y)/m2;
  M.yx = 2.0*(q.r*q.z + q.x*q.y)/m2;
  M.yy = (q.r*q.r - q.x*q.x + q.y*q.y - q.z*q.z)/m2;
  M.yz = 2.0*(q.y*q.z - q.r*q.x)/m2;
  M.zx = 2.0*(q.x*q.z - q.r*q.y)/m2;
  M.zy = 2.0*(q.r*q.x + q.y*q.z)/m2;
  M.zz = (q.r*q.r - q.x*q.x - q.y*q.y + q.z*q.z)/m2;

  return(M);
}


Matrix matrix_multiply(Matrix M1, Matrix M2)
{
  Matrix M;

  M.xx = M1.xx*M2.xx + M1.xy*M2.yx + M1.xz*M2.zx;
  M.xy = M1.xx*M2.xy + M1.xy*M2.yy + M1.xz*M2.zy;
  M.xz = M1.xx*M2.xz + M1.xy*M2.yz + M1.xz*M2.zz;

  M.yx = M1.yx*M2.xx + M1.yy*M2.yx + M1.yz*M2.zx;
  M.yy = M1.yx*M2.xy + M1.yy*M2.yy + M1.yz*M2.zy;
  M.yz = M1.yx*M2.xz + M1.yy*M2.yz + M1.yz*M2.zz;

  M.zx = M1.zx*M2.xx + M1.zy*M2.yx + M1.zz*M2.zx;
  M.zy = M1.zx*M2.xy + M1.zy*M2.yy + M1.zz*M2.zy;
  M.zz = M1.zx*M2.xz + M1.zy*M2.yz + M1.zz*M2.zz;

  return(M);
}

Matrix matrix_transpose(Matrix M1)
{
  Matrix M;

  M.xx = M1.xx;
  M.xy = M1.yx;
  M.xz = M1.zx;

  M.yx = M1.xy;
  M.yy = M1.yy;
  M.yz = M1.zy;

  M.zx = M1.xz;
  M.zy = M1.yz;
  M.zz = M1.zz;

  return(M);
}

Vector matrix_times_vector(Matrix M, Vector v)
{
  Vector r;

  r.x = M.xx*v.x + M.xy*v.y + M.xz*v.z;
  r.y = M.yx*v.x + M.yy*v.y + M.yz*v.z;
  r.z = M.zx*v.x + M.zy*v.y + M.zz*v.z;

  return(r);
}

Vector quaternion_rotate_vector(Quaternion Q, Vector v)
{
  Vector r;
  Quaternion qv, qr;
  
  qv.r = 0.0;
  qv.x = v.x;
  qv.y = v.y;
  qv.z = v.z;

  
  qr = quaternion_multiply(Q, quaternion_multiply(qv, quaternion_conjugate(Q)));

  r.x = qr.x;
  r.y = qr.y;
  r.z = qr.z;

  return(r);
}

QWState quaternion_state_to_qw_state(QuaternionState qs)
{
  QWState qw;
  double m2;

  m2 = qs.r*qs.r + qs.x*qs.x + qs.y*qs.y + qs.z*qs.z;

  qw.t = qs.t;

  qw.q.r = qs.r;
  qw.q.x = qs.x;
  qw.q.y = qs.y;
  qw.q.z = qs.z;

  qw.wa = 2.0*(qs.r*qs.xdot - qs.x*qs.rdot - qs.y*qs.zdot + qs.z*qs.ydot)/m2;
  qw.wb = 2.0*(qs.r*qs.ydot - qs.y*qs.rdot + qs.x*qs.zdot - qs.z*qs.xdot)/m2;
  qw.wc = 2.0*(qs.r*qs.zdot - qs.z*qs.rdot - qs.x*qs.ydot + qs.y*qs.xdot)/m2;

  return(qw);
}

QuaternionState qw_state_to_quaternion_state(QWState qw)
{
  Quaternion q;
  QuaternionState qs;

  qs.t = qw.t;

  q = qw.q;

  qs.r = q.r;
  qs.x = q.x;
  qs.y = q.y;
  qs.z = q.z;

  qs.rdot = (- qw.wa*q.x - qw.wb*q.y - qw.wc*q.z)/2.0;
  qs.xdot = (  qw.wa*q.r - qw.wb*q.z + qw.wc*q.y)/2.0;
  qs.ydot = (  qw.wa*q.z + qw.wb*q.r - qw.wc*q.x)/2.0;
  qs.zdot = (- qw.wa*q.y + qw.wb*q.x + qw.wc*q.r)/2.0;

  return(qs);
}


Matrix Rz(double angle)
{
  Matrix M;
  double ca, sa;

  ca = cos(angle);
  sa = sin(angle);

  M.xx = ca;
  M.xy = -sa;
  M.xz = 0.0;
  M.yx = sa;
  M.yy = ca;
  M.yz = 0.0;
  M.zx = 0.0;
  M.zy = 0.0;
  M.zz = 1.0;

  return(M);
}

Matrix Rx(double angle)
{
  Matrix M;
  double ca, sa;

  ca = cos(angle);
  sa = sin(angle);

  M.xx = 1.0;
  M.xy = 0.0;
  M.xz = 0.0;
  M.yx = 0.0;
  M.yy = ca;
  M.yz = -sa;
  M.zx = 0.0;
  M.zy = sa;
  M.zz = ca;

  return(M);
}

Matrix Ry(double angle)
{
  Matrix M;
  double ca, sa;

  ca = cos(angle);
  sa = sin(angle);

  M.xx = ca;
  M.xy = 0.0;
  M.xz = sa;
  M.yx = 0.0;
  M.yy = 1.0;
  M.yz = 0.0;
  M.zx = -sa;
  M.zy = 0.0;
  M.zz = ca;

  return(M);
}

