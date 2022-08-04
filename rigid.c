/*

Copyright (C) 2022 Jack Wisdom

This file is part of nbody-satellites.  nbody-satellites is software supporting the work
of Jack Wisdom .

nbody-satellite is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

nbody-satellites is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with nbody-satellites.  If not, see <https://www.gnu.org/licenses/>.

*/

QWState axisymmetric_rigid(double AA, double CC, double dt, QWState *rr)
{
  Quaternion q, r1, r2, qp;
  QWState rrp;
  double wa, wb, wc;
  double a, b, wm, c, s, ca, sa;

  q = rr->q;
  wa = rr->wa;
  wb = rr->wb;
  wc = rr->wc;

  a = wc * (1.0 - CC/AA) * dt;
  b = wc * (CC/AA);

  r1.r = cos(a/2.);
  r1.x = 0.0;
  r1.y = 0.0;
  r1.z = sin(a/2.0);

  wm = sqrt(wa*wa + wb*wb + b*b);
  c = wm * dt;

  s = sin(c/2.0)/wm;
  r2.r = cos(c/2.0);
  r2.x = s * wa;
  r2.y = s * wb;
  r2.z = s * b;

  rrp.t = rr->t;

  ca = cos(a);
  sa = sin(a);
  rrp.wa = ca * wa + sa * wb;
  rrp.wb = ca * wb - sa * wa;
  rrp.wc = wc;

  qp = quaternion_multiply(quaternion_multiply(q, r2), r1);
  rrp.q = qp;

  return(rrp);
}
