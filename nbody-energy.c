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


double energy()
{
  double dx, dy, dz, rij2;
  double x, y, z, r, P2, P3, P4, P5, P6, gamma, gamma2;
  double wa, wb, wc;
  int i, j;
  Vector pole;
  Vector zhat;
  double e0, e1, e2, e3;

  e0 = 0.0;
  e1 = 0.0;
  e2 = 0.0;
  e3 = 0.0;

  wa = rr.wa;
  wb = rr.wb;
  wc = rr.wc;

  e0 += 0.5*(AA*wa*wa + BB*wb*wb + CC*wc*wc)*G;

  zhat.x = 0.;
  zhat.y = 0.;
  zhat.z = 1.;
  pole = quaternion_rotate_vector(rr.q, zhat);

  for(i=0; i<nbodies; i++) {
    e1 += 0.5*GM[i]*
      (p[i].xd*p[i].xd + p[i].yd*p[i].yd + p[i].zd*p[i].zd);
  }

  for(i=0; i<nbodies-1; i++) {
    for(j=i+1; j<nbodies; j++) {
      dx = p[i].x - p[j].x;
      dy = p[i].y - p[j].y;
      dz = p[i].z - p[j].z;
      rij2 = dx*dx + dy*dy + dz*dz;
      e1 -= GM[i]*GM[j]/sqrt(rij2);
    }
  }

  for(i=0; i<nsatellites; i++) {
    x = p[i].x - p[PLANET].x;
    y = p[i].y - p[PLANET].y;
    z = p[i].z - p[PLANET].z;
    r = sqrt(x*x + y*y  +z*z);
    gamma = (pole.x*x + pole.y*y + pole.z*z)/r;
    gamma2 = gamma*gamma;

    P2 = 1.5*gamma2 - 0.5;
    P3 = (2.5*gamma2 - 1.5)*gamma;
    /* (n+1)P_{n+1} = (2n+1) x P_n - n P_{n-1} */
    P4 = 1.75*gamma*P3 - 0.75*P2;
    P5 = 1.8*gamma*P4 - 0.8*P3;
    P6 = (11.0*gamma*P5 - 5.0*P4)/6.0;

    e2 += GM[i]*GM[PLANET]*(Re*Re/(r*r*r))*(J2*P2 + J3*P3*(Re/r) + J4*P4*(Re/r)*(Re/r) + J5*P5*(Re/r)*(Re/r)*(Re/r) + J6*P6*(Re/r)*(Re/r)*(Re/r)*(Re/r));
  }

  x = p[SUN].x - p[PLANET].x;
  y = p[SUN].y - p[PLANET].y;
  z = p[SUN].z - p[PLANET].z;
  r = sqrt(x*x + y*y  +z*z);
  gamma = (pole.x*x + pole.y*y + pole.z*z)/r;
  P2 = 1.5*gamma*gamma - 0.5;
  e3 += GM[SUN]*GM[PLANET]*((Re*Re/(r*r*r))*J2*P2);

  /* printf("%.16le %.16le %.16le %.16le\n", e0, e1, e2, e3); */

  return(e0 + e1 + e2 + e3);
}

double E()
{
  double e;
  void dolinks(), undolinks();

  undolinks();
  e = energy();
  dolinks();
  
  return(e);
}
