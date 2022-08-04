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

/*
void leap(double dt, int n)
{
  void B(double);
  void A(double);
  int i;
  
  A(dt/2.0);
  
  for(i=0; i<n-1; i++) {
    B(dt);
    A(dt);
  }
  B(dt);
  A(dt/2);
}
*/

void A(double dt)
{
  int i;
  void kepler_step();
  QWState rrp;

  for(i=0; i<nbodies-1; i++) {

    State ans;

    kepler_step(kc[i], dt, &(pj[i]), &ans);

    pj[i].x  = ans.x;      pj[i].y  = ans.y;      pj[i].z  = ans.z;
    pj[i].xd = ans.xd;     pj[i].yd = ans.yd;     pj[i].zd = ans.zd;
    
  }
  rrp = axisymmetric_rigid(AA, CC, dt, &rr);
  rr.q = rrp.q;
  rr.wa = rrp.wa; rr.wb = rrp.wb; rr.wc = rrp.wc;

  time += dt;
}

double compute_a(double gm, State *s0)
{
  double r0, v0s, a;

  r0 = sqrt(s0->x*s0->x + s0->y*s0->y + s0->z*s0->z);
  v0s = s0->xd*s0->xd + s0->yd*s0->yd + s0->zd*s0->zd;
  a = 1.0/(2.0/r0 - v0s/gm);

  return(a);
}

double compute_n(double gm, State *s0)
{
  double r0, v0s, a, n;

  r0 = sqrt(s0->x*s0->x + s0->y*s0->y + s0->z*s0->z);
  v0s = s0->xd*s0->xd + s0->yd*s0->yd + s0->zd*s0->zd;
  a = 1.0/(2.0/r0 - v0s/gm);
  n = sqrt(gm/(a*a*a));

  return(n);
}

void B(double dt)
{
  void dolinks(), dolinks_acc_only(), undolinks();
  void accelerations(), kicks(double);
  void tides(double);

  undolinks();
  accelerations(dt);
  dolinks_acc_only();
  kicks(dt);
  tides(dt);
}

void kicks(double dt)
{
  double rp2, rp3;

  int i;

  for(i=0; i<nbodies-1; i++) {
    rp2 = pj[i].x*pj[i].x + pj[i].y*pj[i].y + pj[i].z*pj[i].z;
    rp3 = rp2*sqrt(rp2);
    pj[i].xd += dt*(pj[i].xdd + GMjtotal[i]*pj[i].x/rp3);
    pj[i].yd += dt*(pj[i].ydd + GMjtotal[i]*pj[i].y/rp3);
    pj[i].zd += dt*(pj[i].zdd + GMjtotal[i]*pj[i].z/rp3);
  }

}

void tides(double dt)
{
  int i, ij;
  double gm, vs, rdot, pvhat, parameter, e2, ta, tb, tg, r;
  Vector pos, vel, rxv, xyzhat, kspace;
  Vector zhat;
  Vector delta_L, delta_Lbody;

  delta_L.x = 0.0;
  delta_L.y = 0.0;
  delta_L.z = 0.0;

  zhat.x = 0.; zhat.y = 0.; zhat.z = 1.;
  pole = quaternion_rotate_vector(rr.q, zhat);

  for(i=0; i<nsatellites; i++) {

    ij = link[i].i;

    gm = kc[ij];
    pos.x = pj[ij].x;
    pos.y = pj[ij].y;
    pos.z = pj[ij].z;
    vel.x = pj[ij].xd;
    vel.y = pj[ij].yd;
    vel.z = pj[ij].zd;
    vs = vel.x*vel.x + vel.y*vel.y + vel.z*vel.z;
    r = sqrt(pos.x*pos.x + pos.y*pos.y + pos.z*pos.z);
    xyzhat.x = pos.x/r;
    xyzhat.y = pos.y/r;
    xyzhat.z = pos.z/r;
    rdot = xyzhat.x*vel.x + xyzhat.y*vel.y + xyzhat.z*vel.z;
    pvhat = (pole.x*vel.x + pole.y*vel.y + pole.z*vel.z)/sqrt(vs);
    rxv.x = pos.y*vel.z - pos.z*vel.y;
    rxv.y = pos.z*vel.x - pos.x*vel.z;
    rxv.z = pos.x*vel.y - pos.y*vel.x;
    parameter = (rxv.x*rxv.x + rxv.y*rxv.y + rxv.z*rxv.z)/gm;
    e2 = (parameter/r - 1.0)*(parameter/r - 1.0) + rdot*rdot*parameter/gm;
    /* printf("%d %.16le %.16le %.16le\n", i, e2, tide_alpha0[i], tide_alpha1[i]); */
    ta = tide_alpha0[i] + e2*tide_alpha1[i];
    tb = tide_beta0[i] + e2*4.0*tide_alpha1[i];
    tg = tide_gamma[i];
    kspace.x = (tb*rdot*(xyzhat.x - rdot*vel.x/vs) + ta*vel.x + tg*pvhat*(xyzhat.y*vel.z - xyzhat.z*vel.y));
    kspace.y = (tb*rdot*(xyzhat.y - rdot*vel.y/vs) + ta*vel.y + tg*pvhat*(xyzhat.z*vel.x - xyzhat.x*vel.z));
    kspace.z = (tb*rdot*(xyzhat.z - rdot*vel.z/vs) + ta*vel.z + tg*pvhat*(xyzhat.x*vel.y - xyzhat.y*vel.x));


    delta_L.x += dt*GMj[ij]*(pj[ij].y*kspace.z - pj[ij].z*kspace.y)/G;
    delta_L.y += dt*GMj[ij]*(pj[ij].z*kspace.x - pj[ij].x*kspace.z)/G;
    delta_L.z += dt*GMj[ij]*(pj[ij].x*kspace.y - pj[ij].y*kspace.x)/G;

    pj[ij].xd += dt*kspace.x; 
    pj[ij].yd += dt*kspace.y; 
    pj[ij].zd += dt*kspace.z; 
  }

  delta_Lbody = quaternion_rotate_vector(quaternion_conjugate(rr.q), delta_L);
  rr.wa -= delta_Lbody.x/AA;
  rr.wb -= delta_Lbody.y/BB;
  rr.wc -= delta_Lbody.z/CC;

}

void accelerations(double dt)
{
  double dx, dy, dz, tx, ty, tz, rij2, fij;
  int i, j;
  double Re2;
  Vector zhat;

  for(i=0; i<nbodies; i++) {
    p[i].xdd = 0.0; p[i].ydd = 0.0; p[i].zdd = 0.0;
  }

  for(i=0; i<nbodies-1; i++)
    for(j=i+1; j<nbodies; j++) {
      dx = p[i].x - p[j].x;
      dy = p[i].y - p[j].y;
      dz = p[i].z - p[j].z;
      rij2 = dx*dx + dy*dy + dz*dz;
      fij = 1.0/(rij2*sqrt(rij2));
      tx = dx*fij; ty = dy*fij; tz = dz*fij;
      p[i].xdd -= GM[j]*tx; p[i].ydd -= GM[j]*ty; p[i].zdd -= GM[j]*tz;
      p[j].xdd += GM[i]*tx; p[j].ydd += GM[i]*ty; p[j].zdd += GM[i]*tz;
    }
  
  zhat.x = 0.; zhat.y = 0.; zhat.z = 1.;
  pole = quaternion_rotate_vector(rr.q, zhat);
  Re2 = Re*Re;

  for(i=0; i<nsatellites; i++) {
    double alpha, beta, gamma, coef, factor;
    Vector acc;

    double x, y, z, r2, r, P2;
    Vector xyzhat;
    Vector abc;
    double P3, P4, P5, P6;
    double DP2, DP3, DP4, DP5, DP6;
    double rfactor, gamma2;
    double aa, bb;

    /* PLANET figure to satellites */
    x = p[i].x - p[PLANET].x;
    y = p[i].y - p[PLANET].y;
    z = p[i].z - p[PLANET].z;
    r2 = x*x + y*y + z*z;
    r = sqrt(r2);
    xyzhat.x = x/r;
    xyzhat.y = y/r;
    xyzhat.z = z/r;

    abc = quaternion_rotate_vector(quaternion_conjugate(rr.q), xyzhat);
    alpha = abc.x;
    beta = abc.y;
    gamma = abc.z;
    gamma2 = gamma*gamma;
    P2 = 1.5*gamma2 - 0.5;
    DP2 = 3.0*gamma;
    P3 = (2.5*gamma2 - 1.5)*gamma;
    DP3 = 7.5*gamma2 - 1.5;
    /* (n+1)P_{n+1} = (2n+1) x P_n - n P_{n-1} */
    P4 = 1.75*gamma*P3 - 0.75*P2;
    DP4 = 1.75*(P3 + gamma*DP3) - 0.75*DP2;
    P5 = 1.8*gamma*P4 - 0.8*P3;
    DP5 = 1.8*(P4 + gamma*DP4) - 0.8*DP3;
    P6 = (11.0*gamma*P5 - 5.0*P4)/6.0;
    DP6 = (11.0*(P5 + gamma*DP5) - 5.0*DP4)/6.0;
    rfactor = Re2/(r2*r2);
    aa = rfactor*(J2*DP2 + (J3*DP3 + (J4*DP4 + (J5*DP5 + J6*DP6*(Re/r))*(Re/r))*(Re/r))*(Re/r));
    bb = rfactor*(J2*(gamma*DP2 + 3.0*P2) +
		  (J3*(gamma*DP3 + 4.0*P3) +
		   (J4*(gamma*DP4 + 5.0*P4) +
		    (J5*(gamma*DP5 + 6.0*P5) +
		     J6*(gamma*DP6 + 7.0*P6)*(Re/r))*(Re/r))*(Re/r))*(Re/r));
    acc.x = bb*xyzhat.x - aa*pole.x;
    acc.y = bb*xyzhat.y - aa*pole.y;
    acc.z = bb*xyzhat.z - aa*pole.z;

    coef = (GM[PLANET]/G)*GM[i]*Re2/(r*r2);
    factor = coef*(J2*DP2 + J3*DP3*(Re/r) + J4*(Re2/r2)*DP4 + J6*(Re2/r2)*(Re2/r2)*DP6)*dt;
    rr.wa += factor*beta/AA;
    rr.wb -= factor*alpha/AA;
    /* wc unchanged */

    p[i].xdd += GM[PLANET]*acc.x;
    p[PLANET].xdd -= GM[i]*acc.x;
    p[i].ydd += GM[PLANET]*acc.y;
    p[PLANET].ydd -= GM[i]*acc.y;
    p[i].zdd += GM[PLANET]*acc.z;
    p[PLANET].zdd -= GM[i]*acc.z;
  }

  /* PLANET figure to SUN */
  {
    Vector abc, xyzhat, acc;
    double x, y, z, r, r2, aa, bb;
    double alpha, beta, gamma;
    double coef, factor, P2, DP2, rfactor;

    x = p[PLANET].x - p[SUN].x;
    y = p[PLANET].y - p[SUN].y;
    z = p[PLANET].z - p[SUN].z;

    r2 = x*x + y*y + z*z;
    r = sqrt(r2);
    xyzhat.x = x/r;
    xyzhat.y = y/r;
    xyzhat.z = z/r;
    abc = quaternion_rotate_vector(quaternion_conjugate(rr.q), xyzhat);
    alpha = abc.x;
    beta = abc.y;
    gamma = abc.z;

    /* kick from sun on planet figure */
    coef = 3.0*(GM[PLANET]/G)*GM[SUN]*Re2/(r*r2);
    factor = gamma*coef*J2*dt;
    rr.wa += factor*beta/AA;
    rr.wb -= factor*alpha/AA;

    rfactor = J2*Re2/(r2*r2);
    P2 = 1.5*gamma*gamma - 0.5;
    DP2 = 3.0*gamma;
    aa = rfactor*DP2;
    bb = rfactor*(gamma*DP2 + 3.0*P2);
    acc.x = bb*xyzhat.x - aa*pole.x;
    acc.y = bb*xyzhat.y - aa*pole.y;
    acc.z = bb*xyzhat.z - aa*pole.z;

    p[SUN].xdd -= GM[PLANET]*acc.x;
    p[PLANET].xdd += GM[SUN]*acc.x;
    p[SUN].ydd -= GM[PLANET]*acc.y;
    p[PLANET].ydd += GM[SUN]*acc.y;
    p[SUN].zdd -= GM[PLANET]*acc.z;
    p[PLANET].zdd += GM[SUN]*acc.z;
  }
}
