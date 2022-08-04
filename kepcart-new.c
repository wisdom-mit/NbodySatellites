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

void rotate_z(double theta, State state, State *rstate)
{
  double c, s;

  c = cos(theta);
  s = sin(theta);


  rstate->x = state.x*c - state.y*s;
  rstate->y = state.x*s + state.y*c;
  rstate->z = state.z;


  rstate->xd = state.xd*c - state.yd*s;
  rstate->yd = state.xd*s + state.yd*c;
  rstate->zd = state.zd;
}

void rotate_x(double theta, State state, State *rstate)
{
  double c, s;

  c = cos(theta);
  s = sin(theta);

  rstate->x = state.x;
  rstate->y = state.y*c - state.z*s;
  rstate->z = state.y*s + state.z*c;

  rstate->xd = state.xd;
  rstate->yd = state.yd*c - state.zd*s;
  rstate->zd = state.yd*s + state.zd*c;
}


void keplerian(double gm, State state, 
	  double *a, double *e, double *i, double *longnode, double *argperi, double *meananom)
{
  double rxv_x, rxv_y, rxv_z, hs, h, parameter;
  double r, vs, rdotv, rdot, ecostrueanom, esintrueanom, cosnode, sinnode;
  double rcosu, rsinu, u, trueanom, eccanom;
  State rstate, rstate2;
  double cosi, delta_i;

  /* find direction of angular momentum vector */
  rxv_x = state.y * state.zd - state.z * state.yd;
  rxv_y = state.z * state.xd - state.x * state.zd;
  rxv_z = state.x * state.yd - state.y * state.xd;
  hs = rxv_x * rxv_x + rxv_y * rxv_y + rxv_z * rxv_z;
  h = sqrt(hs);

  r = sqrt(state.x * state.x + state.y * state.y + state.z * state.z);
  vs = state.xd * state.xd + state.yd * state.yd + state.zd * state.zd;
  /* v = sqrt(vs);  unnecessary */
  rdotv = state.x * state.xd + state.y * state.yd + state.z * state.zd;
  rdot = rdotv / r;
  parameter = hs / gm;

  cosi = rxv_z / h;
  *i = acos(cosi);

  if(rxv_x!=0.0 || rxv_y!=0.0) {
    *longnode = atan2(rxv_x, -rxv_y);
  } else {
    *longnode = 0.0;
  }
  *longnode -= (2.0*M_PI)*floor(*longnode/(2.0*M_PI));

  *a = 1.0 / (2.0 / r - vs / gm);

  ecostrueanom = parameter / r - 1.0;
  esintrueanom = rdot * h / gm;
  *e = sqrt(ecostrueanom * ecostrueanom + esintrueanom * esintrueanom);

  if(esintrueanom!=0.0 || ecostrueanom!=0.0) {
    trueanom = atan2(esintrueanom, ecostrueanom);
  } else {
    trueanom = 0.0;
  }

  cosnode = cos(*longnode);
  sinnode = sin(*longnode);

  if(fabs(cosi) > 1.e-3) {
    /* u is the argument of latitude */
    rcosu = state.x * cosnode + state.y * sinnode;
    rsinu = (state.y * cosnode - state.x * sinnode)/cos(*i);

  } else {
    delta_i = 1.0;
    rotate_z(- *longnode, state, &rstate);
    rotate_x(delta_i, rstate, &rstate2);
    rotate_z(*longnode, rstate2, &state);
    
    rcosu = state.x * cosnode + state.y * sinnode;
    rsinu = (state.y * cosnode - state.x * sinnode)/cos(*i + delta_i);
  }

  if(rsinu!=0.0 || rcosu!=0.0) {
    u = atan2(rsinu, rcosu);
  } else {
    u = 0.0;
  }

  *argperi = u - trueanom;
  *argperi -= (2.0*M_PI)*floor(*argperi/(2.0*M_PI));
  

  eccanom = 2.0 * atan(sqrt((1.0 - *e)/(1.0 + *e)) * tan(trueanom/2.0));
  *meananom = eccanom - *e * sin(eccanom);

  return;
}

/*
double Kepler_equation(double e, double M)
{
  double E0, E1, E2, den;
  double machine_epsilon=1.6e-16;
  int i;

  E0 = M; 
  i = 0;
  do {
    if(i++ > 1000) {
      printf("too many iters\n");
      exit(-1);
    }
    E1 = M + e * sin(E0);
    E2 = M + e * sin(E1);
    if(fabs(E0-E2) < 10.0*machine_epsilon) return(E2);

    den = E2 - 2.0*E1 + E0;
    if(fabs(den) > machine_epsilon) {
      E0 = E0 - (E1-E0)*(E1-E0)/den;
    }
    else {
      E0 = E2;
      E2 = E1;
    }
  } while(fabs(E0-E2) > 2.0*machine_epsilon);
  return(E0);
}
*/


double Kepler_equation(double e, double M)
{
  double E, err, nerr;
  double esinE, ecosE;
  double f, fp, fpp, fppp;
  double deltaN, deltaH, deltaD;

  nerr = 1.0;
  E = M;
  E = M + e*sin(E);
  E = M + e*sin(E);
  do {
    err = nerr;

    esinE = e*sin(E);
    ecosE = e*cos(E);

    f = E - esinE - M;
    fp = 1 - ecosE;
    fpp = esinE;
    fppp = ecosE;

    deltaN = -f/fp;
    deltaH = -f/(fp + deltaN*fpp/2.0);
    deltaD = -f/(fp + deltaH*(fpp + deltaH*fppp/3.0)/2.0);

    nerr = fabs(deltaD);
    E += deltaD;

  } while(err > 1.e-6);
    
  return(E);
}


void cartesian(double gm, 
	       double a, double e, double i, double longnode, double argperi, double meananom, 
	       State *state)
{
  double meanmotion, cosE, sinE, foo;
  double x, y, z, xd, yd, zd;
  double xp, yp, zp, xdp, ydp, zdp;
  double cosw, sinw, cosi, sini, cosnode, sinnode;
  double E0;

  E0 = Kepler_equation(e, meananom);

  cosE = cos(E0);
  sinE = sin(E0);

  /* compute unrotated positions and velocities */
  foo = sqrt(1.0 - e*e);
  meanmotion = sqrt(gm/(a*a*a));
  x = a * (cosE - e);
  y = foo * a * sinE;
  z = 0.0;
  xd = -a * meanmotion * sinE / (1.0 - e * cosE);
  yd = foo * a * meanmotion * cosE / (1.0 - e * cosE);
  zd = 0.0;

  /* rotate by argument of perihelion in orbit plane*/
  cosw = cos(argperi);
  sinw = sin(argperi);
  xp = x * cosw - y * sinw;
  yp = x * sinw + y * cosw;
  zp = z;
  xdp = xd * cosw - yd * sinw;
  ydp = xd * sinw + yd * cosw;
  zdp = zd;

  /* rotate by inclination about x axis */
  cosi = cos(i);
  sini = sin(i);
  x = xp;
  y = yp * cosi - zp * sini;
  z = yp * sini + zp * cosi;
  xd = xdp;
  yd = ydp * cosi - zdp * sini;
  zd = ydp * sini + zdp * cosi;

  /* rotate by longitude of node about z axis */
  cosnode = cos(longnode);
  sinnode = sin(longnode);
  state->x = x * cosnode - y * sinnode;
  state->y = x * sinnode + y * cosnode;
  state->z = z;
  state->xd = xd * cosnode - yd * sinnode;
  state->yd = xd * sinnode + yd * cosnode;
  state->zd = zd;

  return;
}

