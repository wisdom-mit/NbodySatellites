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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846
#define TWOPI 6.283185307179586476925287

#ifndef _STATE_
#define _STATE_
typedef struct {
  double x, y, z, xd, yd, zd, xdd, ydd, zdd;
} State;
#endif

typedef struct {
  double x, y, z;
} Vector;

typedef struct {
  int body, linki, linkj;
  double time, dt;
  State state;
  Vector w;
} BinaryData;

typedef struct {
  int i, j;
  double weight;
} Link;

#define MAX_N_BODIES 20
int nbodies, nsatellites;

#define MAX_N_PARTICLES 20

double G;
double GMsun, GM[MAX_N_BODIES], GMj[MAX_N_BODIES];
double GMjprevious[MAX_N_BODIES], GMjtotal[MAX_N_BODIES];
State p[MAX_N_BODIES], pj[MAX_N_BODIES];
double fi[MAX_N_BODIES], fj[MAX_N_BODIES];
double kc[MAX_N_BODIES];
double factor1[MAX_N_BODIES], factor2[MAX_N_BODIES];
double time;
double dt;
int PLANET, SUN, CM;
double J2, J3, J4, J5, J6, Re;

Vector pole;
double dadt_Re[MAX_N_BODIES], dedt_over_e[MAX_N_BODIES], didt_over_i[MAX_N_BODIES];
double tide_alpha0[MAX_N_BODIES], tide_alpha1[MAX_N_BODIES], tide_beta0[MAX_N_BODIES], tide_gamma[MAX_N_BODIES];

/* parameters for planet rotation */
double AA, BB, CC, C_minus_B_over_A, A_minus_C_over_B, B_minus_A_over_C;
QWState rr;
double planetary_angular_momentum, omega0;

Link link[MAX_N_BODIES*MAX_N_BODIES];
int nlinks;




