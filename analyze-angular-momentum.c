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

#include "quaternion.h"
#include "nbody.h"

#include "quaternion.c"
#include "vector-stuff.c"
#include "kepcart-new.c"
#include "rigid.c"
#include "nbody-map.c"
#include "nbody-io.c"
#include "jacobi.c"
#include "nbody-energy.c"
#include "nbody-angular-momentum.c"
#include "corrector.c"
#include "universal.c"

void print_states()
{
  int j;

  undolinks();

  printf("%.16le\n", time/365.25);
  printf("%.16le %.16le %.16le %.16le\n", rr.q.r, rr.q.x, rr.q.y, rr.q.z);
  printf("%.16le %.16le %.16le\n", rr.wa, rr.wb, rr.wc);
  for(j=0; j<nbodies; j++) {
    printf("%d\n", j);
    printf("%.16le %.16le %.16le\n", p[j].x , p[j].y , p[j].z);
    printf("%.16le %.16le %.16le\n", p[j].xd , p[j].yd , p[j].zd);
    printf("%.16le %.16le %.16le\n", p[j].xdd , p[j].ydd , p[j].zdd);
  }
  printf("----------------\n");
  for(j=0; j<nbodies; j++) {
    printf("%d\n", j);
    printf("%.16le %.16le %.16le\n", pj[j].x , pj[j].y , pj[j].z);
    printf("%.16le %.16le %.16le\n", pj[j].xd , pj[j].yd , pj[j].zd);
    printf("%.16le %.16le %.16le\n", pj[j].xdd , pj[j].ydd , pj[j].zdd);
  }
  printf("****************************************************************\n");
}

double pv0(double angle)
{
  angle = angle - TWOPI*floor(angle/TWOPI);
  if(angle > PI) angle -= TWOPI;
  return(angle);
}

double pv(double angle)
{
  angle = angle - TWOPI*floor(angle/TWOPI);
  return(angle);
}

int main(int argc, char** argv)
{
  int step, n_per_output, n_outputs, step0;
  void dolinks(), undolinks(), accelerations(), kicks(double), drifts(double);
  void compute_corrector_coefficients();
  void print_states();
  void domasses(), undomasses(), compute_weights();
  void read_ics(), print_links();
  int write_body_states();
  double energy();
  char filename[100], icfilename[100], output_filename[100];
  FILE *file, *outputfile;

  if(argc != 2) {
    printf("nbody setupfile\n");
    exit(-1);
  }

  sscanf(argv[1], "%s", filename);
  file = fopen(filename, "r");
  if(file == NULL){
    printf("setupfile does not exist: %s\n", filename);
    exit(-1);
  }
  if(fscanf(file, "%s", icfilename) == 0) {
    exit(3);
  }
  if(fscanf(file, "%lf", &dt) == 0) {
    exit(4);
  }
  if(fscanf(file, "%d", &n_per_output) == 0) {
    exit(5);
  }
  if(fscanf(file, "%d", &n_outputs) == 0) {
    exit(6);
  }

  read_ics(icfilename);
  subtract_cm();

  strcpy(output_filename, filename);
  strcat(output_filename, ".mb");
  outputfile = fopen(output_filename, "r");
  step0 = 0;
  if(outputfile != NULL){
    step0 = fast_last(outputfile, 1, nbodies*sizeof(BinaryData) + sizeof(QWState));
    fclose(outputfile);
  }
  if(step0<0) {
    printf("empty output file\n");
    exit(-1);
  }
  outputfile = fopen(output_filename, "r");

  compute_masses();
  compute_corrector_coefficients(dt);

  for(step=0; step<step0+1; step++) {

    Vector L, L0;
    Vector J, J0;
    double Lmag, Jmag;

    read_body_states(outputfile, step);

    A(-dt/2.0);
    map_to_real();
    
    undolinks();

    J0.x = 0.0; J0.y = 0.0; J0.z = 0.0;
    L0.x = 0.0; L0.y = 0.0; L0.z = 0.0;
    if(step == 0) {
      L0 = computeL();
      Lmag = sqrt(L0.x*L0.x + L0.y*L0.y + L0.z*L0.z);
      printf("; L: %.16le %.16le %.16le %.16le\n", L0.x, L0.y, L0.z, Lmag);
      J0 = compute_jacobi_L();
      Jmag = sqrt(J0.x*J0.x + J0.y*J0.y + J0.z*J0.z);
      printf("; J: %.16le %.16le %.16le %.16le\n", J0.x, J0.y, J0.z, Jmag);
      printf("(define %s", filename);
      printf("-L\n");
      printf(" '(\n");
    } else {
      L = computeL();
      J = compute_jacobi_L();
      printf("(%lf %.16le %.16le %.16le %.16le %.16le %.16le)\n", time/365.25, 
	     (L.x - L0.x)/Lmag, (L.y - L0.y)/Lmag, (L.z - L0.z)/Lmag, 
	     (J.x - J0.x)/Jmag, (J.y - J0.y)/Jmag, (J.z - J0.z)/Jmag);
    }
  }
  printf("))\n");
}
