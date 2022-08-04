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
#include "rigid.c"
#include "nbody-map.c"
#include "nbody-io.c"
#include "jacobi.c"
#include "nbody-energy.c"
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
  }
  printf("----------------\n");
  for(j=0; j<nbodies; j++) {
    printf("%d\n", j);
    printf("%.16le %.16le %.16le\n", pj[j].x , pj[j].y , pj[j].z);
    printf("%.16le %.16le %.16le\n", pj[j].xd , pj[j].yd , pj[j].zd);
  }
  printf("****************************************************************\n");
}

int main(int argc, char** argv)
{
  int i, step, n_per_output, n_outputs, step0;
  void print_states();
  void compute_corrector_coefficients();
  void domasses(), undomasses(), compute_weights();
  int write_body_states();
  void B(double);
  void A(double);
  double energy();
  char filename[100], icfilename[100], output_filename[100];
  FILE *file, *outputfile;
  
  if(argc != 2) {
    printf("nbody setupfile\n");
    exit(-1);
  }

  if(sscanf(argv[1], "%s", filename) == 0) {
    exit(1);
  }
  file = fopen(filename, "r");
  if(file == NULL){
    printf("setupfile does not exist: %s\n", filename);
    exit(2);
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
  if(outputfile != NULL){
    step0 = fast_last(outputfile, 1, nbodies*sizeof(BinaryData) + sizeof(QWState));
    printf("step0: %d\n", step0);
    fclose(outputfile);
  } else {
    step0 = 0;
  }
  printf("STEP0: %d\n", step0);
  if(step0<=0) {
    outputfile = fopen(output_filename, "w");
  } else {
    outputfile = fopen(output_filename, "r+");
  }

  compute_masses();
  compute_corrector_coefficients(dt);
  dolinks();
  compute_tidal_constants();

  print_states();

  if(step0 <= 0) {
    step0 = 0;
    time = 0.0;
    dolinks();
    print_states();
    real_to_map();
    A(dt/2.0);
    write_body_states(outputfile, 0);
  } else {
    double dt_local;
    dt_local = read_body_states(outputfile, step0);
    printf("dt_local: %.16le\n", dt_local);
  }

  for(step=step0+1; step<n_outputs; step++) {

    for(i=0; i<n_per_output; i++) {
      B(dt);
      A(dt);
    }
    write_body_states(outputfile, step);
  }
}

