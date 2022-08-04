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

int write_body_states(FILE *file, int step)
{
  void copy_state();
  BinaryData bd;
  int body;

  fseek(file, step*(nbodies*sizeof(BinaryData) + sizeof(QWState)), 0);
  
  if(fwrite((char *)&rr, sizeof(QWState), 1, file) == 0) return(0);
  for(body=0; body<nbodies; body++) {
    bd.body = body;
    bd.linki = link[body].i;
    bd.linkj = link[body].j;
    bd.time = time;
    bd.dt = dt;
    copy_state(&pj[body], &bd.state);
    if(fwrite((char *)&bd, sizeof(BinaryData), 1, file) == 0) return(0);
  }
  fflush(file);
  return(1);
}

void copy_state(State *s1, State *s2)
{
  s2->x = s1->x;
  s2->y = s1->y;
  s2->z = s1->z;
  s2->xd = s1->xd;
  s2->yd = s1->yd;
  s2->zd = s1->zd;
  s2->xdd = s1->xdd;
  s2->ydd = s1->ydd;
  s2->zdd = s1->zdd;

  return;
}

double read_body_states(FILE *file, int step)
{
  BinaryData bd;
  int body;
  double dt_local, time_local;

  time_local = -666.0;
  dt_local = -666.0;

  fseek(file, step*(nbodies*sizeof(BinaryData) + sizeof(QWState)), 0);

  if(fread((char *)&rr, sizeof(QWState), 1, file) == EOF) {
    printf("oops\n");
  }
  for(body=0; body<nbodies; body++) {
    if(fread((char *)&bd, sizeof(BinaryData), 1, file) == EOF) {
      printf("oops %d\n", body);
    }
    copy_state(&bd.state, &pj[body]);

    if(body != bd.body) {
      printf("body states mismatch: %d %d\n", bd.body, body);
      exit(-1);
    }
    link[body].i = bd.linki;
    link[body].j = bd.linkj;

    if(body == 0) {
      time_local = bd.time;
      dt_local = bd.dt;
    } else {
      if(time_local != bd.time) {
	printf("time mismatch: %.16le %.16le %d\n", time_local, bd.time, body);
	exit(17);
      }
      if(dt_local != bd.dt) {
	printf("time step mismatch: %.16le %.16le %d\n", dt_local, bd.dt, body);
	exit(18);
      }
    }

    time = time_local;
  }
  return(dt_local);
}

int fast_last(FILE *file, int n_objects, int size_object)
{
  int last;

  fseek(file, 0, SEEK_END);
  last = ftell(file);
  if(last != -1) last /= (n_objects*size_object);
  return(last-1);
}

void read_ics(char* filename)
{
  FILE *file;
  int i, lasti;
  double mass;
  double km_AU, Re_km, polar_moment_MRe2, omega_rad_sec, sec_day;
  double AA_MRe2, BB_MRe2, CC_MRe2;
  double omega;
  double amag, angle;
  Vector axis;
  int version;

  file = fopen(filename, "r");
  if(file == NULL){
    printf("ic file does not exist: %s\n", filename);
    exit(-1);
  }

  if(fscanf(file, "%d", &version) == 0) {
    exit(7);
  }
  if(version != 2) {
    printf("wrong ic file version %d\n", version);
    exit(-1);
  }

  nbodies = 0;
  nsatellites = 0;
  if(fscanf(file, "%lf", &G) == 0) {
    exit(8);
  }
  if(fscanf(file, "%d", &PLANET) == 0) {
    exit(9);
  }
  if(fscanf(file, "%lf %lf %lf", &pole.x, &pole.y, &pole.z) == 0) {
    exit(10);
  }
  if(fscanf(file, "%lf %lf", &polar_moment_MRe2, &omega_rad_sec) == 0) {
    exit(11);
  }
  if(fscanf(file, "%lf %lf %lf %lf %lf %lf", &J2, &J3, &J4, &J5, &J6, &Re_km) == 0) {
    exit(12);
  }

  while(fscanf(file, "%lf", &mass) != EOF) {
    if(fscanf(file, "%lf", &(p[nbodies].x)) == 0) {
      exit(13);
    }
    if(fscanf(file, "%lf", &(p[nbodies].y)) == 0) {
      exit(14);
    }
    if(fscanf(file, "%lf", &(p[nbodies].z)) == 0) {
      exit(15);
    }
    if(fscanf(file, "%lf", &(p[nbodies].xd)) == 0) {
      exit(16);
    }
    if(fscanf(file, "%lf", &(p[nbodies].yd)) == 0) {
      exit(17);
    }
    if(fscanf(file, "%lf", &(p[nbodies].zd)) == 0) {
      exit(18);
    }

    p[nbodies].xdd = 0.0;
    p[nbodies].ydd = 0.0;
    p[nbodies].zdd = 0.0;

    if(mass<0.0) {
      mass = - mass;
      GM[nbodies] = G * mass;
      if(fscanf(file, "%lf %lf %lf", &dadt_Re[nbodies], &dedt_over_e[nbodies], &didt_over_i[nbodies]) == 0) {
	exit(19);
      }
      nsatellites++;
    } else {
      GM[nbodies] = G * mass;
    }
    nbodies++;
  }

  for(i=0; i<nsatellites; i++) {
    p[i].x += p[PLANET].x;
    p[i].y += p[PLANET].y;
    p[i].z += p[PLANET].z;
    p[i].xd += p[PLANET].xd;
    p[i].yd += p[PLANET].yd;
    p[i].zd += p[PLANET].zd;    
    p[i].xdd = 0.0;
    p[i].ydd = 0.0;
    p[i].zdd = 0.0;
  }

  /* adjust constants */
  km_AU = 1.495978707e8; /* was   km_AU = 1.49597870e8; */
  sec_day = 86400.0;
  Re = Re_km/km_AU;

  /* these are PLANET parameters */
  CC_MRe2 = polar_moment_MRe2;
  AA_MRe2 = polar_moment_MRe2 - J2;
  BB_MRe2 = AA_MRe2;
  C_minus_B_over_A = (CC_MRe2 - BB_MRe2)/AA_MRe2;
  A_minus_C_over_B = (AA_MRe2 - CC_MRe2)/BB_MRe2;
  B_minus_A_over_C = (BB_MRe2 - AA_MRe2)/CC_MRe2;
  AA = AA_MRe2*Re*Re*GM[PLANET]/G;
  BB = BB_MRe2*Re*Re*GM[PLANET]/G;
  CC = CC_MRe2*Re*Re*GM[PLANET]/G;

  /* set PLANET rotation state */
  omega = omega_rad_sec*sec_day;
  /* the axis is proportional to zhat cross pole */
  axis.x = -pole.y;
  axis.y =  pole.x;
  axis.z = 0.0;
  amag = sqrt(axis.x*axis.x + axis.y*axis.y + axis.z*axis.z);
  axis.x /= amag;
  axis.y /= amag;
  axis.z /= amag;
  angle = acos(pole.z);
  rr.q.r = cos(angle/2.0);
  rr.q.x = sin(angle/2.0)*axis.x;
  rr.q.y = sin(angle/2.0)*axis.y;
  rr.q.z = sin(angle/2.0)*axis.z; 
  rr.wa = 0.0;
  rr.wb = 0.0;
  rr.wc = omega;

  omega0 = omega;
  planetary_angular_momentum = polar_moment_MRe2*omega0;

  lasti = PLANET;
  SUN = nsatellites;

  for(i=0; i<nsatellites; i++) {
    link[i].i = lasti;
    link[i].j = i;
    lasti = i;
  }
  for(i=nsatellites; i<nbodies-1; i++) {
    link[i].i = i;
    link[i].j = i+1;
  }
  if(nsatellites > 0) {
    link[PLANET-1].j = nsatellites-1;
    link[PLANET].i = nsatellites-1;
  }
  nlinks = nbodies - 1;

  /*
    printf("links   PLANET=%d\n", PLANET);
    for(i=0; i<nlinks; i++) {
    printf("%d %d %d\n", i, link[i].i, link[i].j);
    }
   */

  CM = link[nlinks-1].j;
  if(CM != (nbodies-1)) {
    printf("CM = %d   nbodies = %d\n", CM, nbodies);
    exit(-1);
  }    
}

void compute_tidal_constants()
{
  int i, ij;
  double x, y, z, vx, vy, vz;
  double r, vs, a, factor;
  double dadt;

  for(i=0; i<nsatellites; i++) {
    ij = link[i].i;
    dadt = dadt_Re[i]*Re;
    x = pj[ij].x;
    y = pj[ij].y;
    z = pj[ij].z;
    vx = pj[ij].xd;
    vy = pj[ij].yd;
    vz = pj[ij].zd;
    r = sqrt(x*x + y*y + z*z);
    vs = vx*vx + vy*vy + vz*vz;
    a = 1./(2./r - vs/kc[ij]);
    /* printf("a/Re[%d] %.16le\n", i, a/Re); */
    factor = (2.0*a*a*(2./r - 1./a));
    tide_alpha0[i] = dadt/factor;
    /* with Kaula term */
    tide_alpha1[i] = 2.0*a*(dedt_over_e[i] - dadt/(4.0*a))/factor;
    /* tide_beta0[i] = 2.0*(dedt_over_e[i] - dadt/(4.0*a) + 2.*tide_alpha0[i]);  why is the last term there? */
    tide_beta0[i] = 2.0*(dedt_over_e[i] - dadt/(4.0*a));
    tide_gamma[i] = 2.0*(didt_over_i[i] - dadt/(4.0*a));
    
    /* without Kaula term 
    tide_alpha1[i] = 2.0*a*dedt_over_e[i]/factor;
    tide_beta0[i] = 2.0*dedt_over_e[i];
    tide_gamma[i] = 2.0*didt_over_i[i];
     */

  }
}

void subtract_cm()
{
  int i;
  State XX;
  double M;

  XX.x = 0.0;  XX.y = 0.0;  XX.z = 0.0;
  XX.xd = 0.0;  XX.yd = 0.0;  XX.zd = 0.0;
  M = 0;

  for(i=0; i<nbodies; i++) {
    XX.x += GM[i]*p[i].x;
    XX.y += GM[i]*p[i].y;
    XX.z += GM[i]*p[i].z;
    XX.xd += GM[i]*p[i].xd;
    XX.yd += GM[i]*p[i].yd;
    XX.zd += GM[i]*p[i].zd;
    M += GM[i];
  }

  XX.x /= M;
  XX.y /= M;
  XX.z /= M;
  XX.xd /= M;
  XX.yd /= M;
  XX.zd /= M;

  for(i=0; i<nbodies; i++) {
    p[i].x -= XX.x;
    p[i].y -= XX.y;
    p[i].z -= XX.z;

    p[i].xd -= XX.xd;
    p[i].yd -= XX.yd;
    p[i].zd -= XX.zd;
  }
}





