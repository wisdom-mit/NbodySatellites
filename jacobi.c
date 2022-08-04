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


void dolinks()
{
  int l, i, j;
  State relstate;

  for(i=0; i<nbodies; i++) {
    pj[i].x = p[i].x;
    pj[i].y = p[i].y;
    pj[i].z = p[i].z;
    pj[i].xd = p[i].xd;
    pj[i].yd = p[i].yd;
    pj[i].zd = p[i].zd;
  }

  for(l=0; l<nlinks; l++) {
    i = link[l].i;
    j = link[l].j;

    relstate.x = pj[j].x - pj[i].x;
    relstate.y = pj[j].y - pj[i].y;
    relstate.z = pj[j].z - pj[i].z;
    relstate.xd = pj[j].xd - pj[i].xd;
    relstate.yd = pj[j].yd - pj[i].yd;
    relstate.zd = pj[j].zd - pj[i].zd;

    /* cmstate */
    pj[j].x = fi[l]*pj[i].x + fj[l]*pj[j].x;
    pj[j].y = fi[l]*pj[i].y + fj[l]*pj[j].y;
    pj[j].z = fi[l]*pj[i].z + fj[l]*pj[j].z;
    pj[j].xd = fi[l]*pj[i].xd + fj[l]*pj[j].xd;
    pj[j].yd = fi[l]*pj[i].yd + fj[l]*pj[j].yd;
    pj[j].zd = fi[l]*pj[i].zd + fj[l]*pj[j].zd;

    pj[i].x = relstate.x; pj[i].y = relstate.y; pj[i].z = relstate.z;
    pj[i].xd = relstate.xd; pj[i].yd = relstate.yd; pj[i].zd = relstate.zd;
  }
}

void dolinks_acc_only()
{
  int l, i, j;
  State relstate;

  for(i=0; i<nbodies; i++) {
    pj[i].xdd = p[i].xdd;
    pj[i].ydd = p[i].ydd;
    pj[i].zdd = p[i].zdd;
  }

  for(l=0; l<nlinks; l++) {
    i = link[l].i;
    j = link[l].j;

    relstate.xdd = pj[j].xdd - pj[i].xdd;
    relstate.ydd = pj[j].ydd - pj[i].ydd;
    relstate.zdd = pj[j].zdd - pj[i].zdd;

    /* cmstate */
    pj[j].xdd = fi[l]*pj[i].xdd + fj[l]*pj[j].xdd;
    pj[j].ydd = fi[l]*pj[i].ydd + fj[l]*pj[j].ydd;
    pj[j].zdd = fi[l]*pj[i].zdd + fj[l]*pj[j].zdd;

    pj[i].xdd = relstate.xdd; 
    pj[i].ydd = relstate.ydd; 
    pj[i].zdd = relstate.zdd;
  }
}

void undolinks()
{
  int l, i, j;
  State statei;

  for(i=0; i<nbodies; i++) {
    p[i].x = pj[i].x;
    p[i].y = pj[i].y;
    p[i].z = pj[i].z;
    p[i].xd = pj[i].xd;
    p[i].yd = pj[i].yd;
    p[i].zd = pj[i].zd;
  }

  for(l=0; l<nlinks; l++) {
    i = link[nlinks-l-1].i;
    j = link[nlinks-l-1].j;

    statei.x = p[j].x - factor2[l]*p[i].x;
    statei.y = p[j].y - factor2[l]*p[i].y;
    statei.z = p[j].z - factor2[l]*p[i].z;
    statei.xd = p[j].xd - factor2[l]*p[i].xd;
    statei.yd = p[j].yd - factor2[l]*p[i].yd;
    statei.zd = p[j].zd - factor2[l]*p[i].zd;

    p[j].x = p[j].x + factor1[l]*p[i].x;
    p[j].y = p[j].y + factor1[l]*p[i].y;
    p[j].z = p[j].z + factor1[l]*p[i].z;
    p[j].xd = p[j].xd + factor1[l]*p[i].xd;
    p[j].yd = p[j].yd + factor1[l]*p[i].yd;
    p[j].zd = p[j].zd + factor1[l]*p[i].zd;

    p[i].x = statei.x; p[i].y = statei.y; p[i].z = statei.z;
    p[i].xd = statei.xd; p[i].yd = statei.yd; p[i].zd = statei.zd;
  }
}


void domasses()
{
  int l, i, j, n;
  double totalmass, reducedmass;

  for(n=0; n<nbodies; n++) {
    GMj[n] = GM[n];
  }

  for(l=0; l<nlinks; l++) {
    i = link[l].i;
    j = link[l].j;
    /* printf("%d %d\n", i, j); */

    totalmass = GMj[i] + GMj[j];
    reducedmass = GMj[i]*GMj[j]/totalmass;
    fi[l] = GMj[i]/totalmass;
    fj[l] = GMj[j]/totalmass;

    GMjprevious[i] = GMj[i];
    GMj[i] = reducedmass;
    GMjtotal[i] = totalmass;
    GMj[j] = totalmass;
    kc[i] = totalmass;
  }
}

void undomasses()
{
  int l, n, i, j;
  double totalmass, m1, m2;

  for(n=0; n<nbodies; n++) {
    GM[n] = GMj[n];
  }


  for(l=0; l<nlinks; l++) {
    i = link[nlinks-l-1].i;
    j = link[nlinks-l-1].j;

    m1 = GMjprevious[i];
    totalmass = GM[j];
    m2 = totalmass - m1;

    factor1[l] = m1/totalmass;
    factor2[l] = m2/totalmass;

    GM[i] = m1;
    GM[j] = m2;
  }
}

void compute_masses()
{
  void domasses(), undomasses();

  domasses();
  undomasses();
}
