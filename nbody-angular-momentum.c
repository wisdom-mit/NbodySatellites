Vector angmom()
{
  double wa, wb, wc;
  int i;

  Vector L;
  Vector Lbody;


  wa = rr.wa;
  wb = rr.wb;
  wc = rr.wc;

  Lbody.x = G*AA*wa;
  Lbody.y = G*BB*wb;
  Lbody.z = G*CC*wc;

  L = quaternion_rotate_vector(rr.q, Lbody);

  for(i=0; i<nbodies; i++) {
    L.x += GM[i]*(p[i].y*p[i].zd - p[i].z*p[i].yd);
    L.y += GM[i]*(p[i].z*p[i].xd - p[i].x*p[i].zd);
    L.z += GM[i]*(p[i].x*p[i].yd - p[i].y*p[i].xd);
  }

  return(L);
}

Vector jacobi_angmom()
{
  double wa, wb, wc;
  int i;

  Vector L;
  Vector Lbody;

  wa = rr.wa;
  wb = rr.wb;
  wc = rr.wc;

  Lbody.x = G*AA*wa;
  Lbody.y = G*BB*wb;
  Lbody.z = G*CC*wc;

  L = quaternion_rotate_vector(rr.q, Lbody);

  for(i=0; i<nbodies-1; i++) {
    L.x += GMj[i]*(pj[i].y*pj[i].zd - pj[i].z*pj[i].yd);
    L.y += GMj[i]*(pj[i].z*pj[i].xd - pj[i].x*pj[i].zd);
    L.z += GMj[i]*(pj[i].x*pj[i].yd - pj[i].y*pj[i].xd);
  }

  return(L);
}

Vector jacobi_satellite_angmom()
{
  double wa, wb, wc;
  int i;

  Vector L;
  Vector Lbody;

  wa = rr.wa;
  wb = rr.wb;
  wc = rr.wc;

  Lbody.x = G*AA*wa;
  Lbody.y = G*BB*wb;
  Lbody.z = G*CC*wc;

  L = quaternion_rotate_vector(rr.q, Lbody);

  for(i=0; i<nsatellites; i++) {
    L.x += GMj[i]*(pj[i].y*pj[i].zd - pj[i].z*pj[i].yd);
    L.y += GMj[i]*(pj[i].z*pj[i].xd - pj[i].x*pj[i].zd);
    L.z += GMj[i]*(pj[i].x*pj[i].yd - pj[i].y*pj[i].xd);
  }

  return(L);
}

Vector computeL()
{
  Vector L;
  void dolinks(), undolinks();

  undolinks();
  L = angmom();
  dolinks();
  
  return(L);
}

Vector compute_jacobi_L()
{
  Vector L;
  void dolinks(), undolinks();

  undolinks();
  L = jacobi_angmom();
  dolinks();
  
  return(L);
}

Vector compute_satellite_jacobi_L()
{
  Vector L;
  void dolinks(), undolinks();

  undolinks();
  L = jacobi_satellite_angmom();
  dolinks();
  
  return(L);
}

