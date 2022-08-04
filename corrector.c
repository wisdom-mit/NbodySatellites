/*

Copyright (C) 1996-2022 Jack Wisdom

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

double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16;
double b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16;
double alpha, beta;

void compute_corrector_coefficients(double dt)
{
  alpha = sqrt(7.0/40.0);
  beta = 1./(48.0*alpha);

  a1  = alpha * (  1.0 );
  a2  = alpha * ( -1.0 );
  a3  = alpha * (  2.0 );
  a4  = alpha * ( -2.0 );
  a5  = alpha * (  3.0 );
  a6  = alpha * ( -3.0 );
  a7  = alpha * (  4.0 );
  a8  = alpha * ( -4.0 );
  a9  = alpha * (  5.0 );
  a10 = alpha * ( -5.0 );
  a11 = alpha * (  6.0 );
  a12 = alpha * ( -6.0 );
  a13 = alpha * (  7.0 );
  a14 = alpha * ( -7.0 );
  a15 = alpha * (  8.0 );
  a16 = alpha * ( -8.0 );


  b1  = beta * (  1.8685517340134143 );
  b2  = beta * ( -1.8685517340134143 );
  b3  = beta * ( -1.3090623112714728 );
  b4  = beta * (  1.3090623112714728 );
  b5  = beta * (  .6510325862986641 );
  b6  = beta * ( -.6510325862986641 );
  b7  = beta * ( -.24239903351841396 );
  b8  = beta * (  .24239903351841396 );
  b9  = beta * (  .06652968674402478 );
  b10 = beta * ( -.06652968674402478 );
  b11 = beta * ( -1.2770775246667285e-2 );
  b12 = beta * (  1.2770775246667285e-2 );
  b13 = beta * (  1.5348298318361457e-3 );
  b14 = beta * ( -1.5348298318361457e-3 );
  b15 = beta * ( -8.704091947232721e-5 );
  b16 = beta * (  8.704091947232721e-5 );

  a1 *= dt;
  a2 *= dt;
  a3 *= dt;
  a4 *= dt;
  a5 *= dt;
  a6 *= dt;
  a7 *= dt;
  a8 *= dt;
  a9 *= dt;
  a10 *= dt;
  a11 *= dt;
  a12 *= dt;
  a13 *= dt;
  a14 *= dt;
  a15 *= dt;
  a16 *= dt;

  b1 *= dt;
  b2 *= dt;
  b3 *= dt;
  b4 *= dt;
  b5 *= dt;
  b6 *= dt;
  b7 *= dt;
  b8 *= dt;
  b9 *= dt;
  b10 *= dt;
  b11 *= dt;
  b12 *= dt;
  b13 *= dt;
  b14 *= dt;
  b15 *= dt;
  b16 *= dt;
}

void real_to_map()
{
  void Z();

  Z(a16, b16);
  Z(a14, b14);
  Z(a12, b12);
  Z(a10, b10);
  Z(a8, b8);
  Z(a6, b6);
  Z(a4, b4);
  Z(a2, b2);
  Z(a1, b1);
  Z(a3, b3);
  Z(a5, b5);
  Z(a7, b7);
  Z(a9, b9);
  Z(a11, b11);
  Z(a13, b13);
  Z(a15, b15);
}


void map_to_real()
{
  void Z();

  Z(a15, -b15);
  Z(a13, -b13);
  Z(a11, -b11);
  Z(a9, -b9);
  Z(a7, -b7);
  Z(a5, -b5);
  Z(a3, -b3);
  Z(a1, -b1);
  Z(a2, -b2);
  Z(a4, -b4);
  Z(a6, -b6);
  Z(a8, -b8);
  Z(a10, -b10);
  Z(a12, -b12);
  Z(a14, -b14);
  Z(a16, -b16);
}

void C(a, b)
     double a, b;
{
  void A();
  void B();

  A(-a);
  B(b);
  A(a);
}

void Z(double a, double b)
{
  void C();

  C(-a, -b);
  C(a, b);
}

