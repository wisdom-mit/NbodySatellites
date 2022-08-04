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

double dot_product(Vector u, Vector v)
{
  return(u.x*v.x + u.y*v.y + u.z*v.z);
}

Vector cross_product(Vector u, Vector v)
{
  Vector w;

  w.x = u.y*v.z - u.z*v.y;
  w.y = u.z*v.x - u.x*v.z;
  w.z = u.x*v.y - u.y*v.x;

  return(w);
}

void rotate_vector_z(double theta, Vector v, Vector *rv)
{
  double c, s;

  c = cos(theta);
  s = sin(theta);


  rv->x = v.x*c - v.y*s;
  rv->y = v.x*s + v.y*c;
  rv->z = v.z;
}

void rotate_vector_x(double theta, Vector v, Vector *rv)
{
  double c, s;

  c = cos(theta);
  s = sin(theta);

  rv->x = v.x;
  rv->y = v.y*c - v.z*s;
  rv->z = v.y*s + v.z*c;
}

void rotate_vector_y(double theta, Vector v, Vector *rv)
{
  double c, s;

  c = cos(theta);
  s = sin(theta);

  rv->x = v.z*s + v.x*c;
  rv->y = v.y;
  rv->z = v.z*c - v.x*s;
}

