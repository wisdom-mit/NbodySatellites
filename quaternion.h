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

typedef struct {
  double r, x, y, z;
} Quaternion;

typedef struct {
  double t;
  double r, x, y, z;
  double rdot, xdot, ydot, zdot;
} QuaternionState;

typedef struct {
  double t;
  Quaternion q;
  double wa, wb, wc;
} QWState;

/*
typedef struct {
  double x, y, z;
} Vector;
*/

typedef struct {
  double xx, xy, xz;
  double yx, yy, yz;
  double zx, zy, zz;
} Matrix;

