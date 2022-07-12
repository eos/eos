/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Viktor Kuschke
 *
 * This code is based on code published by Hjalte Frellesvig, Damiano Tommasini
 * and Christopher Wever, which is in the public domain.
 * For further details see ArXiv:1601.02649.
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef EOS_GUARD_EOS_MATHS_MULTIPLEPOLYLOG_LI22_CONST_HH
#define EOS_GUARD_EOS_MATHS_MULTIPLEPOLYLOG_LI22_CONST_HH 1

#include <eos/maths/complex.hh>

namespace eos
{
    namespace li22_impl
    {
        extern const double ccli2logA10;
        extern const double ccli2logA12;
        extern const double ccli2logA1log;
        extern const double ccli3logA10;
        extern const double ccli3logA11;
        extern const double ccli3logA13;
        extern const double ccli3logA1log;
        extern const double ccli4logA10;
        extern const double ccli4logA12;
        extern const double ccli4logA14;
        extern const double ccli4logA1log;

        // FIXME
        extern const double ccli2logA0[100];
        extern const double ccli2logA1[100];
        extern const double ccli3logA1[100];
        extern const double ccli4logA1[100];

        // BELOW ARE THE COEFFICIENTS FOR LI22FAST
        extern const double ccli221[100];
        extern const double ccli222[100];

        // BELOW ARE THE COEFFICIENTS FOR LI22STUFFLE (almost the same as above)
        extern const double ccli22stuffle1[100];
        extern const double ccli22stuffle2[100];

        // BELOW ARE THE CONSTANTS FOR LI22INV

        extern const double ccli234fast2o[50];
        extern const double ccli234fast3o[50];
        extern const double ccli234fast4o[50];

        extern const double ccli234fast2e;
        extern const double ccli234fast3e[50];
        extern const double ccli234fast4e[50];

        extern const double ccli22invc1;
        extern const double ccli22invc2;
        extern const double ccli22invc3;
        extern const double ccli22invc4;
        extern const double ccli22invc5;
        extern const double ccli22invc6;
        extern const double ccli22invc7;
        extern const double ccli22invc8;
        extern const double ccli22invc9;
        extern const double ccli22invc10;

        // BELOW ARE THE CONSTANTS USED FOR LI22BERNOULLI
        extern const double ccli22logA0a[101][101];
        extern const double ccli22logA0b[101][101];
        extern const double ccli22logA0c[101];
        extern const double ccli22logA0d[101];

        // BELOW ARE THE CONSTANTS USED FOR LI22 logA1 (FF2)
        extern const double ccli22logA1k1[100][100];
        extern const double ccli22logA1k2[100][100];
        extern const double ccli22logA1k3[100][100];
        extern const double ccli22logA1k4[100];
        extern const double ccli22logA1k5[100];

        // BELOW ARE THE CONSTANTS USED FOR logA1 ff1
        extern const double cclogA1ff11[100];
        extern const double cclogA1ff12[100];
        extern const double cclogA1ff1c1;
        extern const double cclogA1ff1c2;

        /* Constants for the Li functions used in Holder. */

        extern const double ccholderfo1[100];
        extern const double ccholderfo2[100];
        extern const double ccholderfo3[100];
        extern const double ccholderfo4[100];
        extern const double ccholder11[100];
        extern const double ccholder12[100];
        extern const double ccholder32[100];

        extern const double *ccholder21;
        extern const double *ccholder22;
        extern const double *ccholder23;
        extern const double *ccholder24;
        extern const double *ccholder31;
        extern const double *ccholder41;
        extern const double *ccholder42;
        extern const double *ccholder51;
        extern const double *ccholder52;
        extern const double *ccholder53;

        extern const double ccholderc11;
        extern const double ccholderc12;
        extern const double ccholderc21;
        extern const double ccholderc22;
        extern const double ccholderc31;
        extern const double ccholderc32;
        extern const double ccholderc41;
        extern const double ccholderc42;
        extern const double ccholderc51;
        extern const double ccholderc52;
        extern const double ccholderclog2;

        extern const double ccli22diagonallog[150];
        extern const double ccli22diagonallit[150];
        extern const double ccli22diagonalpow[150][150];
    }
}

#endif
