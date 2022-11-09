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

#include <eos/maths/complex.hh>
#include <eos/maths/power-of.hh>
#include <eos/maths/polylog.hh>
#include <eos/maths/multiplepolylog-li22.hh>
#include <eos/maths/multiplepolylog-li22-const.hh>
#include <eos/utils/exception.hh>

//#include <stdio.h>
#include <iostream>

namespace eos
{
    namespace li22_impl
    {
        // #define USEDEBUG
        #ifdef USEDEBUG
        const char *debugstring[100];
        #define debug(x) *debugstring = x;
        #else
        #define debug(x)
        #endif

        #define LI2LOGA0MAX 8
        //Checked
        #define LINLOGA1MAX 8
        //Checked

        complex<double> li2logA0rec(complex<double> next, int nc, const complex<double> &lsq){
            if(nc<LI2LOGA0MAX){
                return next*ccli2logA0[nc] + li2logA0rec(next*lsq, nc+1, lsq);
            }
            return next*ccli2logA0[nc];
        }

        complex<double> li2logA0(const complex<double> &x){
            const complex<double> mlog = -std::log(1.0-x);
            const complex<double> mlsq = mlog*mlog;
            return mlog - 0.25*mlsq + li2logA0rec(mlog*mlsq, 1, mlsq);
        }

        complex<double> li22rec(complex<double> next, complex<double> extra, int nc,const complex<double> &a,
            const complex<double> &ab, int &ncmax)
        {
            if(nc<ncmax){
                return next + li22rec(ccli221[nc]*a*next + ccli222[nc]*extra, extra*ab, nc+1, a, ab, ncmax);
            }
            else{
                return next + ccli221[nc]*a*next + ccli222[nc]*extra;
            }
        }

        complex<double> li22fast(const complex<double> &a, const complex<double> &b, int &ncmax)
        {
            const complex<double> ab = a*b;
            return li22rec(0.25*a*ab, a*ab*ab, 0, a, ab, ncmax);
        }

        // #define LI22STUFFLEMAX LI22FASTMAX+1

        complex<double> li22stufflerec(complex<double> next, complex<double> extra, int nc, const complex<double> &b,
            const complex<double> &ab, int &ncmax)
        {
            if(nc<ncmax){
                return next + li22stufflerec(ccli22stuffle1[nc]*b*next + ccli22stuffle2[nc]*extra,
                    extra*ab, nc+1, b, ab, ncmax);
            }
            return next + ccli22stuffle1[nc]*b*next + ccli22stuffle2[nc]*extra;
        }

        complex<double> li22stuffle(const complex<double> &a, const complex<double> &b, int &ncmax)
        {
            const complex<double> ab = a*b;
            return -li22stufflerec(ab-b*dilog(a), ab*ab, 0, b, ab, ncmax);
        }

        // This calculates li22(a,b) + li4(a*b) + q*li2(a)
        complex<double> li224q(const complex<double> &a, const complex<double> &b, complex<double> q, int &ncmax)
        {
            const complex<double> ab = a*b;
            return li22stufflerec(ab+a*q, ab*ab, 0, a, ab, ncmax);
        }

        void threelicrrec(complex<double> next, int nc, const complex<double> &xsq, complex<double> *li2var,
            complex<double> *li3var, complex<double> *li4var)
        {
            *li2var += ccli2logA1[nc]*next;
            *li3var += ccli3logA1[nc]*next;
            *li4var += ccli4logA1[nc]*next;

            if(nc<LINLOGA1MAX)
            {
                threelicrrec(next*xsq, nc+1, xsq, li2var, li3var, li4var);
            }
        }

        complex<double> threelicr(const complex<double> &a, complex<double> c2, complex<double> c3){

            const complex<double> x = -std::log(a);
            const complex<double> xsq = x*x;
            const complex<double> logx = std::log(x);
            complex<double> li2var = ccli2logA10 + ccli2logA12*xsq + x*logx;
            complex<double> li3var = ccli3logA11 + xsq*ccli3logA13 + ccli3logA1log*x*logx;
            complex<double> li4var = ccli4logA10 + xsq*(ccli4logA12 + xsq*ccli4logA14 + x*ccli4logA1log*logx);

            threelicrrec(x, 0, xsq, &li2var, &li3var, &li4var);

            li3var*=x;
            li3var+=ccli3logA10;
            return c2*li2var + c3*li3var + 3.*li4var;
        }

        #define LIALLBERNMAX 9

        void threelibeoddrec(complex<double> next, int nc, const complex<double> &xsq, complex<double> *li2var,
            complex<double> *li3var, complex<double> *tli4var)
        {
            *li2var += ccli234fast2o[nc]*next;
            *li3var += ccli234fast3o[nc]*next;
            *tli4var += ccli234fast4o[nc]*next;

            if(nc<LIALLBERNMAX)
            {
                threelibeoddrec(next*xsq, nc+1, xsq, li2var, li3var, tli4var);
            }
        }

        void threelibeevenrec(complex<double> next, int nc, const complex<double> &xsq, complex<double> *li3var,
            complex<double> *tli4var)
        {
            *li3var += ccli234fast3e[nc]*next;
            *tli4var += ccli234fast4e[nc]*next;

            if(nc<LIALLBERNMAX)
            {
                threelibeevenrec(next*xsq, nc+1, xsq, li3var, tli4var);
            }
        }

        complex<double> threelibe(const complex<double> &a, complex<double> c2, complex<double> c3)
        {
            const complex<double> x = -std::log(1.-a);
            const complex<double> xsq = x*x;
            complex<double> li2var = x + ccli234fast2e*xsq;
            complex<double> li3var = x;
            complex<double> tli4var = 3.*x;

            threelibeevenrec(xsq,0,xsq,&li3var,&tli4var);
            threelibeoddrec(x*xsq,1,xsq,&li2var,&li3var,&tli4var);

            return c2*li2var+c3*li3var+tli4var;
        }

        void sixlimmmbeoddrec(complex<double> nextx, complex<double> nexty, int nc, const complex<double> &xsq,
            const complex<double> &ysq, complex<double> *li2dif, complex<double> *li3dif, complex<double> *tli4dif)
        {
            const complex<double> dif = nextx-nexty;
            *li2dif += ccli234fast2o[nc]*dif;
            *li3dif += ccli234fast3o[nc]*dif;
            *tli4dif += ccli234fast4o[nc]*dif;

            if(nc<LIALLBERNMAX){
                sixlimmmbeoddrec(nextx*xsq, nexty*ysq, nc+1, xsq, ysq, li2dif, li3dif, tli4dif);
            }
        }

        void sixlimmmbeevenrec(complex<double> nextx, complex<double> nexty, int nc, const complex<double> &xsq,
            const complex<double> &ysq, complex<double> *li3dif, complex<double> *tli4dif)
        {
            const complex<double> dif = nextx-nexty;
            *li3dif += ccli234fast3e[nc]*dif;
            *tli4dif += ccli234fast4e[nc]*dif;

            if(nc<LIALLBERNMAX)
            {
                sixlimmmbeevenrec(nextx*xsq, nexty*ysq, nc+1, xsq, ysq, li3dif, tli4dif);
            }
        }

        complex<double> sixlimmmbe(const complex<double> &a, const complex<double> &b, complex<double> c2, complex<double> c3)
        {
            const complex<double> x = -std::log(1.-a);
            const complex<double> xsq = x*x;
            const complex<double> y = -std::log(1.-b);
            const complex<double> ysq = y*y;
            const complex<double> xydif = x-y;
            complex<double> li2dif = xydif + ccli234fast2e*(xsq-ysq);
            complex<double> li3dif = xydif;
            complex<double> tli4dif = 3.*xydif;

            sixlimmmbeevenrec(xsq,ysq,0,xsq,ysq,&li3dif,&tli4dif);
            sixlimmmbeoddrec(x*xsq,y*ysq,1,xsq,ysq,&li2dif,&li3dif,&tli4dif);

            return c2*li2dif+c3*li3dif+tli4dif;
        }


        void sixlipmpbeoddrec(complex<double> nextx, complex<double> nexty, int nc, const complex<double> &xsq,
            const complex<double> &ysq, complex<double> *li2sum, complex<double> *li3dif, complex<double> *tli4sum)
        {
            const complex<double> sum = nextx+nexty;
            *li2sum += ccli234fast2o[nc]*sum;
            *li3dif += ccli234fast3o[nc]*(nextx-nexty);
            *tli4sum += ccli234fast4o[nc]*sum;

            if(nc<LIALLBERNMAX)
            {
                sixlipmpbeoddrec(nextx*xsq, nexty*ysq, nc+1, xsq, ysq, li2sum, li3dif, tli4sum);
            }
        }

        void sixlipmpbeevenrec(complex<double> nextx, complex<double> nexty, int nc, const complex<double> &xsq,
            const complex<double> &ysq, complex<double> *li3dif, complex<double> *tli4sum)
        {
            *li3dif += ccli234fast3e[nc]*(nextx-nexty);
            *tli4sum += ccli234fast4e[nc]*(nextx+nexty);

            if(nc<LIALLBERNMAX)
            {
                sixlipmpbeevenrec(nextx*xsq, nexty*ysq, nc+1, xsq, ysq, li3dif, tli4sum);
            }
        }

        //This function calculates c2*(li2(a)+li2(b)) + c3*(li3(a)-li3(b)) - 3*(li4(a)+li4(b))  and calls the two above.
        complex<double> sixlipmpbe(const complex<double> &a, const complex<double> &b, complex<double> c2, complex<double> c3)
        {
            const complex<double> x = -std::log(1.-a);
            const complex<double> xsq = x*x;
            const complex<double> y = -std::log(1.-b);
            const complex<double> ysq = y*y;
            const complex<double> xysum = x+y;
            complex<double> li2sum = xysum + ccli234fast2e*(xsq+ysq);
            complex<double> li3dif = x-y;
            complex<double> tli4sum = 3.*xysum;

            sixlipmpbeevenrec(xsq,ysq,0,xsq,ysq,&li3dif,&tli4sum);
            sixlipmpbeoddrec(x*xsq,y*ysq,1,xsq,ysq,&li2sum,&li3dif,&tli4sum);

            return c2*li2sum+c3*li3dif-tli4sum;
        }

        complex<double> sixlimmm(const complex<double> &a, const complex<double> &b, complex<double> c2, complex<double> c3)
        {
            if(a.real()<0.5){
                if(b.real()<0.5){
                    return sixlimmmbe(a,b,c2,c3);
                }
                else{
                    return (threelibe(a,c2,c3)-threelicr(b,c2,c3));
                }
            }
            else{
                if(b.real()<0.5){
                    return (threelicr(a,c2,c3)-threelibe(b,c2,c3));
                }
                else{
                    return (threelicr(a,c2,c3)-threelicr(b,c2,c3));
                }
            }
        }

        complex<double> sixlipmp(const complex<double> &a, const complex<double> &b, complex<double> c2, complex<double> c3)
        {
            const complex<double> mc2 = -c2;
            const complex<double> mc3 = -c3;

            if(a.real()<0.5){
                if(b.real()<0.5){
                    return sixlipmpbe(a,b,c2,c3);
                }
                else{
                    return -(threelibe(a,mc2,mc3)+threelicr(b,mc2,c3));
                }
            }
            else{
                if(b.real()<0.5){
                    return -(threelicr(a,mc2,mc3)+threelibe(b,mc2,c3));
                }
                else{
                    return -(threelicr(a,mc2,mc3)+threelicr(b,mc2,c3));
                }
            }
        }

        //  This is the inversion formula for li22(x,y)->li22(1/x,1/y)+...
        complex<double> li22inv(const complex<double> &x, const complex<double> &y, int &ncmax)
        {
            const complex<double> logmxy = std::log(-x*y);
            const complex<double> logmx = std::log(-x);
            const complex<double> logmxysq = logmxy*logmxy;
            const complex<double> logmxsq = logmx*logmx;
            const complex<double> ix = 1./x;
            const complex<double> iy = 1./y;

            return li224q(ix,iy,ccli22invc10+ccli22invc9*logmxsq,ncmax) +
                sixlimmm(ix, iy, ccli22invc9*(logmxysq-logmxsq), ccli22invc8*logmxy) +
                ccli22invc1 + logmxy*(ccli22invc4*logmx + ccli22invc6*logmxy) +
                logmxsq*(ccli22invc2 + ccli22invc3*logmxsq + ccli22invc5*logmx*logmxy + ccli22invc7*logmxysq);
        }


        // This function calculates li22(a,b) + q*li2(a)
        complex<double> li22qfast(const complex<double> &a, const complex<double> &b, complex<double> q, int &ncmax)
        {
            complex<double> ab = a*b;
            return li22rec(0.25*a*a*(b+q), a*ab*ab, 0, a, ab, ncmax) + a*q;
        }


        const double ccli22invstuffle1 = 0.5;
        const double ccli22invstuffle2 = 1.64493406684822643647;     // pi^2/6
        const double ccli22invstuffle3 = 1.64493406684822643647;     // pi^2/6
        const double ccli22invstuffle4 = 0.5;
        const double ccli22invstuffle5 = 2.;
        const double ccli22invstuffle6 = 0.0416666666666666666667;   // 15/360
        const double ccli22invstuffle7 = -0.25;                      // -15*6/360
        const double ccli22invstuffle8 = 0.25;                       //  15*6/360
        const double ccli22invstuffle9 = 0.333333333333333333333;    // 8*15/360
        const double ccli22invstuffle10 = -0.125;                    // -15*3/360
        const double ccli22invstuffle11 = 0.822467033424113218236;   // 30*pi^2/360
        const double ccli22invstuffle12 = 3.28986813369645287294;    // 30*4*pi^2/360
        const double ccli22invstuffle13 = -2.46740110027233965471;   // -30*3*pi^2/360
        const double ccli22invstuffle14 = -3.78813131798898367031;   // -14*pi^4/360

        // This is inversion and stuffle for Li22 (for |y|>1).
        complex<double> li22invstuffle(const complex<double> &x, const complex<double> &y, int &ncmax)
        {
        // To be used when x>1 and y>1

            const complex<double> logmx = std::log(-x);
            const complex<double> logmy = std::log(-y);
            const complex<double> logmxy = std::log(-x*y);
            const complex<double> logmxsq = logmx*logmx;
            const complex<double> logmysq = logmy*logmy;
            const complex<double> logmxysq = logmxy*logmxy;
            const complex<double> ix = 1./x;
            const complex<double> iy = 1./y;
        //  double pi = 3.141592653589793238462643;

        //   return -li22qfast(iy,ix,li2(x),ncmax) + sixlimmm(ix,iy,pi*pi/6.+0.5*logmxysq,2.*logmxy) + (15.*(logmxysq*logmxysq-6.*logmxysq*logmysq+(6.*logmxsq+8.*logmxy*logmy-3.*logmysq)*logmysq)+30.*(logmxsq+4.*logmxy*logmy-3.*logmysq)*pi*pi-14.*pi*pi*pi*pi)/360.;

            return -li22qfast(iy,ix,-dilog(ix),ncmax) + dilog(iy)*(ccli22invstuffle1*logmxsq+ccli22invstuffle2) +
                sixlimmm(ix, iy, ccli22invstuffle3 + ccli22invstuffle4*logmxysq, ccli22invstuffle5*logmxy) +
                ccli22invstuffle6*logmxysq*logmxysq + ccli22invstuffle7*logmxysq*logmysq +
                (ccli22invstuffle8*logmxsq + ccli22invstuffle9*logmxy*logmy+ccli22invstuffle10*logmysq)*logmysq +
                ccli22invstuffle11*logmxsq + ccli22invstuffle12*logmxy*logmy + ccli22invstuffle13*logmysq +
                ccli22invstuffle14;
        }

        const double ccli22invsp1 = 1.64493406684822643647;      // Pi^2/6
        const double ccli22invsp2 = 0.5;
        const double ccli22invsp3 = 0.5;
        const double ccli22invsp4 = -2.;
        const double ccli22invsp5 = 0.822467033424113218236;     // Pi^2/12
        const double ccli22invsp6 = 0.0416666666666666666667;    // 1/24
        const double ccli22invsp7 = 1.89406565899449183515;      // 7*Pi^4/360

        // This is Li22 inversion for |x|<1
        complex<double> li22invspecial(const complex<double> &x, const complex<double> &y, int &ncmax)
        {
            const complex<double> logmxy = std::log(-x*y);
            const complex<double> logmx = std::log(-x);
            const complex<double> logmxysq = logmxy*logmxy;
            const complex<double> logmxsq = logmx*logmx;
            const complex<double> ix = 1./x;
            const complex<double> iy = 1./y;
        //  double pi = 3.141592653589793238462643;

        //  return li224q(ix,iy,pi*pi/6.+0.5*logmxsq,ncmax) - sixlipmp(ix, y, 0.5*(logmxsq-logmxysq), -2.*logmxy) + (pi*pi/(12.) + 1/(24.)*logmxysq)*logmxysq + 7.*pi*pi*pi*pi/(360.);

            return li224q(ix,iy,ccli22invsp1+ccli22invsp2*logmxsq,ncmax) -
                sixlipmp(ix, y, ccli22invsp3*(logmxsq-logmxysq), ccli22invsp4*logmxy) +
                (ccli22invsp5 + ccli22invsp6*logmxysq)*logmxysq + ccli22invsp7;
        }


        const double ccli22invstufflesp1 = 0.5;
        const double ccli22invstufflesp2 = 1.64493406684822643647;   // pi^2/6
        const double ccli22invstufflesp3 = -1.64493406684822643647;  // -pi^2/6
        const double ccli22invstufflesp4 = -0.5;
        const double ccli22invstufflesp5 = 2.;
        const double ccli22invstufflesp6 = -0.125;                   // -3/24
        const double ccli22invstufflesp7 = 0.333333333333333333333;  // 8/24
        const double ccli22invstufflesp8 = -0.25;                    // -6/24
        const double ccli22invstufflesp9 = 0.0416666666666666666667; // 1/24
        const double ccli22invstufflesp10 = 0.25;                    //  6/24
        const double ccli22invstufflesp11 = -0.25;                   // -6/24
        const double ccli22invstufflesp12 = 0.333333333333333333333; // 8/24
        const double ccli22invstufflesp13 = -0.125;                  // -3/24
        const double ccli22invstufflesp14 = -2.46740110027233965471; //  -2*3*pi^2/24
        const double ccli22invstufflesp15 = -0.822467033424113218236;//  -2*pi^2/24
        const double ccli22invstufflesp16 = 3.28986813369645287294;  // 2*4*pi^2/24
        const double ccli22invstufflesp17 = -2.46740110027233965471; //  -2*3*pi^2/24
        const double ccli22invstufflesp18 = -12.1761363792503046546; //  -3*pi^4/24

        // This is li22 inv+stuffle for |y|>1
        complex<double> li22invstufflespecial(const complex<double> &x, const complex<double> &y, int &ncmax)
        {
        // To be used when x<1 and y>1

            const complex<double> logmx = std::log(-x);
            const complex<double> logmy = std::log(-y);
            const complex<double> logmxy = std::log(-x*y);
            const complex<double> logmxsq = logmx*logmx;
            const complex<double> logmysq = logmy*logmy;
            const complex<double> logmxysq = logmxy*logmxy;
            const complex<double> ix = 1./x;
            const complex<double> iy = 1./y;
        //  double pi = 3.141592653589793238462643;

        //   return -li22qfast(iy,ix,li2(x),ncmax) + sixlipmp(x,iy,-pi*pi/6.-0.5*logmxysq,2.*logmxy)+(-3.*logmxsq*logmxsq+8.*logmx*logmxsq*logmxy-6.*logmxsq*logmxysq+logmxysq*logmxysq+6.*logmxsq*logmysq-6.*logmxysq*logmysq+8.*logmxy*logmy*logmysq-3.*logmysq*logmysq-2.*(3.*logmxsq+logmxysq-4.*logmxy*(logmx+logmy)+3.*logmysq)*pi*pi-3.*pi*pi*pi*pi)/24.;

        //   return -li22qfast(iy,ix,-li2(ix),ncmax) + li2(iy)*(0.5*logmxsq+pi*pi/6.) + sixlipmp(x,iy,-pi*pi/6.-0.5*logmxysq,2.*logmxy)+(-3.*logmxsq*logmxsq+8.*logmx*logmxsq*logmxy-6.*logmxsq*logmxysq+logmxysq*logmxysq+6.*logmxsq*logmysq-6.*logmxysq*logmysq+8.*logmxy*logmy*logmysq-3.*logmysq*logmysq-2.*(3.*logmxsq+logmxysq-4.*logmxy*(logmx+logmy)+3.*logmysq)*pi*pi-3.*pi*pi*pi*pi)/24.;

            return -li22qfast(iy,ix,-dilog(ix),ncmax) + dilog(iy)*(ccli22invstufflesp1*logmxsq + ccli22invstufflesp2) +
                sixlipmp(x,iy,ccli22invstufflesp3 + ccli22invstufflesp4*logmxysq, ccli22invstufflesp5*logmxy) +
                ccli22invstufflesp6*logmxsq*logmxsq + ccli22invstufflesp7*logmx*logmxsq*logmxy +
                ccli22invstufflesp8*logmxsq*logmxysq + ccli22invstufflesp9*logmxysq*logmxysq +
                ccli22invstufflesp10*logmxsq*logmysq + ccli22invstufflesp11*logmxysq*logmysq +
                ccli22invstufflesp12*logmxy*logmy*logmysq + ccli22invstufflesp13*logmysq*logmysq +
                ccli22invstufflesp14*logmxsq + ccli22invstufflesp15*logmxysq + ccli22invstufflesp16*logmxy*(logmx+logmy) +
                ccli22invstufflesp17*logmysq + ccli22invstufflesp18;
        }

        complex<double> li22diagonalind(const complex<double> &iab, int &nmo, const complex<double> &logomx, const complex<double> &litx)
        {
            complex<double> res = 0.;
            complex<double> iabp = 1.;
            int i;

            for(i=0;i<=nmo;++i)
            {
                res += ccli22diagonalpow[nmo][i]*iabp;
                iabp *= iab;
            }
            res += ccli22diagonallog[nmo]*logomx*(1.-iabp);
            res += ccli22diagonallit[nmo]*litx*(1.+iabp);

            return res;
        }

        complex<double> li22diagonal(const complex<double> &a, const complex<double> &b, char ncmax)
        {
            complex<double> res = 0;
            const complex<double> ab = a*b;
            const complex<double> iab = 1./ab;
            const complex<double> logomx = std::log(1.-ab);
            const complex<double> litx = dilog(ab);
            complex<double> ap = 1.;
            int i;

            for(i=0;i<ncmax;++i)
            {
                ap *= a;
                res += ap*li22diagonalind(iab,i,logomx,litx);
            }

            return res;
        }

        // This is Li22(x/2,y)
        complex<double> holderf1rec(complex<double> part1, complex<double> part2, const complex<double> &x,
            const complex<double> &xy, int n, int &nn)
        {
            if(n<nn){
                return part1 + holderf1rec(part1*x*ccholder11[n]+part2, part2*xy*ccholder12[n], x, xy, n+1, nn);
            }
            return part1;
        }

        complex<double> holderf1(const complex<double> &x, const complex<double> &xy, int &nn)
        {
            return holderf1rec(x*xy*ccholderc11, x*xy*xy*ccholderc12, x, xy, 0, nn);
        }

        // This is Li1111(x/2, 1/x, z, 1/z)
        complex<double> holderf2rec(complex<double> part1, complex<double> part2, complex<double> part3,
            complex<double> part4, const complex<double> &x, const complex<double> &z, int n, int &nn)
        {
            if(n<nn){
                return part1 + holderf2rec(part1*x*ccholder21[n]+part2, part2*ccholder22[n]+part3,
                    part3*z*ccholder23[n]+part4, part4*ccholder24[n], x, z, n+1, nn);
            }
            return part1;
        }

        complex<double> holderf2(const complex<double> &x, const complex<double> &z, int &nn)
        {
            return holderf2rec(0., 0., z*x*ccholderc21, z*x*ccholderc22, x, z, 0, nn);
        }


        // This is Li12(x/2,y)
        complex<double> holderf3rec(complex<double> part1, complex<double> part2, const complex<double> &x,
            const complex<double> &xy, int n, int &nn)
        {
            if(n<nn){
                return part1 + holderf3rec(part1*x*ccholder31[n]+part2, part2*xy*ccholder32[n], x, xy, n+1, nn);
            }
            return part1;
        }

        complex<double> holderf3(const complex<double> &x, const complex<double> &xy, int &nn)
        {
            return holderf3rec(x*xy*ccholderc31, x*xy*xy*ccholderc32, x, xy, 0, nn);
        }

        // This is Li11(x/2, 1/x)
        complex<double> holderf4rec(complex<double> part1, complex<double> part2, const complex<double> &x, int n, int &nn)
        {
            if(n<nn){
                return part1 + holderf4rec(part1*x*ccholder41[n]+part2, part2*ccholder42[n], x, n+1, nn);
            }
            return part1;
        }

        complex<double> holderf4(const complex<double> &x, int &nn)
        {
            return holderf4rec(x*ccholderc41, x*ccholderc42, x, 0, nn);
        }

        // This is -Li111(1/2,x,1/x)
        complex<double> holderf5rec(complex<double> part1, complex<double> part2, complex<double> part3,
            const complex<double> &y, int n, int &nn)
        {
            if(n<nn){
                return part1 + holderf5rec(part1*ccholder51[n]+part2, y*part2*ccholder52[n]+part3,
                    part3*ccholder53[n], y, n+1, nn);
            }
            return part1;
        }

        complex<double> holderf5(const complex<double> &y, int &nn)
        {
            return -holderf5rec(0., y*ccholderc51, y*ccholderc52, y, 0, nn);
        }

        // This is the Holder relation with q=2 used for othervise slowly convergent regions.
        complex<double> li22holder(const complex<double> &x, const complex<double> &y, int ncmax)
        {
            const complex<double> xy = x*y;
            const complex<double> r1 = xy/(xy-1.);
            const complex<double> r2 = x/(x-1.);
            const complex<double> xyh = xy*0.5;

            return holderf1(x,xy,ncmax) + holderf2(r1,r2,ncmax) + ccholderclog2*holderf3(x,xy,ncmax) -
                holderf4(r2,ncmax)*li2logA0(xyh) - holderf5(r2,ncmax)*std::log(1.-xyh);
        }

        complex<double> li22holderstuffle(const complex<double> &x, const complex<double> &y, int ncmax)
        {
            const complex<double> xy = x*y;

            return dilog(x)*dilog(y)-li22holder(y,x,ncmax)-quadlog(xy);
        }

        complex<double> li22logA0recint(complex<double> next, int nc, int &nco, const complex<double> &p,
            const complex<double> &logapp, int &ncmax)
        {
            if(nc<ncmax){
                return next*(ccli22logA0a[nc][nco] + ccli22logA0b[nc][nco]*logapp) +
                    li22logA0recint(next*p, nc+1, nco, p, logapp, ncmax);
            }
            return next*(ccli22logA0a[nc][nco] + ccli22logA0b[nc][nco]*logapp);
        }

        complex<double> li22logA0rec(complex<double> next, int nc, const complex<double> &t,
            const complex<double> &p, const complex<double> &logapp, const complex<double> &logompt, const complex<double> &litt, int &ncmax)
        {
            if(nc<ncmax){
                return next*(ccli22logA0c[nc]*logompt + ccli22logA0d[nc]*litt) +
                    li22logA0recint(next*p, 1, nc, p, logapp, ncmax) +
                    li22logA0rec(next*t, nc+1, t, p, logapp, logompt, litt, ncmax);
            }

            return next*(ccli22logA0c[nc]*logompt+ccli22logA0d[nc]*litt) +
                li22logA0recint(next*p, 1, nc, p, logapp, ncmax);
        }

        // This is the log expansion of li22 usable around (0., 0.).
        complex<double> li22logA0(const complex<double> &x, const complex<double> &y, int ncmax)
        {
            const complex<double> xy = x*y;
            const complex<double> a = 1./xy;
            const complex<double> p = -std::log(1.-xy);
            const complex<double> t = -std::log(1.-y);
            const complex<double> logapp = std::log(a) + std::log(p);
            const complex<double> logompt = std::log(1.-p/t);
            const complex<double> pot = p/t;
            const complex<double> litt = dilog(pot) + logapp*logompt;

            return li22logA0rec(1., 0, t, p, logapp, logompt, litt, ncmax);
        }


        const complex<double> cctwopii(0.,6.28318530717958647692528676656);
        const double ccxizero = 1.;
        const double ccxzero = exp(-ccxizero);

        complex<double> li22logA1recint(complex<double> next, int nc, int &nco,
            const complex<double> &b, const complex<double> &logb, const complex<double> &logbxi, int &ncmax)
        {
            if(nc<ncmax){
                return next*(ccli22logA1k1[nc][nco] + ccli22logA1k2[nc][nco]*logb +
                    ccli22logA1k3[nc][nco]*logbxi) + li22logA1recint(next*b, nc+1, nco, b, logb, logbxi,ncmax);
            }

            return next*(ccli22logA1k1[nc][nco] + ccli22logA1k2[nc][nco]*logb + ccli22logA1k3[nc][nco]*logbxi);
        }

        complex<double> li22logA1rec(complex<double> next, int nc, const complex<double> &a,
            const complex<double> &b, const complex<double> &logb, const complex<double> &logbxi, const complex<double> &logabdif, int &ncmax)
        {
            if(nc<ncmax){
                return next*(ccli22logA1k4[nc] + b*ccli22logA1k5[nc])*logabdif +
                li22logA1recint(next, 0, nc, b, logb, logbxi, ncmax) +
                li22logA1rec(next*a, nc+1, a, b, logb, logbxi, logabdif, ncmax);
            }

            return next*(ccli22logA1k4[nc] + b*ccli22logA1k5[nc])*logabdif +
            li22logA1recint(next, 0, nc, b, logb, logbxi, ncmax);
        }

        inline double signum(double x)
        {
            return ((x > 0.) ? 1. : ((x < 0.) ? -1. : 0.));
        }

        complex<double> li22logA1ff2(const complex<double> &x, const complex<double> &y, int &ncmax)
        {
            const complex<double> a = std::log(y);
            const complex<double> b = -std::log(x*y);
            const complex<double> apb = a+b;
            const complex<double> logb = std::log(b);
            const complex<double> logbxi = std::log(b+ccxizero);
            const complex<double> logabdif = std::log(apb) - std::log(apb+ccxizero);
            const complex<double> li2a1 = -(b+ccxizero)/a;
            const complex<double> li2a2 = -b/a;
            complex<double> res;

            res = li22logA1rec(1.,0,a,b,logb,logbxi,logabdif, ncmax) + apb*a*( dilog(li2a1) - dilog(li2a2) +
                logbxi*std::log((apb+ccxizero)/a) - logb*std::log(apb/a));

            if((abs(arg(b+ccxizero))<abs(arg(-a))) && (abs(arg(-a)) < abs(arg(b)))){
                if(arg(b)*arg(-a) > 0. && 1. < abs((b.imag())/(a.imag()))){
                    res += cctwopii*signum(a.imag())*apb*a*log(-a);
            }}
            return res;
        }

        complex<double> crli12rec(complex<double> part1, complex<double> part2, const complex<double> &x,
            const complex<double> &xy, int n, int &nn)
        {
            if(n<nn){
                return part1 + crli12rec(part1*x*cclogA1ff11[n]+part2, part2*xy*cclogA1ff12[n], x, xy, n+1, nn);
            }
            return part1;
        }

        complex<double> crli12(const complex<double> &x, const complex<double> &y, int nn)
        {
            const complex<double> xy = x*y;
            return crli12rec(x*xy*cclogA1ff1c1, x*xy*xy*cclogA1ff1c2, x, xy, 0, nn);
        }

        complex<double> li22logA1ff1(const complex<double> &x, const complex<double> &y)
        {
            const complex<double> xp = x*ccxzero;
            int nmax = 90;
        //   return li22(xp,y) + ccxizero*crli12(xp,y,100);
            return li22fast(xp,y,nmax) + ccxizero*crli12(xp,y,100);
        }

        // This is the log expansion of li22 usable around (1., 1.) .

        complex<double> li22logA1(const complex<double> &x, const complex<double> &y, int ncmax)
        {
            if(ncmax>99){printf("WRONG!!");}
            return li22logA1ff1(x,y) + li22logA1ff2(x,y,ncmax);
        }

        // Log expansion around 1 with stuffle.
        complex<double> li22logA1stuffle(const complex<double> &x, const complex<double> &y, int ncmax)
        {
            const complex<double> xy = x*y;
            return dilog(x)*dilog(y) - li22logA1(y,x,ncmax) - quadlog(xy);
        }

        complex<double> li22smalla(const complex<double> &a, const complex<double> &b)
        {
            const double absa = abs(a);

            if(absa<0.1){
                return li22diagonal(a,b,18);
            }
            if(absa<0.2){
                return li22diagonal(a,b,25);
            }
            if(absa<0.3){
                return li22diagonal(a,b,33);
            }
            if(absa<0.4){
                return li22diagonal(a,b,43);
            }

            if(abs(arg(a))>0.4 && abs(arg(a*b))>0.9){
                return li22holder(a,b,95);
            }

            if(absa<0.5){
                return li22diagonal(a,b,56);
            }
            if(absa<0.6){
                return li22diagonal(a,b,77);
            }
            if(absa<0.7){
                return li22diagonal(a,b,112);
            }

            throw InternalError("Incorrect use of li22smalla");
            return 0.;
        }

        complex<double> li22smallastuffle(const complex<double> &a, const complex<double> &b)
        {
            const complex<double> ab = a*b;

            return dilog(a)*dilog(b) - quadlog(ab) - li22smalla(b,a);
        }

        complex<double> li22smallainv(const complex<double> &x, const complex<double> &y)
        {
            const complex<double> xy = x*y;
            const complex<double> ix = 1./x;
            const complex<double> iy = 1./y;
            const complex<double> logmxy = std::log(-xy);
            const complex<double> logmx = std::log(-x);
            const double pisqos = 1.64493406684822643647;

            return li22smalla(ix,iy) - quadlog(xy) + 3.*(quadlog(y) + quadlog(ix)) + 2.*(trilog(ix) - trilog(y))*logmxy +
                dilog(ix)*(pisqos + logmxy*logmxy/2.) + dilog(y)*(logmxy*logmxy - logmx*logmx)/2.;
        }

        complex<double> li22smallainvstuffle(const complex<double> &x, const complex<double> &y)
        {
            const complex<double> xy = x*y;
            const complex<double> ix = 1./x;
            const complex<double> iy = 1./y;
            const complex<double> logmxy = std::log(-xy);
            const complex<double> logmx = std::log(-x);
            const double pisqos = 1.64493406684822643647;

            return li22smallastuffle(ix,iy) - quadlog(xy)+3.*(quadlog(y) + quadlog(ix))+2.*(trilog(ix) - trilog(y))*logmxy+
                dilog(ix)*(pisqos+logmxy*logmxy/2.) + dilog(y)*(logmxy*logmxy-logmx*logmx)/2.;
        }

        //This is Li22(x,x)
        complex<double> li22xx(complex<double> x)
        {
            const complex<double> xsq = x*x;
            const complex<double> li2c = dilog(x);

            return 0.5*(li2c*li2c - quadlog(xsq));
        }

        //This is Li22(1/y,y)
        const double li22iyyc1 = 9.86960440108935861883;  // Pi^2
        const double li22iyyc2 = -1.08232323371113819152;  // -Pi^4/90

        complex<double> li22iyy(complex<double> y)
        {
            const complex<double> logc = std::log(-y);
            const complex<double> li2c = dilog(y);

            return 3.*quadlog(y) - 0.5*li2c*(li2c + logc*logc + li22iyyc1) + li22iyyc2;
        }

        //This is Li22(x,1)
        const double li22x1c1 = -2.;
        const double li22x1c2 = 2.;
        const double li22x1c3 = 0.5;
        const double li22x1c4 = -0.166666666666666666667;
        const double li22x1c5 = 0.333333333333333333333;
        const double li22x1c6 = 1.64493406684822643647;  // Pi^2/6
        const double li22x1c7 = -2.40411380631918857080; // -2*Zeta(3)
        const double li22x1c8 = 2.16464646742227638303;  // Pi^4/45

        complex<double> li22x1(complex<double> x)
        {
            const complex<double> omx = 1.-x;
            const complex<double> logomx = std::log(omx);
            const complex<double> logomxsq = logomx*logomx;
            const complex<double> arg1 = 1./omx;
            const complex<double> arg2 = -x/omx;
            const complex<double> li2x = dilog(x);

            return li22x1c1*(quadlog(arg1) + quadlog(arg2) + quadlog(x)) + li22x1c2*trilog(x)*logomx +
                li22x1c3*li2x*li2x + logomxsq*( li22x1c4*logomxsq + li22x1c5*logomx*log(-x) + li22x1c6) +
                li22x1c7*logomx + li22x1c8;
        }

        //This is Li22(1,y)
        const double li221yc1 = 2.;
        const double li221yc2 = -2.;
        const double li221yc3 = -0.5;
        const double li221yc4 = 0.166666666666666666667;
        const double li221yc5 = -0.333333333333333333333;
        const double li221yc6 = -1.64493406684822643647; // -Pi^2/6
        const double li221yc7 = 2.40411380631918857080;  // 2*Zeta(3)
        const double li221yc8 = -2.16464646742227638303; // -Pi^4/45
        const double li221yc9 = 1.64493406684822643647;  // Pi^2/6

        complex<double> li221y(complex<double> y)
        {
            const complex<double> omy = 1.-y;
            const complex<double> logomy = std::log(omy);
            const complex<double> logomysq = logomy*logomy;
            const complex<double> arg1 = 1./omy;
            const complex<double> arg2 = -y/omy;
            const complex<double> li2y = dilog(y);

            return li221yc1*(quadlog(arg1) + quadlog(arg2)) + quadlog(y) + li221yc2*trilog(y)*logomy +
                (li221yc9 + li221yc3*li2y)*li2y + logomysq*( li221yc4*logomysq + li221yc5*logomy*std::log(-y) +
                li221yc6) + li221yc7*logomy + li221yc8;
        }


        const double epsdif = 5e-14;

        complex<double> li22basic(const complex<double> &x, const complex<double> &y)
        {
            const double absxy = abs(x*y);
            const double absx = abs(x);
            const double absy = abs(y);
            int ncmax;

            if(absxy < 0.7)
            {

                if(absxy<0.3){
                    ncmax=20;
                }
                else{if(absxy<0.5){
                        ncmax=35;
                    }
                    else{
                        ncmax=65;
                    }
                }

                if(absy>1.15||absx<0.25){
                    debug("A (Li22fast)");
                    return li22fast(x,y,ncmax);
                }
                else{
                    debug("B (Stuffle)");
                    return li22stuffle(x,y,ncmax);
                }

            }

            // The special formulae
            if(abs(y-x) < epsdif){
                debug("X (x=y)");
                return li22xx(x);
            }

            if(abs(y-1./x) < epsdif){
                debug("X (y=1/x)");
                return li22iyy(y);
            }

            if(abs(y-1.) < epsdif){
                debug("X (y=1)");
                return li22x1(x);
            }

            if(abs(x-1.) < epsdif){
                debug("X (x=1)");
                return li221y(y);
            }

            if(1.0 / absxy < 0.7)
            {

                if(absx<0.3){
                    debug("D/Ea (diagonal extra)");
                    return li22smalla(x,y);
                }

                if(absxy>3.33333){
                    ncmax=20;
                }
                else{if(absxy>2.0){
                        ncmax=35;
                    }
                    else{
                        ncmax=65;
                    }
                }

                if(absy<0.87||absx>3.5){
                    if(absy<1){
                        debug("C (Inversion (special))");
                        return li22invspecial(x,y,ncmax);
                    }
                    else{
                        debug("C (Inversion)");
                        return li22inv(x,y,ncmax);
                    }
                }
                else{
                    if(absx<1){
                        debug("D (InvStuffle (special))");
                        return li22invstufflespecial(x,y,ncmax);
                    }
                    else{
                        debug("D (InvStuffle)");
                        return li22invstuffle(x,y,ncmax);
                    }
                }
            }


            //  One of the points are on the bad circle
            if(absx<0.7){
                debug("Ea (diagonal)");
                return li22smalla(x,y);
            }

            if(absy<0.7){
                debug("Eb (diagonal+stuffle)");
                return li22smallastuffle(x,y);
            }

            if(1.0 / absx < 0.7){
                debug("Ec (diagonal+inversion)");
                return li22smallainv(x,y);
            }

            if(1.0 / absy < 0.7){
                debug("Ed (diagonal+inv+st)");
                return li22smallainvstuffle(x,y);
            }

            const double aax = abs(arg(x));
            const double aaxy = abs(arg(x*y));

            //   if(abs(x-y)<1e-14){
            //       debug("X (x=y (internal))");
            //       return li22xx((x+y)/2.);
            //   }

            // 1.04 = pi/3
            if(aax>1.04 && aaxy>1.04){
                debug("F (Holder)");
                return li22holder(x,y,100);
            }

            if(aaxy > 2.2){
                debug("Fb (Holder+stuffle)");
                return li22holderstuffle(x,y,100);
            }
            if(aax<1.15 && aaxy<1.15){
                debug("H (LogA1 center)");
                return li22logA1(x,y,40);
            }
            if(aax<0.25 || aaxy<0.25){
                debug("H (LogA1 edge)");
                return li22logA1(x,y,60);
            }
            debug("G (LogA0)");
            return li22logA0(x,y,50);
        }

        const double epsilon = 1e-14;
        const complex<double> imagi(0.,1.);
        /*
        PROBABLY NOT NEEDED SINCE POLYLOGS ALREADY IMPLEMENTED IN EOS
        complex<double> li(int n, const complex<double> &x)
        {
            complex<double> x2;
            if(abs(x.imag()) < epsilon && x.real() > 1.){
                x2 = x*(1.-5.*epsilon*imagi);
            }
            else{
                x2 = x;
            }
            return libasic(n,x2);
        }
        */

        void writecomplex(complex<double> x)
        {
            const double im = x.imag();
            char c;
            if(im<0){
                c='-';
            }
            else{
                c='+';
            }
            printf("%.20f%c%.20f*I\n",x.real(),c,abs(im));
        }

    }

    complex<double> li22(const complex<double> & x, const complex<double> & y)
    {

        complex<double> xy2, x2, y2;

        if(abs(x.imag()) < li22_impl::epsilon && abs(x-1.) > li22_impl::epsilon){
            x2 = x*(1.-5*li22_impl::epsilon*li22_impl::imagi);
        }
        else{
            x2 = x;
        }

        if(abs(y.imag()) < li22_impl::epsilon && abs(y-1.) > li22_impl::epsilon){
            y2 = y*(1.-5*li22_impl::epsilon*li22_impl::imagi);
        }
        else{
            y2 = y;
        }

        xy2 = x2*y2;

        if(abs(xy2.imag()) < li22_impl::epsilon && abs(x*y-1.) > li22_impl::epsilon){
            x2 = x2*(1.-5*li22_impl::epsilon*li22_impl::imagi);
        }

        return li22_impl::li22basic(x2,y2);
    }
}
