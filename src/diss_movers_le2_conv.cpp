#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"

template<typename T>
T Movers(T Fr, T Fl, T Ur, T Ul, T L_max, T L_min);

template<typename T>
T MaxEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny);

template<typename T>
T MinEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny);

template<typename T>
T Minmod_Slope(T del_p, T del_m, T dist);


void
Diss_MOVERS_LE2_Conv(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
                     double ***&sj, double ***&diss) {

    int im1, jm1, ip1, jp1;
    double nx, ny, ds, rhoinv, rhol, ul, vl, rhor, ur, vr, pl, pr, hl, hr, El, Er, gam1, ggm1, Vnl, Vnr;
    double ar, al, max_eig, min_eig, dx, dy;
    double *fd, *fr, *fl, *Ur, *Ul, *deltl, *deltr, *mdiss;

    //double area_avg;

    fd = new double[nconv];
    fr = new double[nconv];
    fl = new double[nconv];
    Ur = new double[nconv];
    Ul = new double[nconv];
    deltl = new double[nconv];
    deltr = new double[nconv];
    mdiss = new double[nconv];

    gam1 = Gamma - 1.0;
    ggm1 = Gamma / gam1;

    //compute in i-direction
    for (int j = 2; j <= jb; j++) {
        for (int i = 2; i <= id1; i++) {

            im1 = i - 1;
            ip1 = i + 1;

            ds = sqrt(si[0][i][j] * si[0][i][j] + si[1][i][j] * si[1][i][j]);
            nx = si[0][i][j] / ds;
            ny = si[1][i][j] / ds;
            dx = 0.5 * (x[i + 1][j] - x[i - 1][j]);

            //area_avg = 0.5*(area[im1][j]+area[i][j]);
            //left and right states

            rhoinv = 1.0 / cv[0][im1][j];
            rhol = cv[0][im1][j];
            ul = cv[1][im1][j] * rhoinv;
            vl = cv[2][im1][j] * rhoinv;
            pl = dv[0][im1][j];
            hl = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            El = pl / (gam1) + 0.5 * rhol * (ul * ul + vl * vl);
            al = dv[2][im1][j];
            Vnl = fabs(ul * nx + vl * ny);

            rhoinv = 1.0 / cv[0][i][j];
            rhor = cv[0][i][j];
            ur = cv[1][i][j] * rhoinv;
            vr = cv[2][i][j] * rhoinv;
            pr = dv[0][i][j];
            hr = ggm1 * pr / rhor + 0.5 * (ur * ur + vr * vr);
            Er = pr / (gam1) + 0.5 * rhor * (ur * ur + vr * vr);
            ar = dv[2][i][j];
            Vnr = fabs(ur * nx + vr * ny);

            max_eig = MaxEigVal(ur, ul, vr, vl, ar, al, nx, ny);
            min_eig = MinEigVal(ur, ul, vr, vl, ar, al, nx, ny);

            fl[0] = (rhol * ul) * nx + (rhol * vl) * ny;
            fl[1] = (rhol * ul * ul + pl) * nx + (rhol * ul * vl) * ny;
            fl[2] = (rhol * ul * vl) * nx + (rhol * vl * vl + pl) * ny;
            fl[3] = (rhol * ul * hl) * nx + (rhol * vl * hl) * ny;

            fr[0] = (rhor * ur) * nx + (rhor * vr) * ny;
            fr[1] = (rhor * ur * ur + pr) * nx + (rhor * ur * vr) * ny;
            fr[2] = (rhor * ur * vr) * nx + (rhor * vr * vr + pr) * ny;
            fr[3] = (rhor * ur * hr) * nx + (rhor * vr * hr) * ny;

            Ur[0] = rhor;
            Ur[1] = rhor * ur;
            Ur[2] = rhor * vr;
            Ur[3] = Er;

            Ul[0] = rhol;
            Ul[1] = rhol * ul;
            Ul[2] = rhol * vl;
            Ul[3] = El;

            mdiss[0] = Movers(fr[0], fl[0], Ur[0], Ul[0], max_eig, min_eig);
            mdiss[1] = Movers(fr[1], fl[1], Ur[1], Ul[1], max_eig, min_eig);
            mdiss[2] = Movers(fr[2], fl[2], Ur[2], Ul[2], max_eig, min_eig);
            mdiss[3] = Movers(fr[3], fl[3], Ur[3], Ul[3], max_eig, min_eig);

            deltl[0] = 0.5 * dx * Minmod_Slope(dui[0][i][j], dui[0][im1][j], dx);
            deltr[0] = 0.5 * dx * Minmod_Slope(dui[0][ip1][j], dui[0][i][j], dx);
            deltl[1] = 0.5 * dx * Minmod_Slope(dui[1][i][j], dui[1][im1][j], dx);
            deltr[1] = 0.5 * dx * Minmod_Slope(dui[1][ip1][j], dui[1][i][j], dx);
            deltl[2] = 0.5 * dx * Minmod_Slope(dui[2][i][j], dui[2][im1][j], dx);
            deltr[2] = 0.5 * dx * Minmod_Slope(dui[2][ip1][j], dui[2][i][j], dx);
            deltl[3] = 0.5 * dx * Minmod_Slope(dui[3][i][j], dui[3][im1][j], dx);
            deltr[3] = 0.5 * dx * Minmod_Slope(dui[3][ip1][j], dui[3][i][j], dx);

            Ul[0] = cv[0][im1][j] + deltl[0];
            Ul[1] = cv[1][im1][j] + deltl[1];
            Ul[2] = cv[2][im1][j] + deltl[2];
            Ul[3] = cv[3][im1][j] + deltl[3];

            Ur[0] = cv[0][i][j] - deltr[0];
            Ur[1] = cv[1][i][j] - deltr[1];
            Ur[2] = cv[2][i][j] - deltr[2];
            Ur[3] = cv[3][i][j] - deltr[3];

            fd[0] = 0.5 * fabs(mdiss[0] +
                               MinMod(dui[0][ip1][j], dui[0][i][j], dui[0][im1][j]) * (std::max(Vnr, Vnl) - mdiss[0])) *
                    (Ur[0] - Ul[0]);
            fd[1] = 0.5 * fabs(mdiss[1] +
                               MinMod(dui[1][ip1][j], dui[1][i][j], dui[1][im1][j]) * (std::max(Vnr, Vnl) - mdiss[1])) *
                    (Ur[1] - Ul[1]);
            fd[2] = 0.5 * fabs(mdiss[2] +
                               MinMod(dui[2][ip1][j], dui[2][i][j], dui[2][im1][j]) * (std::max(Vnr, Vnl) - mdiss[2])) *
                    (Ur[2] - Ul[2]);
            fd[3] = 0.5 * fabs(mdiss[3] +
                               MinMod(dui[3][ip1][j], dui[3][i][j], dui[3][im1][j]) * (std::max(Vnr, Vnl) - mdiss[3])) *
                    (Ur[3] - Ul[3]);

            //final dissipation terms

            diss[0][i][j] = diss[0][i][j] - fd[0] * ds;
            diss[1][i][j] = diss[1][i][j] - fd[1] * ds;
            diss[2][i][j] = diss[2][i][j] - fd[2] * ds;
            diss[3][i][j] = diss[3][i][j] - fd[3] * ds;

            diss[0][im1][j] = diss[0][im1][j] + fd[0] * ds;
            diss[1][im1][j] = diss[1][im1][j] + fd[1] * ds;
            diss[2][im1][j] = diss[2][im1][j] + fd[2] * ds;
            diss[3][im1][j] = diss[3][im1][j] + fd[3] * ds;
        }
    }

    //compute in j-direction
    for (int i = 2; i <= ib; i++) {
        for (int j = 2; j <= jd1; j++) {

            jm1 = j - 1;
            jp1 = j + 1;

            ds = sqrt(sj[0][i][j] * sj[0][i][j] + sj[1][i][j] * sj[1][i][j]);
            nx = sj[0][i][j] / ds;
            ny = sj[1][i][j] / ds;
            dy = 0.5 * (y[i][j + 1] - y[i][j - 1]);

            // area_avg = 0.5*(area[i][jm1]+area[i][j]);

            //left and right states


            rhoinv = 1.0 / cv[0][i][jm1];
            rhol = cv[0][i][jm1];
            ul = cv[1][i][jm1] * rhoinv;
            vl = cv[2][i][jm1] * rhoinv;
            pl = dv[0][i][jm1];
            hl = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            El = pl / (gam1) + 0.5 * rhol * (ul * ul + vl * vl);
            al = dv[2][i][jm1];
            Vnl = fabs(ul * nx + vl * ny);

            rhoinv = 1.0 / cv[0][i][j];
            rhor = cv[0][i][j];
            ur = cv[1][i][j] * rhoinv;
            vr = cv[2][i][j] * rhoinv;
            pr = dv[0][i][j];
            hr = ggm1 * pr / rhor + 0.5 * (ur * ur + vr * vr);
            Er = pr / (gam1) + 0.5 * rhor * (ur * ur + vr * vr);
            ar = dv[2][i][j];
            Vnr = fabs(ur * nx + vr * ny);

            fl[0] = (rhol * ul) * nx + (rhol * vl) * ny;
            fl[1] = (rhol * ul * ul + pl) * nx + (rhol * ul * vl) * ny;
            fl[2] = (rhol * ul * vl) * nx + (rhol * vl * vl + pl) * ny;
            fl[3] = (rhol * ul * hl) * nx + (rhol * vl * hl) * ny;

            fr[0] = (rhor * ur) * nx + (rhor * vr) * ny;
            fr[1] = (rhor * ur * ur + pr) * nx + (rhor * ur * vr) * ny;
            fr[2] = (rhor * ur * vr) * nx + (rhor * vr * vr + pr) * ny;
            fr[3] = (rhor * ur * hr) * nx + (rhor * vr * hr) * ny;

            Ur[0] = rhor;
            Ur[1] = rhor * ur;
            Ur[2] = rhor * vr;
            Ur[3] = Er;

            Ul[0] = rhol;
            Ul[1] = rhol * ul;
            Ul[2] = rhol * vl;
            Ul[3] = El;

            mdiss[0] = Movers(fr[0], fl[0], Ur[0], Ul[0], max_eig, min_eig);
            mdiss[1] = Movers(fr[1], fl[1], Ur[1], Ul[1], max_eig, min_eig);
            mdiss[2] = Movers(fr[2], fl[2], Ur[2], Ul[2], max_eig, min_eig);
            mdiss[3] = Movers(fr[3], fl[3], Ur[3], Ul[3], max_eig, min_eig);

            deltl[0] = 0.5 * dy * Minmod_Slope(duj[0][i][j], duj[0][i][jm1], dy);
            deltr[0] = 0.5 * dy * Minmod_Slope(duj[0][i][jp1], duj[0][i][j], dy);
            deltl[1] = 0.5 * dy * Minmod_Slope(duj[1][i][j], duj[1][i][jm1], dy);
            deltr[1] = 0.5 * dy * Minmod_Slope(duj[1][i][jp1], duj[1][i][j], dy);
            deltl[2] = 0.5 * dy * Minmod_Slope(duj[2][i][j], duj[2][i][jm1], dy);
            deltr[2] = 0.5 * dy * Minmod_Slope(duj[2][i][jp1], duj[2][i][j], dy);

            Ul[0] = cv[0][i][jm1] + deltl[0];
            Ul[1] = cv[1][i][jm1] + deltl[1];
            Ul[2] = cv[2][i][jm1] + deltl[2];
            Ul[3] = cv[3][i][jm1] + deltl[3];

            Ur[0] = cv[0][i][j] - deltr[0];
            Ur[1] = cv[1][i][j] - deltr[1];
            Ur[2] = cv[2][i][j] - deltr[2];
            Ur[3] = cv[3][i][j] - deltr[3];

            fd[0] = 0.5 * fabs(mdiss[0] +
                               MinMod(duj[0][i][jp1], duj[0][i][j], duj[0][im1][j]) * (std::max(Vnr, Vnl) - mdiss[0])) *
                    (Ur[0] - Ul[0]);
            fd[1] = 0.5 * fabs(mdiss[1] +
                               MinMod(duj[1][i][jp1], duj[1][i][j], duj[1][im1][j]) * (std::max(Vnr, Vnl) - mdiss[1])) *
                    (Ur[1] - Ul[1]);
            fd[2] = 0.5 * fabs(mdiss[2] +
                               MinMod(duj[2][i][jp1], duj[2][i][j], duj[2][im1][j]) * (std::max(Vnr, Vnl) - mdiss[2])) *
                    (Ur[2] - Ul[2]);
            fd[3] = 0.5 * fabs(mdiss[3] +
                               MinMod(duj[3][i][jp1], duj[3][i][j], duj[3][im1][j]) * (std::max(Vnr, Vnl) - mdiss[3])) *
                    (Ur[3] - Ul[3]);

            //final dissipation terms

            diss[0][i][j] = diss[0][i][j] - fd[0] * ds;
            diss[1][i][j] = diss[1][i][j] - fd[1] * ds;
            diss[2][i][j] = diss[2][i][j] - fd[2] * ds;
            diss[3][i][j] = diss[3][i][j] - fd[3] * ds;

            diss[0][i][jm1] = diss[0][i][jm1] + fd[0] * ds;
            diss[1][i][jm1] = diss[1][i][jm1] + fd[1] * ds;
            diss[2][i][jm1] = diss[2][i][jm1] + fd[2] * ds;
            diss[3][i][jm1] = diss[3][i][jm1] + fd[3] * ds;
        }
    }

    /*for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            std::cout<<"the computed dissipation values="<<diss[2][i][j]<<std::endl;
        }
    }*/

    delete[] fd;
    delete[] fr;
    delete[] fl;
    delete[] Ur;
    delete[] Ul;
    delete[] deltl;
    delete[] deltr;
    delete[] mdiss;
}


template<typename T>
T Movers(T Fr, T Fl, T Ur, T Ul, T L_max, T L_min) {

    double S;
    const double epsilon = 1e-4;

    if (fabs(Fr - Fl) < epsilon) return 0.0;
    else if (fabs(Ur - Ul) < epsilon) return L_min;
    else if (fabs(Ur - Ul) > epsilon and fabs(Fr - Fl) > epsilon) S = fabs(((Fr - Fl) / (Ur - Ul)));
    else S = L_min;

    if (S < epsilon) return 0;
    else if ((S) >= L_max) return (L_max);
    else if ((S) <= L_min) return (L_min);
    else return (S);
    return S;
}

template<typename T>
T MaxEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny) {
    double L1r, L2r, L3r, L1l, L2l, L3l;
    L1r = fabs(ur * nx + vr * ny + ar);
    L1l = fabs(ul * nx + vl * ny + al);
    L2r = fabs(ur * nx + vr * ny);
    L2l = fabs(ul * nx + vl * ny);
    L3r = fabs(ur * nx + vr * ny - ar);
    L3l = fabs(ul * nx + vl * ny - al);

    //return Max3(max(L1l,L1r), max(L2l,L2r), max(L3l,L3r));
    return Max2(Max3(L1r, L2r, L3r), Max3(L1l, L2l, L3l));
}

template<typename T>
T MinEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny) {
    double L1r, L2r, L3r, L1l, L2l, L3l;
    L1r = fabs(ur * nx + vr * ny + ar);
    L1l = fabs(ul * nx + vl * ny + al);
    L2r = fabs(ur * nx + vr * ny);
    L2l = fabs(ul * nx + vl * ny);
    L3r = fabs(ur * nx + vr * ny - ar);
    L3l = fabs(ul * nx + vl * ny - al);

    //return Max3(min(L1l,L1r), min(L2l,L2r), min(L3l,L3r));
    return Max2(Min3(L1r, L2r, L3r), Min3(L1l, L2l, L3l));
}

template<typename T>
T Minmod_Slope(T del_p, T del_m, T dist) {

    double r1, r2, phi;

    r1 = del_p / dist;
    r2 = del_m / dist;

    if (fabs(r1) < fabs(r2) and r1 * r2 > 0) phi = r1;
    else if (fabs(r2) < fabs(r1) and r1 * r2 > 0) phi = r2;
    else if (r1 * r2 <= 0) phi = 0.0;

    return phi;
}
