/**********
    R Interface.

    Copyright (C) 2020-2021  Chunlin Li

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <https://www.gnu.org/licenses/>.
**********/

#include "R.h"
#include "Rinternals.h"
#include "glmtlp.h"

extern "C"
{

    SEXP rlasso(SEXP b0,
                SEXP b,
                SEXP r,
                SEXP X,
                SEXP w,
                SEXP rho,
                SEXP lambda,
                SEXP nlambda,
                SEXP n,
                SEXP p,
                SEXP delta,
                SEXP tol,
                SEXP cd_maxit)
    {
        linreg_l1_ssr(REAL(b0),
                      REAL(b),
                      REAL(r),
                      REAL(X),
                      REAL(w),
                      REAL(rho),
                      REAL(lambda),
                      *INTEGER(nlambda),
                      *INTEGER(n),
                      *INTEGER(p),
                      *REAL(delta),
                      *REAL(tol),
                      *INTEGER(cd_maxit));

        return R_NilValue;
    }

    SEXP rtlp(SEXP b0,
              SEXP b,
              SEXP r,
              SEXP X,
              SEXP w,
              SEXP rho,
              SEXP lambda,
              SEXP nlambda,
              SEXP tau,
              SEXP n,
              SEXP p,
              SEXP delta,
              SEXP tol,
              SEXP cd_maxit,
              SEXP dc_maxit)
    {
        linreg_tlp_ssr(REAL(b0),
                       REAL(b),
                       REAL(r),
                       REAL(X),
                       REAL(w),
                       REAL(rho),
                       REAL(lambda),
                       *INTEGER(nlambda),
                       *REAL(tau),
                       *INTEGER(n),
                       *INTEGER(p),
                       *REAL(delta),
                       *REAL(tol),
                       *INTEGER(cd_maxit),
                       *INTEGER(dc_maxit));

        return R_NilValue;
    }

    SEXP rl0(SEXP b0,
             SEXP b,
             SEXP r,
             SEXP X,
             SEXP w,
             SEXP rho,
             SEXP s,
             SEXP ns,
             SEXP lambda,
             SEXP nlambda,
             SEXP tau,
             SEXP n,
             SEXP p,
             SEXP delta,
             SEXP tol,
             SEXP cd_maxit,
             SEXP dc_maxit)
    {
        linreg_l0_ssr(REAL(b0),
                      REAL(b),
                      REAL(r),
                      REAL(X),
                      REAL(w),
                      REAL(rho),
                      INTEGER(s),
                      *INTEGER(ns),
                      REAL(lambda),
                      *INTEGER(nlambda),
                      *REAL(tau),
                      *INTEGER(n),
                      *INTEGER(p),
                      *REAL(delta),
                      *REAL(tol),
                      *INTEGER(cd_maxit),
                      *INTEGER(dc_maxit));

        return R_NilValue;
    }

    SEXP rlogistic_l1(SEXP b0,
                      SEXP b,
                      SEXP y,
                      SEXP X,
                      SEXP rho,
                      SEXP lambda,
                      SEXP nlambda,
                      SEXP n,
                      SEXP p,
                      SEXP delta,
                      SEXP tol,
                      SEXP nr_maxit,
                      SEXP cd_maxit)
    {
        logistic_l1_ssr(REAL(b0),
                        REAL(b),
                        REAL(y),
                        REAL(X),
                        REAL(rho),
                        REAL(lambda),
                        *INTEGER(nlambda),
                        *INTEGER(n),
                        *INTEGER(p),
                        *REAL(delta),
                        *REAL(tol),
                        *INTEGER(nr_maxit),
                        *INTEGER(cd_maxit));

        return R_NilValue;
    }

} // extern "C"
