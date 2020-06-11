#include <Rmath.h>
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

int DFS_cycle(const mat &U, int start, int end)
{
    if(start == end)
        return 1;

    int p = (int) A.n_cols;
    int bottom = 0;
    int top = 0;
    int S_len = 1;
    int cycle = 0;

    uvec AN(p);
    AN.fill(0);
    uvec S(p);
    S.fill(0);

    AN(start) = 1;
    S(0) = start;

    while(S_len > 0)
    {
        int i = S(bottom);
        S_len--;
        bottom++;

        for(int j = 0; j < p; ++j)
        {
            if(A(i,j) != 0)
            {
                if(j == end)
                {
                    cycle = 1;
                    break;
                }
                else
                {
                    if(AN(j) == 0)
                    {
                        top++;
                        S(top) = j;
                        S_len++;
                        AN(j) = 1;
                    }
                }
            }
        }

        if(cycle == 1)
            break;
    }

    return cycle;
}


// [[Rcpp::export]]
List intdag_ans(mat V, mat D, mat CorrMat) 
{
    int p = (int) V.n_cols;
    int q = (int) V.n_rows;

    uvec x_idx(q);
    for(int l = 0; l < q; ++l)
        x_idx(l) = l;
    uvec y_idx(p);
    for(int j = 0; j < p; ++j)
        y_idx(j) = j;

    mat Pi = zeros(p,p);
    mat Phi = zeros(q,p);

    mat V_nz(q,p);
    for(int l = 0; l < q; ++l) 
        for(int j = 0; j < p; ++j) 
            V.nz(l,j) = (V(l,j) == 0, 0, 1);
    
    L0_col = 










    

    return List::create(Named("Pi") = Pi, Named("Phi") = Phi);
}





