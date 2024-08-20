#include <Rcpp.h>
using namespace Rcpp; 

// [[Rcpp::export]]
NumericMatrix one_hot(CharacterVector& x){
    const int N_CHAR = x[0].size();
    const int N_STRINGS = x.size();
    const int N_BASES = 4;
    // initialized to 0
    NumericMatrix out(N_STRINGS, N_CHAR * N_BASES);
    int col_offset = 0;

    for (int i = 0; i < N_STRINGS; i++) {
        for(int j = 0; j < N_CHAR; j++) {
            if(x[i][j] == 'A') {
                col_offset = 0;
            } else if(x[i][j] == 'T') {
                col_offset = 1;
            } else if(x[i][j] == 'G') {
                col_offset = 2;
            } else {
                col_offset = 3;
            }

            out(i, j * N_BASES + col_offset) = 1;
        }
    }
    return out;
}
