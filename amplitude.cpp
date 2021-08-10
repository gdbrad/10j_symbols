#include <iostream>
#include <cmath>
using namespace std;

class DFKR{
    public:
    double j_f[5] = {0,0.5,1,1.5,2};
    int f;
    int e;
    int v;
    
    _Complex int A_f;
    _Complex int A_e;
    _Complex int A_v;

    double j_e[4];      //  spins labelling 4 faces incident to each edge    //
    double j_v[10];     //  spins labelling 10 faces incident to each vertex //

    // construct adjacency matrix //
    
    const int adjacencies[5][4] = {
    { 1, 2, 3, 4},
    { 0, 2, 3, 4},
    { 0, 1, 3, 4},
    { 0, 1, 2, 4},
    { 0, 1, 2, 3}

    };
    // gluings invovled in the 4-sphere triangulation along simplex boundary //
    const int gluings[6][5][5] = {
    { { 0, 1, 2, 3, 4}, {1,0,2,3,4}{}}    
    { { 0, 1, 2, 3 }, { 1, 0, 2, 3 }, { 1, 2, 0, 3 }, { 1, 2, 3, 0 } },
    { { 0, 1, 2, 3 }, { 0, 1, 2, 3 }, { 0, 2, 1, 3 }, { 0, 2, 3, 1 } },
    { { 1, 0, 2, 3 }, { 0, 1, 2, 3 }, { 0, 1, 2, 3 }, { 0, 1, 3, 2 } },
    { { 2, 0, 1, 3 }, { 0, 2, 1, 3 }, { 0, 1, 2, 3 }, { 0, 1, 2, 3 } },
    { { 3, 0, 1, 2 }, { 0, 3, 1, 2 }, { 0, 1, 3, 2 }, { 0, 1, 2, 3 } }
    _Complex int Z_F (_Complex int A_f, _Complex int A_e, _Complex int A_v){
        // use 1 j_f ; product over 2-simplices -> face amplitude //
        // 4j symbol //
        for(f=0;f<i;f++){
            return (2j(f) + 1) ** 2;

        }
        // use 4 spins j_e ; product over 3-simplices -> edge amplitude //
        for(e=0;e<i;e++){

        }
        // use 10 spins j_v ; product over 4-simplices -> vertex amplitude //
        for(v=0;v<i;v++){

        }

    };

};
