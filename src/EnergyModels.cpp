#include "EnergyModels.hpp"
#include <assert.h>

double GCWeight = 1.;

double TEMP = 37.;

double COMPATIBLE_BP_MATRIX[4][4] =
        /* 3' NT:     N_A,     N_U,     N_G,     N_C */
        {
        /*5' NT: N_A*/
               {  DBL_MAX,      0., DBL_MAX, DBL_MAX},
        /*5' NT: N_U*/
               {       0., DBL_MAX,      0., DBL_MAX},
        /*5' NT: N_G*/
               {  DBL_MAX,      0., DBL_MAX,      0.},
        /*5' NT: N_C*/
               {  DBL_MAX, DBL_MAX,      0., DBL_MAX}
        };

double NUSSINOV_BP_MATRIX[4][4] =
        /* 3' NT:     N_A,     N_U,     N_G,     N_C */
        {
        /*5' NT: N_A*/
               {  DBL_MAX,     -2., DBL_MAX, DBL_MAX},
        /*5' NT: N_U*/
               {      -2., DBL_MAX,     -1., DBL_MAX},
        /*5' NT: N_G*/
               {  DBL_MAX,     -1., DBL_MAX,     -3.},
        /*5' NT: N_C*/
               {  DBL_MAX, DBL_MAX,     -3., DBL_MAX}
        };

double GC_TERM = -0.09070;
double AU_TERM = 1.26630;
double GU_TERM = 0.78566;
double GC_IN = -2.10208;
double AU_IN = -0.52309;
double GU_IN = -0.88474;

double FITTING_BP_INNER_MATRIX[4][4] =
        /* 3' NT:     N_A,     N_U,     N_G,     N_C */
        {
        /*5' NT: N_A*/
               {  DBL_MAX,   AU_IN, DBL_MAX, DBL_MAX},
        /*5' NT: N_U*/
               {    AU_IN, DBL_MAX,   GU_IN, DBL_MAX},
        /*5' NT: N_G*/
               {  DBL_MAX,   GU_IN, DBL_MAX,   GC_IN},
        /*5' NT: N_C*/
               {  DBL_MAX, DBL_MAX,   GC_IN, DBL_MAX}
        };

double FITTING_BP_TERMINAL_MATRIX[4][4] =
        /* 3' NT:     N_A,     N_U,     N_G,     N_C */
        {
        /*5' NT: N_A*/
               {  DBL_MAX, AU_TERM, DBL_MAX, DBL_MAX},
        /*5' NT: N_U*/
               {  AU_TERM, DBL_MAX, GU_TERM, DBL_MAX},
        /*5' NT: N_G*/
               {  DBL_MAX, GU_TERM, DBL_MAX, GC_TERM},
        /*5' NT: N_C*/
               {  DBL_MAX, DBL_MAX, GC_TERM, DBL_MAX}
        };

double STACK_AUAU = -0.18826;
double STACK_AUCG = -1.13291;
double STACK_AUGC = -1.09787;
double STACK_AUGU = -0.38606;
double STACK_AUUA = -0.26510;
double STACK_AUUG = -0.62086;
double STACK_CGAU = -1.11752;
double STACK_CGCG = -2.23740;
double STACK_CGGC = -1.89434;
double STACK_CGGU = -1.22942;
double STACK_CGUA = -1.10548;
double STACK_CGUG = -1.44085;
double STACK_GUAU = -0.55066;
double STACK_GUCG = -1.26209;
double STACK_GUGC = -1.58478;
double STACK_GUGU = -0.72185;
double STACK_GUUA = -0.49625;
double STACK_GUUG = -0.68876;

double FITTING_STACKS_MATRIX[4][4][4][4];

void initStackingModel()
{
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            for(int k=0;k<4;k++){
                for(int l=0;l<4;l++){
                    FITTING_STACKS_MATRIX[(Nucleotide) i][(Nucleotide) j][(Nucleotide) k][(Nucleotide) l] = DBL_MAX;
                }
            }
        }
    }
    FITTING_STACKS_MATRIX[N_A][N_U][N_A][N_U] = STACK_AUAU;
    FITTING_STACKS_MATRIX[N_U][N_A][N_A][N_U] = STACK_AUAU;
    FITTING_STACKS_MATRIX[N_A][N_U][N_C][N_G] = STACK_AUCG;
    FITTING_STACKS_MATRIX[N_U][N_A][N_C][N_G] = STACK_AUCG;
    FITTING_STACKS_MATRIX[N_A][N_U][N_G][N_C] = STACK_AUGC;
    FITTING_STACKS_MATRIX[N_U][N_A][N_G][N_C] = STACK_AUGC;
    FITTING_STACKS_MATRIX[N_A][N_U][N_G][N_U] = STACK_AUGU;
    FITTING_STACKS_MATRIX[N_U][N_A][N_G][N_U] = STACK_AUGU;
    FITTING_STACKS_MATRIX[N_A][N_U][N_U][N_A] = STACK_AUUA;
    FITTING_STACKS_MATRIX[N_U][N_A][N_U][N_A] = STACK_AUUA;
    FITTING_STACKS_MATRIX[N_A][N_U][N_U][N_G] = STACK_AUUG;
    FITTING_STACKS_MATRIX[N_U][N_A][N_U][N_G] = STACK_AUUG;

    FITTING_STACKS_MATRIX[N_C][N_G][N_A][N_U] = STACK_CGAU;
    FITTING_STACKS_MATRIX[N_G][N_C][N_A][N_U] = STACK_CGAU;
    FITTING_STACKS_MATRIX[N_C][N_G][N_C][N_G] = STACK_CGCG;
    FITTING_STACKS_MATRIX[N_G][N_C][N_C][N_G] = STACK_CGCG;
    FITTING_STACKS_MATRIX[N_C][N_G][N_G][N_C] = STACK_CGGC;
    FITTING_STACKS_MATRIX[N_G][N_C][N_G][N_C] = STACK_CGGC;
    FITTING_STACKS_MATRIX[N_C][N_G][N_G][N_U] = STACK_CGGU;
    FITTING_STACKS_MATRIX[N_G][N_C][N_G][N_U] = STACK_CGGU;
    FITTING_STACKS_MATRIX[N_C][N_G][N_U][N_A] = STACK_CGUA;
    FITTING_STACKS_MATRIX[N_G][N_C][N_U][N_A] = STACK_CGUA;
    FITTING_STACKS_MATRIX[N_C][N_G][N_U][N_G] = STACK_CGUG;
    FITTING_STACKS_MATRIX[N_G][N_C][N_U][N_G] = STACK_CGUG;

    FITTING_STACKS_MATRIX[N_G][N_U][N_A][N_U] = STACK_GUAU;
    FITTING_STACKS_MATRIX[N_U][N_G][N_A][N_U] = STACK_GUAU;
    FITTING_STACKS_MATRIX[N_G][N_U][N_C][N_G] = STACK_GUCG;
    FITTING_STACKS_MATRIX[N_U][N_G][N_C][N_G] = STACK_GUCG;
    FITTING_STACKS_MATRIX[N_G][N_U][N_G][N_C] = STACK_GUGC;
    FITTING_STACKS_MATRIX[N_U][N_G][N_G][N_C] = STACK_GUGC;
    FITTING_STACKS_MATRIX[N_G][N_U][N_G][N_U] = STACK_GUGU;
    FITTING_STACKS_MATRIX[N_U][N_G][N_G][N_U] = STACK_GUGU;
    FITTING_STACKS_MATRIX[N_G][N_U][N_U][N_A] = STACK_GUUA;
    FITTING_STACKS_MATRIX[N_U][N_G][N_U][N_A] = STACK_GUUA;
    FITTING_STACKS_MATRIX[N_G][N_U][N_U][N_G] = STACK_GUUG;
    FITTING_STACKS_MATRIX[N_U][N_G][N_U][N_G] = STACK_GUUG;

}

bool isCompatible(Nucleotide n1, Nucleotide n2){
    return COMPATIBLE_BP_MATRIX[n1][n2]!=DBL_MAX;
}

double dGBasePair(Nucleotide n1, Nucleotide n2, bool isTerminal){
    switch (dGModel) {
    case COMPATIBLE_BP_MODEL:
        return COMPATIBLE_BP_MATRIX[n1][n2];
    case NUSSINOV_BP_MODEL:
        return NUSSINOV_BP_MATRIX[n1][n2];
    case FITTED_BP_MODEL:
        if (isTerminal)
            return FITTING_BP_TERMINAL_MATRIX[n1][n2];
        else
            return FITTING_BP_INNER_MATRIX[n1][n2];
    case FITTED_STACKING_PAIRS_MODEL:
        break;
    }
    assert(false && "Energy model not implemented yet!");
    return 0.;
}

double dGStacking(Nucleotide n5A, Nucleotide n5B,Nucleotide n3B, Nucleotide n3A){
    switch (dGModel) {
    case FITTED_STACKING_PAIRS_MODEL:
        return FITTING_STACKS_MATRIX[n5A][n3A][n5B][n3B];
    }
    assert(false && "Energy model not implemented yet!");
    return 0.;
}


double dGStructure(SecondaryStructure * r, const string &seq)
{
    vector<Loop *> loops = r->getLoops();
    double res = 0.;
    for (int i=0;i<loops.size();i++)
    {
       vector<int> indices = loops[i]->getIndices();
       vector<Nucleotide> nts;
       for (int j=0;j<indices.size();j++)
       {
           nts.push_back(char2nt(seq[indices[j]]));
       }
       res += loops[i]->scoreLoop(nts);
    }
    return res;
}


void updateFittedModel(double N_GC_IN, double N_AU_IN, double N_GU_IN, double N_GC_TERM, double N_AU_TERM, double N_GU_TERM){
    GC_IN = N_GC_IN;
    AU_IN = N_AU_IN;
    GU_IN = N_GU_IN;
    GC_TERM = N_GC_TERM;
    AU_TERM = N_AU_TERM;
    GU_TERM = N_GU_TERM;

    FITTING_BP_INNER_MATRIX[N_G][N_C] = GC_IN;
    FITTING_BP_INNER_MATRIX[N_C][N_G] = GC_IN;
    FITTING_BP_INNER_MATRIX[N_A][N_U] = AU_IN;
    FITTING_BP_INNER_MATRIX[N_U][N_A] = AU_IN;
    FITTING_BP_INNER_MATRIX[N_G][N_U] = GU_IN;
    FITTING_BP_INNER_MATRIX[N_U][N_G] = GU_IN;

    FITTING_BP_TERMINAL_MATRIX[N_G][N_C] = GC_TERM;
    FITTING_BP_TERMINAL_MATRIX[N_C][N_G] = GC_TERM;
    FITTING_BP_TERMINAL_MATRIX[N_A][N_U] = AU_TERM;
    FITTING_BP_TERMINAL_MATRIX[N_U][N_A] = AU_TERM;
    FITTING_BP_TERMINAL_MATRIX[N_G][N_U] = GU_TERM;
    FITTING_BP_TERMINAL_MATRIX[N_U][N_G] = GU_TERM;
}


