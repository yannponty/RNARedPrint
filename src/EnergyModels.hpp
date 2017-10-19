#ifndef ENERGYMODELS_H
#define ENERGYMODELS_H

#include <float.h>
#include "Nucleotide.hpp"
#include "RNAStructure.hpp"

using namespace std;

extern double TEMP;
#define RT (0.0019872370936902486 * (273.15 + TEMP) * 100)

extern double GCWeight;

extern double COMPATIBLE_BP_MATRIX[4][4];

extern double NUSSINOV_BP_MATRIX[4][4];

extern double GC_IN, AU_IN, GU_IN, GC_TERM, AU_TERM, GU_TERM;

extern double FITTING_BP_INNER_MATRIX[4][4];
extern double FITTING_BP_TERMINAL_MATRIX[4][4];

extern double FITTING_STACKS_MATRIX[4][4][4][4];


typedef enum {COMPATIBLE_BP_MODEL, NUSSINOV_BP_MODEL, FITTED_BP_MODEL, FITTED_STACKING_PAIRS_MODEL}  EnergyModel;

extern EnergyModel dGModel;

bool isCompatible(Nucleotide n1, Nucleotide n2);

double dGBasePair(Nucleotide n1, Nucleotide n2, bool isTerminal);

double dGStacking(Nucleotide n5A, Nucleotide n5B,Nucleotide n3B, Nucleotide n3A);

double dGStructure(SecondaryStructure * r, const string & seq);

void updateFittedModel(double N_GC_IN, double N_AU_IN, double N_GU_IN, double N_GC_TERM, double N_AU_TERM, double N_GU_TERM);

void initStackingModel();

#endif // ENERGYMODELS_H
