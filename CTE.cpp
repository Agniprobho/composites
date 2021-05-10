#include <iostream>

using namespace std;

#define E1f 220000000000;
#define E2f 14000000000; 
#define G12f 9000000000;
#define u12f 0.2;
#define u23f 0.25;
#define Em 3450000000;
#define PRm 0.35;
#define vf 0.712098;
#define alph_A -0.99;
#define alph_T 16.7;
#define alph_M 50;

void CTE()
{ 
/*
	Most accurate formula for predicting lamina coefficientsof thermal expansion (CTE) was derived by Levin
	Reference: Levin, V. M. Mekhanika Tverdogo Tela ., 88 (1968)
*/
	// Compute homogenized properties using Chamis' model
	double E1 = E1f * vf + Em * (1 - vf);
	double E2 = Em / (1 - (sqrt(vf) * (1 - (Em / E2f))));
	double E3 = E2; //transverse isotropy
	double G12 = Gm / (1 - (sqrt(vf) * (1 - (Gm / G12f))));
	double G13 = G12; //transverse isotropy
	double G23 = Gm / (1 - (sqrt(vf) * (1 - (Gm / G23f))));
	double v12 = u12f * vf + PRm * (1 - vf);
	double v13 = v12; //transverse isotropy
	double v23 = E2 / (2 * G23) - 1;
	
	// CTE Levin equation constants
	// Compliance coefficients
	double S11 = 1 / E1;
	double S22 = 1 / E2;
	double S12 = -v12 / E1;
	double S23 = -v23 / E1;
	
	// Weighted average values
	double a1_hat = alph_A * vf + alph_M * (1 - vf);
	double a2_hat = alph_T * vf + alph_M * (1 - vf);
	double S11_hat = (vf / E1f) + ((1 - vf) / Em);
	double S22_hat = (vf / E2f) + ((1 - vf) / Em);
	double S12_hat = -((u12f / E1f) * vf) - ((PRm / Em) * (1 - vf));
	double S23_hat = -((u23f / E2f) * vf) - ((PRm / Em) * (1 - vf));
	
	// P-coefficients
	double A11 = (1 / E1f) - (1 / Em);
	double A22 = (1 / E2f) - (1 / Em);
	double A12 = -(u12f / E1f) + (PRm / Em);
	double A23 = -(u23f / E2f) + (PRm / Em);
	
	double detA = A11 * (A22 * A22 - A23 * A23) + 2 * A12 * (A12 * A23 - A22 * A12);
	
	double P11 = (A22 * A22 - A23 * A23) / detA;
	double P22 = (A11 * A22 - A12 * A12) / detA;
	double P12 = (A12 * A23 - A22 * A12) / detA;
	double P23 = (A12 * A12 - A11 * A23) / detA;
	
	// Coefficients in the longitudinal and transverse directions
	double alph_1 = a1_hat + (S11 - S11_hat) * ((alph_A - alph_M) * P11 + 2 * (alph_T - alph_M) * P12) + 
					2 * (S12 - S12_hat) * ((alph_A - alph_M) * P12 + (alph_T - alph_M) * (P22 + P23));
	double alph_2 = a2_hat + (S12 - S12_hat) * ((alph_A - alph_M) * P11 + 2 * (alph_T - alph_M) * P12) + 
					(S22 - S22_hat + S23 - S23_hat) * ((alph_A - alph_M) * P12 + (alph_T - alph_M) * (P22 + P23));
	
	
	cout<<"\nalpha1\t"<<alph_1<<"\n"<<"alpha2\t"<<alph_2<<"\n";
}

int main()
{
	CTE();
	return 0;
}