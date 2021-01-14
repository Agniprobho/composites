#include <iostream>

using namespace std;

#define E1f 221000000000;
#define E2f 13800000000; 
#define G12f 13800000000;
#define u12f 0.2;
#define u23f 0.25;
#define Em 4400000000;
#define PRm 0.34;
#define vf 0.7;
  
// Function to calculate and store inverse 
void inverse(double A[6][6], double inverse[6][6]) 
{ 
	int n=6;
	double d;
	/************ augmented matrix ***************/
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<n; j++)
			a[i][j] = A[i][j];
		for (int j=n; j<2*n; j++)
		{
			if (i==(j-n))
				a[i][j]=1.0;
			else
				a[i][j]=0.0;
		}
	}
	
	/************** partial pivoting **************/
    for(i=n-1;i>0;i--)
    {
        if(a[i-1][0]<a[i][0])
        for(j=0;j<n*2;j++)
        {
            d=a[i][j];
            a[i][j]=a[i-1][j];
            a[i-1][j]=d;
        }
    }
	
    /********** reducing to diagonal  matrix ***********/

    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        if(j!=i)
        {
            d=a[j][i]/a[i][i];
            for(k=0;k<n*2;k++)
                a[j][k]-=a[i][k]*d;
        }
    }
	
    /************** reducing to unit matrix *************/
    for(i=0;i<n;i++)
    {
        d=a[i][i];
        for(j=0;j<n*2;j++)
            a[i][j]=a[i][j]/d;
    }
	
	/************** inverted matrix *************/
	for(i=0;i<n;i++)
    {
        for(j=n;j<n*2;j++)
            inverse[i][j-n]=a[i][j];
    }
} 

void materialPropertiesPMM()
{ 
/*
	Periodic Microstructure Micromechanics (PMM) transversely isotropic fiber.
	Elastic-only equations taken from: Barbero, E.J. and Luciano, R., "Micromechanical formulas for the
	relaxation tensor of linear viscoelastic composites with transversely isotropic fibers", International Journal of Solids and Structures,
	32(13):1859--1872, 1995. http:%barbero.cadec-online.com/papers/1995/95BarberoLucianoMicromechanicalFormulas.pdf
*/
	double EA = E1f;
	double ET = E2f; 
	double GA = G12f;
	double vA = u12f;
	double vT = u23f;
	double EM = Em;
	double vM = PRm;
	double Vf = vf;
	
	double f, Delta;
	double C_11, C_22, C_33, C_12, C_13, C_23, C_44, C_55, C_66;
	double C11t, C12t, C22t, C23t, C44t, C66t;
	double C11, C22, C33, C12, C13, C23, C44, C55, C66;
	double lam_m, mu_m; //mu_m is Gm
	double S3, S6, S7;
	double a1, a2, a3, a4;
	
	//fiber coefficients
    f =ET/EA;
    Delta = (1-2*(vA*vA)*f-(vT*vT)-2*(vA*vA)*vT*f)/(EA*(ET*ET));
    C_11 = (1-(vT*vT))/((ET*ET)*Delta);
    C_22 = (1-(vA*vA)*f)/(EA*ET*Delta);
    C_33 = C_22;
    C_12 = (vA*f+vA*vT*f)/((ET*ET)*Delta);
    C_13 = C_12;
    C_23 = (vT + (vA*vA)*f)/(EA*ET*Delta);
    C_44 = ET/(2*(1+vT));
    C_55 = GA;
    C_66 = C_55;
	
	//matrix coefficients
    lam_m = (EM*vM)/((1+vM)*(1-2*vM));
    mu_m = EM/(2*(1+vM)); //isotropic matrix
	
    //geometry
    S3 = 0.49247 - 0.47603*Vf - 0.02748*(Vf*Vf);
    S6 = 0.36844 - 0.14944*Vf - 0.27152*(Vf*Vf);
    S7 = 0.12346 - 0.32035*Vf + 0.23517*(Vf*Vf);
	
    //a_values
    a1 = 4*(mu_m*mu_m) - 2*mu_m*C_33 + 6*lam_m*mu_m - 2*C_11*mu_m - 2*mu_m*C_23 + C_23*C_11 + 4*lam_m*C_12 - 2*(C_12*C_12) - lam_m*C_33 -2*C_11*lam_m + C_11*C_33 - lam_m*C_23;
    a2 = 8*(mu_m*mu_m*mu_m)- 8*(mu_m*mu_m)*C_33 + 12*(mu_m*mu_m)*lam_m -4*(mu_m*mu_m)*C_11 - 2*mu_m*(C_23*C_23) + 4*mu_m*lam_m*C_23 + 4*mu_m*C_11*C_33 - 8*mu_m*lam_m*C_33 - 4*mu_m*(C_12*C_12) + 2*mu_m*(C_33*C_33) -4*mu_m*C_11*lam_m + 8*mu_m*lam_m*C_12 + 2*lam_m*C_11*C_33 + 4*C_12*C_23*lam_m - 4*C_12*C_33*lam_m - 2*lam_m*C_11*C_23 - 2*C_23*(C_12*C_12) + (C_23*C_23)*C_11 + 2*C_33*(C_12*C_12) - C_11*(C_33*C_33) + lam_m*(C_33*C_33) - lam_m*(C_23*C_23); 
    a3 = ((4*(mu_m*mu_m) + 4*lam_m*mu_m - 2*C_11*mu_m - 2*mu_m*C_33 - C_11*lam_m - lam_m*C_33 - (C_12*C_12))/a2)  + ((C_11*C_33 + 2*lam_m*C_12)/a2) - ((S3-((S6)/(2-2*vM)))/mu_m) ;
    a4 = -1*((-2*mu_m*C_23 + 2*lam_m*mu_m - lam_m*C_23 - C_11*lam_m - (C_12*C_12) + 2*lam_m*C_12 + C_11*C_23)/a2) + (S7)/(mu_m*(2-2*vM));
    
	//C_values
    C11t = lam_m + 2*mu_m - Vf*(-(a4*a4) + (a3*a3))*(1/(-1*(((2*mu_m + 2*lam_m - C_33 - C_23)*((a4*a4)-(a3*a3)))/a1) + ((2*(a4-a3)*(lam_m-C_12)*(lam_m-C_12))/(a1*a1)))); 
    C12t = lam_m + Vf*(((lam_m-C_12)*(a4-a3))/a1)*(1/(1*(((2*mu_m + 2*lam_m - C_33 - C_23)*((a3*a3)-(a4*a4)))/a1) + ((2*(a4-a3)*(lam_m-C_12)*(lam_m-C_12))/(a1*a1)))); 
    C22t = lam_m + 2*mu_m - Vf*(((2*mu_m + 2*lam_m - C_33 - C_23)*a3/a1) -(((lam_m-C_12)*(lam_m-C_12))/(a1*a1)))*(1/(1*(((2*mu_m + 2*lam_m - C_33 - C_23)*((a3*a3)-(a4*a4)))/a1) + ((2*(a4-a3)*(lam_m-C_12)*(lam_m-C_12))/(a1*a1))));
    C23t = lam_m + Vf*(((2*mu_m + 2*lam_m - C_33 - C_23)*a4/a1) -(((lam_m-C_12)*(lam_m-C_12))/(a1*a1)))*(1/(1*(((2*mu_m + 2*lam_m - C_33 - C_23)*((a3*a3)-(a4*a4)))/a1) + ((2*(a4-a3)*(lam_m-C_12)*(lam_m-C_12))/(a1*a1))));
    C44t = mu_m - Vf*(1/((2/(2*mu_m - C_22 + C_23)) - (1/(mu_m))*(2*S3 - (4*S7/(2-2*vM)))));
    C66t = mu_m -Vf*(1/((1/(mu_m-C_66)) - S3/mu_m));
    
	//C_total
    C11 = C11t;
    C12 = C12t;
    C13 = C12t;
    C22 = (3.0/4.0)*C22t + (1.0/4.0)*C23t +(1.0/4.0)*C44t;
    C33 = C22;
    C23 = (1.0/4.0)*C22t + (3.0/4.0)*C23t -(1.0/4.0)*C44t;
    C55 = C66t;
    C66 = C66t;
    C44 = (C22-C23)/2.0;
	
	//properties
    double C[6][6] = {{C11, C12, C13, 0, 0, 0}, {C12, C22, C23, 0, 0, 0}, {C13, C23, C33, 0, 0, 0}, {0, 0, 0, C44, 0, 0}, {0, 0, 0, 0, C55, 0}, {0, 0, 0, 0, 0, C66}};
	
	double S[6][6];
    inverse(C,S);
    
    for (int i=0; i<6; i++)
	{
	    for (int j=0; j<6; j++)
	        cout<<S[i][j]<<"\t";
	   cout<<"\n\n";
	}
	
	double EXy, EYy, EZy, GXYy, GYZy, GXZy, PRXYy, PRYZy, PRXZy;
	EXy = 1/S[0][0]; //E1
	EYy = 1/S[1][1]; //E2
	EZy = EYy; //E3 = E2
	GXYy = 1/S[5][5]; //G12
	GYZy = 1/S[3][3]; //G23
	GXZy = GXYy; //G13 = G12
	PRXYy = -S[1][0]/S[0][0]; //v12
	PRYZy = -S[2][1]/S[1][1]; //v23
	PRXZy = PRXYy; //v13 = v12
	
	cout<<"E1\t"<<EXy<<"\n"<<"E2\t"<<EYy<<"\n"<<"v12\t"<<PRXYy<<"\n"<<"v23\t"<<PRYZy<<"\n"<<"G12\t"<<GXYy<<"\n"<<"G23\t"<<GYZy<<"\n";
}

int main()
{
	materialPropertiesPMM();
	return 0;
}