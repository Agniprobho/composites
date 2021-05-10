#include <iostream>

using namespace std;

#define E1f 87000000000;
#define E2f 87000000000;
#define G12f 35300000000;
#define G23f 35300000000;
#define u12f 0.23;
#define u23f 0.23;
#define F1f 4800000000;
#define F2f 3500000000;
#define FSf 1100000000;
#define Em 2490000000;
#define Gm 920000000;
#define PRm 0.35;
#define Ftm 57000000;
#define Fcm 77000000;
#define Fsm 45000000;
#define vf 0.7;
  
void materialPropertiesChamis()
{ 
/*
	Chamis' Microstructure model to calculate material properties and strength.
*/
	double EXy = E1f * vf + Em * (1 - vf);
	double EYy = Em / (1 - (sqrt(vf) * (1 - (Em / E2f))));
	double EZy = EYy; //transverse isotropy
	double GXYy = Gm / (1 - (sqrt(vf) * (1 - (Gm / G12f))));
	double GZXy = GXYy; //transverse isotropy
	double GYZy = Gm / (1 - (sqrt(vf) * (1 - (Gm / G23f))));
	double PRXYy = u12f * vf + PRm * (1 - vf);
	double PRZXy = PRXYy; //transverse isotropy
	double PRYZy = EYy / (2 * GYZy) - 1;

	double F1t = F1f * (vf + (Em / E1f) * (1 - vf));	//longitudinal tensile strength
	double F2t = Ftm * (1 + (vf - sqrt(vf)) * (1 - (Em / E2f)));	//transverse tensile strength
	double F3t = F2t;	//transverse isotropy
	double F1c = vf * F2f;
	double F2c = Fcm * (1 + (vf - sqrt(vf)) * (1 - (Em / E2f)));	//transverse compressive strength
	double F3c = F2c; //transverse isotropy
	double F4 = Fsm * ((1 - (sqrt(vf) * (1 - (Gm / G23f)))) / (1 - (vf * (1 - (Gm / G23f)))));
	double F6 = Fsm * (1 + (vf - sqrt(vf)) * (1 - (Gm / G12f)));
	double F5 = F6;	//uniform thickness of unit-cell; plane sections remain plane
	
	cout<<"\nE1\t"<<EXy<<"\n"<<"E2\t"<<EYy<<"\n"<<"v12\t"<<PRXYy<<"\n"<<"v23\t"<<PRYZy<<"\n"<<"G12\t"<<GXYy<<"\n"<<"G23\t"<<GYZy<<"\n";
}

int main()
{
	materialPropertiesChamis();
	return 0;
}