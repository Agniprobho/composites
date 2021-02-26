void laminate_ABDH(double E1, double E2, double G12, double G23, double v12, int N, double THK, double *THETA, double A[3][3], double B[3][3], double D[3][3], double H[2][2])
{
/*
	Author: Agniprobho Mazumder 
	Purpose: Calculates laminate stiffness matrices A,B,D,H, macrolevel analysis
	Reference: Introduction to Composite Materials Design, E.J. Barbero
	Parameters: 
		Input, E1, E2, G12, G23, v12 of transversely isotropic composite laminate
		Input, THK, thickness of each lamina
		Input, N, number of laminas
		Input, THETA[N], Laminate stacking sequence
		Output, A,B,D,H matrices
	NOTE: Code has further room for optimization. This code assumes uniform thickness of all laminas
*/
	double q11, q22, q12, q66, q44, q55;
	double a1, a2, a3, b1, b2, c1, c2, c3, d1, d2, e1, e2, f1, f2, g1, g2, h1, h2, i1, i2;
	double m, n, th;
	double zk, coeff_A, coeff_B, coeff_D, coeff_H;
	double Q11, Q22, Q12, Q16, Q26, Q66, Q44, Q55, Q45;
	double delta;
	
	delta = 1 - (v12*v12*(E2/E1));
	q11 = E1/delta;
	q22 = E2/delta;
	q12 = v12*(E2/delta);
	q66 = G12;
	q44 = G23;
	q55 = G12;
	
	a1 = q11;
	a2 = 2*(q12 + 2*q66);
	a3 = q22;
	
	b1 = q11 + q22 -4*q66;
	b2 = q12;
	
	c1 = a1;
	c2 = a2;
	c3 = a3;
	
	d1 = q11 - q12 - 2*q66;
	d2 = q12 - q22 + 2*q66;
	
	e1 = d1;
	e2 = d2;
	
	f1 = q11 + q22 - 2*q12 - 2*q66;
	f2 = q66;
	
	g1 = q44;
	g2 = q55;
	
	h1 = g1;
	h2 = g2;
	
	i1 = q55 - q44;
	
	for (int k=0; k<N; k++)	// each lamina layer
	{
		th = (3.141593*THETA[k])/180.;
		m = cos(th);
		n = sin(th);
		
		zk = (-N*(THK/2.) + k*THK + THK/2.);
		
		coeff_A = THK;
		coeff_B = -THK*zk;
		coeff_D = (THK*zk*zk + (THK*THK*THK)/12);
		coeff_H = -(5/4)*(THK - (4*coeff_D)/(THK*THK));
		
		Q11 = a1*(m*m*m*m) + a2*(m*m*n*n) + a3*(n*n*n*n);
		Q12 = b1*(n*n*m*m) + b2*(m*m*m*m + n*n*n*n);
		Q22 = c1*(n*n*n*n) + c2*(m*m*n*n) + c3*(m*m*m*m);
		Q16 = d1*(n*m*m*m) + d2*(m*n*n*n);
		Q26 = e1*(n*n*n*m) + e2*(m*m*m*n);
		Q66 = f1*(m*m*n*n) + f2*(m*m*m*m + n*n*n*n);
		Q44 = g1*(m*m) + g2*(n*n);
		Q55 = h1*(n*n) + h2*(m*m);
		Q45 = i1*(m*n);
		
		//Matrix A
		A[0][0] += Q11*coeff_A;	//A11
		A[0][1] += Q12*coeff_A;	//A12
		A[0][2] += Q16*coeff_A;	//A16
		A[1][0] += Q12*coeff_A;	//A12
		A[1][1] += Q22*coeff_A;	//A22
		A[1][2] += Q26*coeff_A;	//A26
		A[2][0] += Q16*coeff_A;	//A16
		A[2][1] += Q26*coeff_A;	//A26
		A[2][2] += Q66*coeff_A;	//A66
		
		//Matrix B
		B[0][0] += Q11*coeff_B;	//B11
		B[0][1] += Q12*coeff_B;	//B12
		B[0][2] += Q16*coeff_B;	//B16
		B[1][0] += Q12*coeff_B;	//B12
		B[1][1] += Q22*coeff_B;	//B22
		B[1][2] += Q26*coeff_B;	//B26
		B[2][0] += Q16*coeff_B;	//B16
		B[2][1] += Q26*coeff_B;	//B26
		B[2][2] += Q66*coeff_B;	//B66
		
		//Matrix D
		D[0][0] += Q11*coeff_D;	//D11
		D[0][1] += Q12*coeff_D;	//D12
		D[0][2] += Q16*coeff_D;	//D16
		D[1][0] += Q12*coeff_D;	//D12
		D[1][1] += Q22*coeff_D;	//D22
		D[1][2] += Q26*coeff_D;	//D26
		D[2][0] += Q16*coeff_D;	//D16
		D[2][1] += Q26*coeff_D;	//D26
		D[2][2] += Q66*coeff_D;	//D66
		
		//Matrix H
		H[0][0] += Q44*coeff_H;	//H44
		H[0][1] += Q45*coeff_H;	//H45
		H[1][0] += Q45*coeff_H;	//H45
		H[1][1] += Q55*coeff_H;	//H55
	}
}

