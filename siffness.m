E1=input('E1: ');
E2=input('E2: ');
E3=E2;
G12=input('G12: ');
G23=input('G23: (enter 0 if unknown)');

if G23==0
	G23=0.67*G12;
end

G13=G12;
p12=input('p12: ');
p23=input('p23: (enter 0 if unknown)');

if p23==0
	p23=(E2/(2*G23))-1;
end

p13=p12;

S1=[1/E1 -p12/E1 -p13/E1 0 0 0;-p12/E1 1/E2 -p23/E2 0 0 0;-p13/E1 -p23/E2 1/E3 0 0 0;0 0 0 1/G23 0 0;0 0 0 0 1/G13 0;0 0 0 0 0 1/G12];
x1=input('orientation');
x=(pi/180)*x1;
sum=[0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0];

for i=1:length(x)
	l1(i)=cos(x(i));
	m1(i)=sin(x(i));
	l2(i)=-sin(x(i));
	m2(i)=cos(x(i));
	n3=1;
	T=[(l1(i))^2 (m1(i))^2 0 0 0 (2*l1(i)*m1(i));(l2(i))^2 (m2(i))^2 0 0 0 (2*l2(i)*m2(i));0 0 1 0 0 0;0 0 0 m2(i) l2(i) 0;0 0 0 m1(i) l1(i) 0;(l1(i)*l2(i)) (m1(i)*m2(i)) 0 0 0 ((l1(i)*m2(i))+(m1(i)*l2(i)))];
	S2=T'*S1*T;
	C1=inv(S2);
	sum=sum+(1/length(x))*C1;
end

S=inv(sum);
Ex=double(1/S(1,1));
Ey=double(1/S(2,2));
Ez=double(1/S(3,3));
Gyz=double(1/S(4,4));
Gxz=double(1/S(5,5));
Gxy=double(1/S(6,6));
pxy=-Ex*S(1,2);
pyz=-Ey*S(2,3);
pzx=-Ex*S(1,3);
fprintf('Ex=%0.0f\nEy=%0.0f\nEz=%0.0f\nGyz=%0.0f\nGxz=%0.0f\nGxy=%0.0f\npxy=%0.3f\npyz=%0.3f\npzx=%0.3f\n',Ex,Ey,Ez,Gyz,Gxz,Gxy,pxy,pyz,pzx);