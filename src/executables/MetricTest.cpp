#include "Metric.h"
#include <iostream>

using namespace std;

void print_test(const double result, const double expected){
	cout << result << " ";
	bool pass = fabs(result-expected) / (fabs(expected)>0 ? fabs(expected) : 1.0) < TINY;
	if(pass) cout << endl;
	else cout << "\tFAIL: expected " << expected << endl;
}

int main(){

	cout << "|==================|" << endl;
	cout << "| 3 Minkowski Test |" << endl;
	cout << "|==================|" << endl;
	ThreeMetric g3;
	g3.data[ixx] = 1;
	g3.data[iyy] = 1;
	g3.data[izz] = 1;
	g3.data[ixy] = 0;
	g3.data[ixz] = 0;
	g3.data[iyz] = 0;
	Tuple<double,3> kup;
	kup[0] = 1;
	kup[1] = 2;
	kup[2] = 3;
	cout << "kup={" << kup[0] << ","<<kup[1]<<","<<kup[2]<<"}" << endl;

	double det = g3.det();
	cout << " * Determinant=";
	print_test(det,1.0);

	Tuple<double,3> klow = g3.lower(kup);
	cout << " * klow[0]=";
	print_test(klow[0],1.);
	cout << " * klow[1]=";
	print_test(klow[1],2.);
	cout << " * klow[2]=";
	print_test(klow[2],3.);

	ThreeMetric g3inv = g3.inverse();
	cout << " * g3inv[xx]=";
	print_test(g3inv.data[ixx],1.0);
	cout << " * g3inv[xy]=";
	print_test(g3inv.data[ixy],0.0);
	cout << " * g3inv[xz]=";
	print_test(g3inv.data[ixz],0.0);
	cout << " * g3inv[yy]=";
	print_test(g3inv.data[iyy],1.0);
	cout << " * g3inv[yz]=";
	print_test(g3inv.data[iyz],0.0);
	cout << " * g3inv[zz]=";
	print_test(g3inv.data[izz],1.0);


	cout << "|======================|" << endl;
	cout << "| 3 Schwarzschild Test |" << endl;
	cout << "|======================|" << endl;
	g3.data[ixx] = 0.5;
	g3.data[iyy] = 1;
	g3.data[izz] = 1;
	g3.data[ixy] = 0;
	g3.data[ixz] = 0;
	g3.data[iyz] = 0;
	cout << "kup={" << kup[0] << ","<<kup[1]<<","<<kup[2]<<"}" << endl;

	det = g3.det();
	cout << " * Determinant=";
	print_test(det,0.5);

	klow = g3.lower(kup);
	cout << " * klow[0]=";
	print_test(klow[0],0.5);
	cout << " * klow[1]=";
	print_test(klow[1],2);
	cout << " * klow[2]=";
	print_test(klow[2],3);

	g3inv = g3.inverse();
	cout << " * g3inv[xx]=";
	print_test(g3inv.data[ixx],2.0);
	cout << " * g3inv[xy]=";
	print_test(g3inv.data[ixy],0.0);
	cout << " * g3inv[xz]=";
	print_test(g3inv.data[ixz],0.0);
	cout << " * g3inv[yy]=";
	print_test(g3inv.data[iyy],1.0);
	cout << " * g3inv[yz]=";
	print_test(g3inv.data[iyz],0.0);
	cout << " * g3inv[zz]=";
	print_test(g3inv.data[izz],1.0);

// Cholesky decomposition requires a positive definite matrix
//	cout << "|==============|" << endl;
//	cout << "| 3 Weird Test |" << endl;
//	cout << "|==============|" << endl;
//	g3.data[ixx] = 0;
//	g3.data[iyy] = 0;
//	g3.data[izz] = 0;
//	g3.data[ixy] = 1;
//	g3.data[ixz] = 1;
//	g3.data[iyz] = 1;
//	cout << "kup={" << kup[0] << ","<<kup[1]<<","<<kup[2]<<"}" << endl;
//
//	det = g3.det();
//	cout << " * Determinant=";
//	print_test(det,2.0);
//
//	klow = g3.lower(kup);
//	cout << " * klow[0]=";
//	print_test(klow[0],5);
//	cout << " * klow[1]=";
//	print_test(klow[1],4);
//	cout << " * klow[2]=";
//	print_test(klow[2],3);
//
//	g3inv = g3.inverse();
//	cout << " * g3inv[xx]=";
//	print_test(g3inv.data[ixx],-0.5);
//	cout << " * g3inv[xy]=";
//	print_test(g3inv.data[ixy],0.5);
//	cout << " * g3inv[xz]=";
//	print_test(g3inv.data[ixz],0.5);
//	cout << " * g3inv[yy]=";
//	print_test(g3inv.data[iyy],-0.5);
//	cout << " * g3inv[yz]=";
//	print_test(g3inv.data[iyz],0.5);
//	cout << " * g3inv[zz]=";
//	print_test(g3inv.data[izz],-0.5);
//
//	kup = g3inv.lower(klow);
//	cout << " * kup[0]=";
//	print_test(kup[0],1);
//	cout << " * kup[1]=";
//	print_test(kup[1],2);
//	cout << " * kup[2]=";
//	print_test(kup[2],3);

	cout << "|==================|" << endl;
	cout << "| 4 Minkowski Test |" << endl;
	cout << "|==================|" << endl;
	Metric g;
	g.gammalow.data[ixx] = 1;
	g.gammalow.data[iyy] = 1;
	g.gammalow.data[izz] = 1;
	g.gammalow.data[ixy] = 0;
	g.gammalow.data[ixz] = 0;
	g.gammalow.data[iyz] = 0;
	for(size_t i=0; i<3; i++) g.betaup[i] = 0;
	g.alpha = 1.0;
	g.update();
	Tuple<double,4> kup4;
	kup4[0] = 0;
	kup4[1] = 1;
	kup4[2] = 2;
	kup4[3] = 3;
	cout << "kup={" << kup4[0] << ","<<kup4[1]<<","<<kup4[2]<<","<<kup4[3]<<"}" << endl;

	Tuple<double,4> klow4 = g.lower<4>(kup4);
	cout << " * klow[0]=";
	print_test(klow4[0],0.);
	cout << " * klow[1]=";
	print_test(klow4[1],1.);
	cout << " * klow[2]=";
	print_test(klow4[2],2.);
	cout << " * klow[3]=";
	print_test(klow4[3],-3.);

	double dotp = g.dot<4>(kup4,kup4);
	cout << " * kup4.kup4=";
	print_test(dotp,-4.);

	Tuple<double,4> tmp;
	for(size_t i=0; i<4; i++) tmp[i] = kup4[i];
	g.normalize(tmp);
	cout << " * norm[0]=";
	print_test(tmp[0],0.);
	cout << " * norm[1]=";
	print_test(tmp[1],1./2.);
	cout << " * norm[2]=";
	print_test(tmp[2],2./2.);
	cout << " * norm[3]=";
	print_test(tmp[3],3./2.);

	for(size_t i=0; i<4; i++) tmp[i] = kup4[i];
	g.normalize_null_changeupt(tmp);
	cout << " * nullnorm[0]=";
	print_test(tmp[0],kup4[0]);
	cout << " * nullnorm[1]=";
	print_test(tmp[1],kup4[1]);
	cout << " * nullnorm[2]=";
	print_test(tmp[2],kup4[2]);
	cout << " * nullnorm[3]=";
	print_test(tmp[3],sqrt(5.));

	Tuple<double,4> u;
	u[0] = 1;
	u[1] = 0;
	u[2] = 0;
	u[3] = 2;
	Tuple<double,4> e;
	e[0] = 1;
	e[1] = 0;
	e[2] = 0;
	e[3] = 0;
	g.orthogonalize<4>(e,u);
	cout << " * e[0]=";
	print_test(e[0],4./3.);
	cout << " * e[1]=";
	print_test(e[1],0);
	cout << " * e[2]=";
	print_test(e[2],0);
	cout << " * e[3]=";
	print_test(e[3],2./3.);
	cout << " * e.u=";
	print_test(g.dot<4>(e,u),0);

	cout << "|==============|" << endl;
	cout << "| 4 Shift Test |" << endl;
	cout << "|==============|" << endl;
	g.gammalow.data[ixx] = 1;
	g.gammalow.data[iyy] = 1;
	g.gammalow.data[izz] = 1;
	g.gammalow.data[ixy] = 0;
	g.gammalow.data[ixz] = 0;
	g.gammalow.data[iyz] = 0;
	g.betaup[0] = 0.1;
	g.betaup[1] = 0.2;
	g.betaup[2] = 0.3;
	g.alpha = 1.0;
	g.update();
	kup4[0] = 0;
	kup4[1] = 1;
	kup4[2] = 2;
	kup4[3] = 3;
	cout << "kup={" << kup4[0] << ","<<kup4[1]<<","<<kup4[2]<<","<<kup4[3]<<"}" << endl;

	cout << " * betalow[0]=";
	print_test(g.betalow[0],0.1);
	cout << " * betalow[1]=";
	print_test(g.betalow[1],0.2);
	cout << " * betalow[2]=";
	print_test(g.betalow[2],0.3);
	dotp = g.dot<3>(g.betaup,g.betaup);
	cout << " * beta.beta=";
	print_test(dotp,.14);

	double gtt = g.gtt;
	cout << " * gtt=";
	print_test(gtt, -0.86);

	klow4 = g.lower<4>(kup4);
	cout << " * klow4[0]=";
	print_test(klow4[0],0.3);
	cout << " * klow4[1]=";
	print_test(klow4[1],1.6);
	cout << " * klow4[2]=";
	print_test(klow4[2],2.9);
	cout << " * klow4[3]=";
	print_test(klow4[3],-1.78);

	dotp = g.dot<4>(kup4,kup4);
	cout << " * kup4.kup4=";
	print_test(dotp,2.06);

	for(size_t i=0; i<4; i++) tmp[i] = kup4[i];
	g.normalize(tmp);
	cout << " * norm[0]=";
	print_test(tmp[0],0./sqrt(dotp));
	cout << " * norm[1]=";
	print_test(tmp[1],1./sqrt(dotp));
	cout << " * norm[2]=";
	print_test(tmp[2],2./sqrt(dotp));
	cout << " * norm[3]=";
	print_test(tmp[3],3./sqrt(dotp));

	for(size_t i=0; i<4; i++) tmp[i] = kup4[i];
	g.normalize_null_changeupt(tmp);
	dotp = g.dot<4>(tmp,tmp);
	cout << " * null.null=";
	print_test(dotp,0);
	cout << "{ ";
	for(size_t i=0; i<4; i++) cout << tmp[i] << " ";
	cout << "}" << endl;

	for(size_t i=0; i<4; i++) tmp[i] = kup4[i];
	g.normalize_null_preserveupt(tmp);
	dotp = g.dot<4>(tmp,tmp);
	cout << " * null.null=";
	print_test(dotp,0);
	cout << "{ ";
	for(size_t i=0; i<4; i++) cout << tmp[i] << " ";
	cout << "}" << endl;

	return 0;
}
