// 三原色＋白色から変換マトリックスを求める方法
// boostを使用

#include <iostream>
#include "boost/numeric/ublas/lu.hpp"

using namespace std;
using namespace boost::numeric::ublas;

//逆行列を求めます
void InvertMatrix(matrix<double>&mat, matrix<double>&matinv)
{
	permutation_matrix<> pm(mat.size1());

	matrix<double> copiyMatrix = matrix<double>(mat);
	matrix<double> invMatrix = identity_matrix<double>(mat.size1());

	lu_factorize(copiyMatrix, pm);
	lu_substitute(copiyMatrix, pm, invMatrix);

	matinv.resize(invMatrix.size1(), invMatrix.size2(), false);
	matinv.assign(invMatrix);
}

//三原色と白色から変換マトリックスを求めます
int CalcXYZMatrix(double red_x, double red_y, double green_x, double green_y, double blue_x, double blue_y, double white_x, double white_y, matrix<double> &mat3x3)
{
	//マトリックスを求めます
	//rRx/Ry + gGx/Gy + bBx/By = Wx/Wy
	//r      + g      + b      = 1
	//rRz/Ry + gGz/Gy + bBz/By = Wz/Wy
	//の連立方程式からr,g,bを求めます
	matrix<double> in_mat(3,4);

	in_mat(0,0) = red_x / red_y;
	in_mat(0,1) = green_x / green_y;
	in_mat(0,2) = blue_x / blue_y;
	in_mat(0,3) = white_x / white_y;
	in_mat(1,0) = 1.;
	in_mat(1,1) = 1.;
	in_mat(1,2) = 1.;
	in_mat(1,3) = 1.;
	in_mat(2,0) = (1 - red_x - red_y) / red_y;
	in_mat(2,1) = (1 - green_x - green_y) / green_y;
	in_mat(2,2) = (1 - blue_x - blue_y) / blue_y;
	in_mat(2,3) = (1 - white_x - white_y) / white_y;

	//Gauss-Jordan
	for(int k = 0; k < in_mat.size1(); k++){
		double p=in_mat(k,k);
		int i;
		for (i = k; i <  in_mat.size2(); i++){
			in_mat(k,i) /= p;
		}
		
		for (i = 0; i < in_mat.size1(); i++){
			if(i != k){
				double d = in_mat(i,k);
				for(int j = k; j < in_mat.size2(); j++)
					in_mat(i,j) -= d * in_mat(k,j);
			}
		}
	}

	matrix<double> out_mat(3,3);
	out_mat(0,0) = in_mat(0,3) * red_x / red_y;
	out_mat(0,1) = in_mat(1,3) * green_x / green_y;
	out_mat(0,2) = in_mat(2,3) * blue_x / blue_y;
	out_mat(1,0) = in_mat(0,3);
	out_mat(1,1) = in_mat(1,3);
	out_mat(1,2) = in_mat(2,3);
	out_mat(2,0) = in_mat(0,3) * (1 - red_x - red_y) / red_y;
	out_mat(2,1) = in_mat(1,3) * (1 - green_x - green_y) / green_y;
	out_mat(2,2) = in_mat(2,3) * (1 - blue_x - blue_y) / blue_y;

	mat3x3 = out_mat;
}

void dump(matrix<double>&in_mat, matrix<double>&out_mat)
{
	cout << "X = " << in_mat(0,0) << " * Red + " << in_mat(0,1) << " * Green + " << in_mat(0,2) << " * Blue" << endl;
	cout << "Y = " << in_mat(1,0) << " * Red + " << in_mat(1,1) << " * Green + " << in_mat(1,2) << " * Blue" << endl;
	cout << "Z = " << in_mat(2,0) << " * Red + " << in_mat(2,1) << " * Green + " << in_mat(2,2) << " * Blue" << endl;
	cout << endl;
	cout << "Red = " << out_mat(0,0) << " * X + " << out_mat(0,1) << " * Y + " << out_mat(0,2) << " * Z" << endl;
	cout << "Green = " << out_mat(1,0) << " * X + " << out_mat(1,1) << " * Y + " << out_mat(1,2) << " * Z" << endl;
	cout << "Blue = " << out_mat(2,0) << " * X + " << out_mat(2,1) << " * Y + " << out_mat(2,2) << " * Z" << endl;
	cout << endl;
}

int main(int argc, char **argv)
{
	double red_x;
	double red_y;
	double green_x;
	double green_y;
	double blue_x;
	double blue_y;
	double white_x;
	double white_y;
	matrix<double> matRGBtoXYZ;
	matrix<double> matXYZtoRGB;

	matRGBtoXYZ = matrix<double>(3,3);
	matXYZtoRGB = matrix<double>(3,3);

	/*
	//sRGB(Illuminant D65)の場合、以下の値になっているはずです
	matRGBtoXYZ(0,0) = .412424;
	matRGBtoXYZ(0,1) = .357579;
	matRGBtoXYZ(0,2) = .180464;
	matRGBtoXYZ(1,0) = .212656;
	matRGBtoXYZ(1,1) = .715158;
	matRGBtoXYZ(1,2) = .0721856;
	matRGBtoXYZ(2,0) = .0193324;
	matRGBtoXYZ(2,1) = .119193;
	matRGBtoXYZ(2,2) = .950444;
	*/

	//sRGB (D65)
	red_x = 0.64;
	red_y = 0.33;
	green_x = 0.30;
	green_y = 0.60;
	blue_x = 0.15;
	blue_y = 0.06;
	white_x = 0.3128;
	white_y = 0.3292;

	CalcXYZMatrix(red_x, red_y, green_x, green_y, blue_x, blue_y, white_x, white_y, matRGBtoXYZ);
	InvertMatrix(matRGBtoXYZ, matXYZtoRGB);

	cout << "sRGB(Illuminant D65)" << endl;
	dump(matRGBtoXYZ, matXYZtoRGB);

	//AdobeRGB
	red_x = 0.64;
	red_y = 0.33;
	green_x = 0.21;
	green_y = 0.71;
	blue_x = 0.15;
	blue_y = 0.06;
	white_x = 0.3127;
	white_y = 0.3290;

	CalcXYZMatrix(red_x, red_y, green_x, green_y, blue_x, blue_y, white_x, white_y, matRGBtoXYZ);
	InvertMatrix(matRGBtoXYZ, matXYZtoRGB);

	cout << "AdobeRGB" << endl;
	dump(matRGBtoXYZ, matXYZtoRGB);

	return 0;
}
