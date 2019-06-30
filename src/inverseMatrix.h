#ifndef LC_INVERSEMATRIX_H
#define LC_INVERSEMATRIX_H



#include <conio.h>
#include <iostream>
#include <math.h>

//----------------------------------------------------------
//read matrix : cin>> in c++ , scanf() in c
int scanf(double a[100][100]) {
	int i, j, n;
	std::cout << "\n Enter Length Of Matrix N*N : ";
	std::cin >> n;
	std::cout << "\n--------------------------\n";
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++) {
			std::cout << " Matrix[" << i + 1 << "][" << j + 1 << "] : ";
			std::cin >> a[i][j];
		}
	std::cout << "\n----------------------------------------------------\n";
	return n;
}

//-----------------------------------------------------
// show matrix : cout<< in c++ , printf() in c
void printf(double a[100][100], int n, int show) {
	int i, j;
	if (show == 1)
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++)
				std::cout << " " << a[i][j] << " \t";
			std::cout << "\n";
		}
	else if (show == 2) {
		std::cout << "\n\n The Inverse Of Matrix Is : \n\n";
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++)
				std::cout << " " << a[i][j] << " \t";
			std::cout << "\n";
		}
	}
}

//---------------------------------------------------
//	calculate minor of matrix OR build new matrix : k-had = minor
void minor(double b[100][100], double a[100][100], int i, int n) {
	int j, l, h = 0, k = 0;
	for (l = 1; l < n; l++)
		for (j = 0; j < n; j++) {
			if (j == i)
				continue;
			b[h][k] = a[l][j];
			k++;
			if (k == (n - 1)) {
				h++;
				k = 0;
			}
		}
}// end function

//---------------------------------------------------
//	calculate determinte of matrix
double det(double a[100][100], int n) {
	int i;
	double b[100][100], sum = 0;
	if (n == 1)
		return a[0][0];
	else if (n == 2)
		return (a[0][0] * a[1][1] - a[0][1] * a[1][0]);
	else
		for (i = 0; i < n; i++) {
			minor(b, a, i, n);	// read function
			sum = (double)(sum + a[0][i] * pow(-1, i) * det(b, (n - 1)));	// read function	// sum = determinte matrix
		}
	return sum;
}// end function

//---------------------------------------------------
//	calculate transpose of matrix
void transpose(double c[100][100], double d[100][100], int n, double det) {
	int i, j;
	double b[100][100];
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			b[i][j] = c[j][i];
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			d[i][j] = b[i][j] / det;	// array d[][] = inverse matrix
}// end function

//---------------------------------------------------
//	calculate cofactor of matrix
void cofactor(double a[100][100], double d[100][100], int n, double determinte) {
	double b[100][100], c[100][100];
	int l, h, m, k, i, j;
	for (h = 0; h < n; h++)
		for (l = 0; l < n; l++) {
			m = 0;
			k = 0;
			for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)
					if (i != h && j != l) {
						b[m][k] = a[i][j];
						if (k < (n - 2))
							k++;
						else {
							k = 0;
							m++;
						}
					}
			c[h][l] = (double)pow(-1, (h + l)) * det(b, (n - 1));	// c = cofactor Matrix
		}
	transpose(c, d, n, determinte);	// read function
}// end function

//---------------------------------------------------
//	calculate inverse of matrix
void inverse(double a[100][100], double d[100][100], int n, double det) {
	if (det == 0)
		std::cout << "\nInverse of Entered Matrix is not possible\n";
	else if (n == 1)
		d[0][0] = 1;
	else
		cofactor(a, d, n, det);	// read function
}// end function

//---------------------------------------------------
//main fuction exe
int main_invmat(void) {
	int i, j, n;
	double a[100][100], d[100][100], deter;
	std::cout << "\n C++ Program To Find Inverse Of Matrix\n\n";
	n = scanf(a);	// read function
	int print_matrix = 1;
	printf(a, n, print_matrix);	// read function
	deter = (double)det(a, n);	// read function
	std::cout << "----------------------------------------------------\n";
	std::cout << "\n\n Determinant of the Matrix : " << deter;
	std::cout << "\n\n-----------------------\n";
	inverse(a, d, n, deter);	// read function
	int print_inverse = 2;
	printf(d, n, print_inverse);	// read function
	std::cout << "\n\n==============================* THE END *==============================\n";
	std::cout << "\n		**** Thanks For Using The Program!!! ****\n";
	//getch();
	return 0;
}// end main












#endif // !LC_INVERSEMATRIX_H
