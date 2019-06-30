/*
=================================================
|   circle_fitter_2    C++                      |
|   Author:            Luca Ciucci              |
|   Version:           0                        |
|   Date:              x / x / 2019             |
|   e-mail:            develop@lucaciucci99.com |
|   licese:            license.md               |
=================================================
*/

template <class T>
T sqr(T x)
{
	return x * x;
}

#define DEBUG
#include "SFML/Graphics.hpp"
#include "LC/LC.h"
//#include "D:\LUCA\codici\vs\algoritmi\general004\general004\LC\LC.h"
#include "haloFitter.h"
#include "parabola_fit.h"
#include "H5Cpp.h"

//using namespace lc::math;

double asd(double x, double y)
{
	return sqr(x - lc::math::m_PI<> * 0.0) + sqr(y - 1.0);
}

std::vector<double> p_tmp = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 10 };

double dummyChi2(double x)
{
	double sum_v = 0.0;
	for (int i = 0; i < p_tmp.size(); i++)
		sum_v += sqr(p_tmp[i] - x);
	return sum_v;
}

double l_i[3] = { 600, 525, 450 };
double dl_i[3] = { 50, 50, 50 };
double n_i[3] = { 1.3297, 1.33148, 1.33633 };
double dn_i[3] = { 1,1,1 };// { 0.00012, 0.00012, 0.00031 };

double cauchy_chi2(double A, double n0)
{
	double sum_v = 0.0;
	for (int i = 0; i < 3; i++)
	{
		sum_v += sqr(n_i[i] - n0 - A / sqr(l_i[i])) / sqr(dn_i[i]);
	}
	return sum_v;
}// http://optics.sgu.ru/_media/optics/staff/bashkatov/bashkatov_spie_03_5068_393.pdf

int main(int argc, char** argv)
{
	/*while (1)
	{
		HaloFitter<> haloFitter("lunar-halo_despeckle.jpg");
		//HaloFitter<> haloFitter;
		haloFitter.takeUserPoints();
		haloFitter.setDRadius(25);
		haloFitter.setMinLevel(00.0);
		haloFitter.refine(0, 100);
		haloFitter.showResult();
	}*/

	for (int i = 1; i < 256; i++)
	{
		std::cout << "myColor {" << inferno[i].r << ", " << inferno[1].g << ", " << inferno[i].b << ", " << (double)i / 255 << "}," << std::endl;
	}
	VectorN<5> v = 0.0;
	std::cout << lc::math::OptimizeFunction(asd, VectorN<2>(0.0), 1.0, 100);
	std::cout << lc::math::getHessianMatrix(asd, VectorN<2>(0.0));

	// https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470824566.app1
	VectorN</*dimensione*/1> params; params[0] = 0.0;
	params = lc::math::OptimizeFunction(dummyChi2, params, 0.5, 1000);
	double minChi2 = dummyChi2(params[0]);
	double maxChi2 = dummyChi2(params[0]) * ((p_tmp.size() - 1) + sqrt(2.0 * (p_tmp.size() - 1))) / (p_tmp.size() - 1);
	double sigma2 = minChi2 / ((double)p_tmp.size() - 1);
	lc::math::Mat<1, 1> HessianMatrix = lc::math::getHessianMatrix(dummyChi2, params);
	HessianMatrix.M[0][0] /= sigma2;
	lc::math::Mat<1, 1> pCov; pCov.M[0][0] = 1.0/HessianMatrix.M[0][0];
	std::cout << pCov << std::endl;
	std::cout << params[0] << "+-" << sqrt(pCov.M[0][0] * 2.0) << std::endl;

	std::cout << "test x2" << std::endl;
	params; params[0] = 0.0;
	params = lc::math::OptimizeFunction(dummyChi2, params, 0.5, 1000);
	double factor = dummyChi2(params[0]) / (p_tmp.size()-1);
	minChi2 = dummyChi2(params[0]) / factor;
	//maxChi2 = dummyChi2(params[0]) / factor * ((p_tmp.size() - 1) + sqrt(2.0 * (p_tmp.size() - 1))) / (p_tmp.size() - 1);
	//sigma2 = maxChi2 - minChi2;
	sigma2 = 2.0 * (p_tmp.size() - 1);
	HessianMatrix = lc::math::getHessianMatrix(dummyChi2, params);
	HessianMatrix.M[0][0] /= factor;
	HessianMatrix.M[0][0] /= sqrt(sigma2);
	pCov.M[0][0] = 1.0 / HessianMatrix.M[0][0];
	std::cout << "x2 = " << minChi2 << std::endl;
	std::cout << "sigma2 = " << sigma2 << std::endl;
	std::cout << "maxX2 = " << minChi2 + sqrt(sigma2) << std::endl;
	std::cout << "x2 @ max params(" << params[0] + sqrt(pCov.M[0][0] * 2.0) << ") = " << dummyChi2(params[0] + sqrt(pCov.M[0][0] * 2.0))/ factor << std::endl;
	std::cout << params[0] << "+-" << sqrt(pCov.M[0][0] * 2.0) << std::endl;


	// test
	int N = 1000;// punti
	int nTests = 1000;// prove
	double centro_dist = 5.0;
	double var_dist = 2.0;
	double totalChi2 = 0.0;
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(/*media*/centro_dist,/*sigma*/ 2.0);
	for (int t = 0; t < nTests; t++)
	{
		p_tmp.clear();
		for (int i = 0; i < N; i++)
		{
			p_tmp.push_back(distribution(generator));
			//p_tmp.push_back
		}
		// avvia fit
		params; params[0] = 0.0;
		params = lc::math::OptimizeFunction(dummyChi2, params, 0.5, 1000);
		//factor = dummyChi2(params[0]) / (p_tmp.size() - 1);
		minChi2 = dummyChi2(params[0]);
		factor = minChi2 / ((double)p_tmp.size() - 1);
		//minChi2 = dummyChi2(params[0]) / factor;
		sigma2 = 2.0 * (p_tmp.size() - 1);
		HessianMatrix = lc::math::getHessianMatrix(dummyChi2, params);
		//HessianMatrix.M[0][0] /= factor;
		//HessianMatrix.M[0][0] /= sqrt(sigma2);
		HessianMatrix.M[0][0] /= factor;
		pCov.M[0][0] = 1.0 / HessianMatrix.M[0][0];
		std::cout << params[0] << "+-" << sqrt(pCov.M[0][0] * 2.0) << std::endl;
		
		totalChi2 += sqr(params[0] - centro_dist) / (pCov.M[0][0] * 2.0);
	}
	std::cout << "chi2 totale = " << totalChi2 << std::endl;
	// !test

	// prova fit cauchy
	VectorN<2> cauchyParams = 0.0;
	cauchyParams = lc::math::OptimizeFunction(cauchy_chi2, cauchyParams, 0.1, 10000);
	std::cout << "parametri di cauchy: " << cauchyParams << std::endl;
	std::cout << "x2: " << cauchy_chi2(cauchyParams[0], cauchyParams[1]) << std::endl;

	// !! PROGRAMMA
	std::string fileName;
	if (argc == 2)
	{
		fileName = argv[1];
	}
	else
	{
		std::cout << "nome file ?" << std::endl;
		std::cin >> fileName;
		if (fileName == "1") fileName = "lunar-halo.jpg";
		if (fileName == "2") fileName = "lunar-halo_despeckle.jpg";
		if (fileName == "3") fileName = "double-rainbow.jpg";
		if (fileName == "4") fileName = "double-rainbow_meno_15px_gaussianBlur.jpg";
		if (fileName == "5") fileName = "lunar-halo_contrast.jpg";
		if (fileName == "6") fileName = "lunar-halo_despeckle_contrast.jpg";
	}
	std::cout << "Color? (0..2 <-> myColor)" << std::endl;
	int color; std::cin >> color;
	std::cout << "Radius tolerance in px? (0 -> default = 25)" << std::endl;
	int dr; std::cin >> dr; if (dr <= 0) dr = 25;
	std::cout << "Round angle subdivision for refinement? (0 -> default = 100)" << std::endl;
	int nSub; std::cin >> nSub; if (nSub <= 0) nSub = 100;

	HaloFitter<> haloFitter(fileName);
	haloFitter.getUserPoints();
	haloFitter.setDRadius(dr);
	//haloFitter.setMinLevel(00.0);
	std::cout << "dati dell'utente risultano nel fit: " << haloFitter.UsrFitCircle << std::endl;
	haloFitter.refine(color, nSub);
	haloFitter.showResult();
	haloFitter.plotRes();
	return 0;
}
// TODO refine fit come fatto in covMatrix


/*
+/- moltiplica il colore, F1 reset colore
hold R/G/B per vedere R/G/B
F2 per disabilitare visualizzazione del cerchio

*/