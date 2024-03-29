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
#ifndef LC_HALOFITTER_H
#define LC_HALOFITTER_H


#include <vector>
#include <cmath>
#include <math.h>
#include "LC/LC.h"
#include <string>

using namespace lc::math;

template <class T, class T1, class T2>
T clamp(T v, T1 min, T2 max)
{
	return (v > max) ? max : ((v < min) ? min : v);
}

template <class T, class T1, class T2>
bool isIninterval(T v, T1 min, T2 max)
{
	return v >= min && v <= max;
}

class Retta {
public:
	// ax + by + c = 0
	double a = 0, b = 0, c = 0;

	Retta()
	{
		a = 0; b = 0; c = 0;
	}
	Retta(double a_, double b_, double c_)
	{
		a = a_; b = b_, c = c_;
	}

	double evalAt(double x, double y)
	{
		return a * x + b * y + c;
	}
};

inline int sign(double x)
{
	return ((x == 0.0) ? 0 : ((x > 0.0) ? 1 : -1));
}


template <class LetVoid = void>
class HaloFitter {
public:
	// ================================
	//            CONSTRUCTOR
	// ================================
	//HaloFitter(std::string imgName = "double-rainbow.jpg");
	HaloFitter(std::string imgName = "double-rainbow_meno_15px_gaussianBlur.jpg");
	~HaloFitter() {};

	// ================================
	//            FUNCTIONS
	// ================================
	void getUserPoints(void);
	void plot(bool plotUser = false);
	void plotUserView(void);
	void plotComputerView(void);
	void updateUserPts(void);
	void refine(int color = 0, int nStep = 200);
	double& setDRadius(double new_dr = 10.0);
	double& setMinLevel(double l = 50)
	{
		return minLevel = l;
	}
	void showResult(void);
	void plotRes(void);

	// tmp
	Vector3<> UsrFitCircle;
private:
	// immagine
	sf::Texture img_texture;
	sf::Texture chi2_texture;
	sf::Uint8* img_pixels;
	sf::Uint8* tmp_pixels;
	sf::Uint8* mask_pixels;
	sf::Uint8* chi2Mask_pixels; int chi2MaskSize; int chi2MaskDivider;
	sf::Sprite img_sprite;
	sf::Sprite chi2_sprite;


	// window
	int wsX, wsY;
	sf::RenderWindow window;
	sf::Event event;
	double currZoom;
	double zoomFactor;
	double colorZoom;
	bool showCircle;
	sf::View initialView;

	// user points
	std::vector<Vector3<>> userPts;
	//Vector3<> UsrFitCircle;

	// refinement
	double dr;
	std::vector<Vector3<>> selectedPxs;
	std::vector<Vector3<>> computerPts;
	Vector3<> computerFitCircle;
	sf::Color maskColor;
	double minLevel;

	// tmp chi2
	double hM[100][100];
	double covM[100][100];
#define MYEPSYLON 0.001
	void hessianMatrix(double xc, double yx, double r, double h = MYEPSYLON);
#undef MYEPSYLON
	void computeCovMatrix(void);

	// FUNCTIONS
	Vector3<> fastCircleFit(std::vector<Vector3<>>& points);
	void resetZoom(void);
	void computeChi2Map(void);
	sf::Color toHeatmap(double value, double min, double Max);

};


double chi2_on_r(double R);

struct myColor
{
	myColor(double _r, double _g, double _b, double _a)
	{
		r = _r;
		g = _g;
		b = _b;
		a = _a;
	}
	double r, g, b, a;
};
std::vector<myColor> inferno{
	myColor {0.002267, 0.00127, 0.01857, 0.00392157},
myColor {0.003299, 0.00127, 0.024239, 0.00784314},
myColor {0.004547, 0.00127, 0.030909, 0.0117647},
myColor {0.006006, 0.00127, 0.038558, 0.0156863},
myColor {0.007676, 0.00127, 0.046836, 0.0196078},
myColor {0.009561, 0.00127, 0.055143, 0.0235294},
myColor {0.011663, 0.00127, 0.06346, 0.027451},
myColor {0.013995, 0.00127, 0.071862, 0.0313725},
myColor {0.016561, 0.00127, 0.080282, 0.0352941},
myColor {0.019373, 0.00127, 0.088767, 0.0392157},
myColor {0.022447, 0.00127, 0.097327, 0.0431373},
myColor {0.025793, 0.00127, 0.10593, 0.0470588},
myColor {0.029432, 0.00127, 0.114621, 0.0509804},
myColor {0.033385, 0.00127, 0.123397, 0.054902},
myColor {0.037668, 0.00127, 0.132232, 0.0588235},
myColor {0.042253, 0.00127, 0.141141, 0.0627451},
myColor {0.046915, 0.00127, 0.150164, 0.0666667},
myColor {0.051644, 0.00127, 0.159254, 0.0705882},
myColor {0.056449, 0.00127, 0.168414, 0.0745098},
myColor {0.06134, 0.00127, 0.177642, 0.0784314},
myColor {0.066331, 0.00127, 0.186962, 0.0823529},
myColor {0.071429, 0.00127, 0.196354, 0.0862745},
myColor {0.076637, 0.00127, 0.205799, 0.0901961},
myColor {0.081962, 0.00127, 0.215289, 0.0941176},
myColor {0.087411, 0.00127, 0.224813, 0.0980392},
myColor {0.09299, 0.00127, 0.234358, 0.101961},
myColor {0.098702, 0.00127, 0.243904, 0.105882},
myColor {0.104551, 0.00127, 0.25343, 0.109804},
myColor {0.110536, 0.00127, 0.262912, 0.113725},
myColor {0.116656, 0.00127, 0.272321, 0.117647},
myColor {0.122908, 0.00127, 0.281624, 0.121569},
myColor {0.129285, 0.00127, 0.290788, 0.12549},
myColor {0.135778, 0.00127, 0.299776, 0.129412},
myColor {0.142378, 0.00127, 0.308553, 0.133333},
myColor {0.149073, 0.00127, 0.317085, 0.137255},
myColor {0.15585, 0.00127, 0.325338, 0.141176},
myColor {0.162689, 0.00127, 0.333277, 0.145098},
myColor {0.169575, 0.00127, 0.340874, 0.14902},
myColor {0.176493, 0.00127, 0.348111, 0.152941},
myColor {0.183429, 0.00127, 0.354971, 0.156863},
myColor {0.190367, 0.00127, 0.361447, 0.160784},
myColor {0.197297, 0.00127, 0.367535, 0.164706},
myColor {0.204209, 0.00127, 0.373238, 0.168627},
myColor {0.211095, 0.00127, 0.378563, 0.172549},
myColor {0.217949, 0.00127, 0.383522, 0.176471},
myColor {0.224763, 0.00127, 0.388129, 0.180392},
myColor {0.231538, 0.00127, 0.3924, 0.184314},
myColor {0.238273, 0.00127, 0.396353, 0.188235},
myColor {0.244967, 0.00127, 0.400007, 0.192157},
myColor {0.25162, 0.00127, 0.403378, 0.196078},
myColor {0.258234, 0.00127, 0.406485, 0.2},
myColor {0.26481, 0.00127, 0.409345, 0.203922},
myColor {0.271347, 0.00127, 0.411976, 0.207843},
myColor {0.27785, 0.00127, 0.414392, 0.211765},
myColor {0.284321, 0.00127, 0.416608, 0.215686},
myColor {0.290763, 0.00127, 0.418637, 0.219608},
myColor {0.297178, 0.00127, 0.420491, 0.223529},
myColor {0.303568, 0.00127, 0.422182, 0.227451},
myColor {0.309935, 0.00127, 0.423721, 0.231373},
myColor {0.316282, 0.00127, 0.425116, 0.235294},
myColor {0.32261, 0.00127, 0.426377, 0.239216},
myColor {0.328921, 0.00127, 0.427511, 0.243137},
myColor {0.335217, 0.00127, 0.428524, 0.247059},
myColor {0.3415, 0.00127, 0.429425, 0.25098},
myColor {0.347771, 0.00127, 0.430217, 0.254902},
myColor {0.354032, 0.00127, 0.430906, 0.258824},
myColor {0.360284, 0.00127, 0.431497, 0.262745},
myColor {0.366529, 0.00127, 0.431994, 0.266667},
myColor {0.372768, 0.00127, 0.4324, 0.270588},
myColor {0.379001, 0.00127, 0.432719, 0.27451},
myColor {0.385228, 0.00127, 0.432955, 0.278431},
myColor {0.391453, 0.00127, 0.433109, 0.282353},
myColor {0.397674, 0.00127, 0.433183, 0.286275},
myColor {0.403894, 0.00127, 0.433179, 0.290196},
myColor {0.410113, 0.00127, 0.433098, 0.294118},
myColor {0.416331, 0.00127, 0.432943, 0.298039},
myColor {0.422549, 0.00127, 0.432714, 0.301961},
myColor {0.428768, 0.00127, 0.432412, 0.305882},
myColor {0.434987, 0.00127, 0.432039, 0.309804},
myColor {0.441207, 0.00127, 0.431594, 0.313725},
myColor {0.447428, 0.00127, 0.43108, 0.317647},
myColor {0.453651, 0.00127, 0.430498, 0.321569},
myColor {0.459875, 0.00127, 0.429846, 0.32549},
myColor {0.4661, 0.00127, 0.429125, 0.329412},
myColor {0.472328, 0.00127, 0.428334, 0.333333},
myColor {0.478558, 0.00127, 0.427475, 0.337255},
myColor {0.484789, 0.00127, 0.426548, 0.341176},
myColor {0.491022, 0.00127, 0.425552, 0.345098},
myColor {0.497257, 0.00127, 0.424488, 0.34902},
myColor {0.503493, 0.00127, 0.423356, 0.352941},
myColor {0.50973, 0.00127, 0.422156, 0.356863},
myColor {0.515967, 0.00127, 0.420887, 0.360784},
myColor {0.522206, 0.00127, 0.419549, 0.364706},
myColor {0.528444, 0.00127, 0.418142, 0.368627},
myColor {0.534683, 0.00127, 0.416667, 0.372549},
myColor {0.54092, 0.00127, 0.415123, 0.376471},
myColor {0.547157, 0.00127, 0.413511, 0.380392},
myColor {0.553392, 0.00127, 0.411829, 0.384314},
myColor {0.559624, 0.00127, 0.410078, 0.388235},
myColor {0.565854, 0.00127, 0.408258, 0.392157},
myColor {0.572081, 0.00127, 0.406369, 0.396078},
myColor {0.578304, 0.00127, 0.404411, 0.4},
myColor {0.584521, 0.00127, 0.402385, 0.403922},
myColor {0.590734, 0.00127, 0.40029, 0.407843},
myColor {0.59694, 0.00127, 0.398125, 0.411765},
myColor {0.603139, 0.00127, 0.395891, 0.415686},
myColor {0.60933, 0.00127, 0.393589, 0.419608},
myColor {0.615513, 0.00127, 0.391219, 0.423529},
myColor {0.621685, 0.00127, 0.388781, 0.427451},
myColor {0.627847, 0.00127, 0.386276, 0.431373},
myColor {0.633998, 0.00127, 0.383704, 0.435294},
myColor {0.640135, 0.00127, 0.381065, 0.439216},
myColor {0.64626, 0.00127, 0.378359, 0.443137},
myColor {0.652369, 0.00127, 0.375586, 0.447059},
myColor {0.658463, 0.00127, 0.372748, 0.45098},
myColor {0.66454, 0.00127, 0.369846, 0.454902},
myColor {0.670599, 0.00127, 0.366879, 0.458824},
myColor {0.676638, 0.00127, 0.363849, 0.462745},
myColor {0.682656, 0.00127, 0.360757, 0.466667},
myColor {0.688653, 0.00127, 0.357603, 0.470588},
myColor {0.694627, 0.00127, 0.354388, 0.47451},
myColor {0.700576, 0.00127, 0.351113, 0.478431},
myColor {0.7065, 0.00127, 0.347777, 0.482353},
myColor {0.712396, 0.00127, 0.344383, 0.486275},
myColor {0.718264, 0.00127, 0.340931, 0.490196},
myColor {0.724103, 0.00127, 0.337424, 0.494118},
myColor {0.729909, 0.00127, 0.333861, 0.498039},
myColor {0.735683, 0.00127, 0.330245, 0.501961},
myColor {0.741423, 0.00127, 0.326576, 0.505882},
myColor {0.747127, 0.00127, 0.322856, 0.509804},
myColor {0.752794, 0.00127, 0.319085, 0.513725},
myColor {0.758422, 0.00127, 0.315266, 0.517647},
myColor {0.76401, 0.00127, 0.311399, 0.521569},
myColor {0.769556, 0.00127, 0.307485, 0.52549},
myColor {0.775059, 0.00127, 0.303526, 0.529412},
myColor {0.780517, 0.00127, 0.299523, 0.533333},
myColor {0.785929, 0.00127, 0.295477, 0.537255},
myColor {0.791293, 0.00127, 0.29139, 0.541176},
myColor {0.796607, 0.00127, 0.287264, 0.545098},
myColor {0.801871, 0.00127, 0.283099, 0.54902},
myColor {0.807082, 0.00127, 0.278898, 0.552941},
myColor {0.812239, 0.00127, 0.274661, 0.556863},
myColor {0.817341, 0.00127, 0.27039, 0.560784},
myColor {0.822386, 0.00127, 0.266085, 0.564706},
myColor {0.827372, 0.00127, 0.26175, 0.568627},
myColor {0.832299, 0.00127, 0.257383, 0.572549},
myColor {0.837165, 0.00127, 0.252988, 0.576471},
myColor {0.841969, 0.00127, 0.248564, 0.580392},
myColor {0.846709, 0.00127, 0.244113, 0.584314},
myColor {0.851384, 0.00127, 0.239636, 0.588235},
myColor {0.855992, 0.00127, 0.235133, 0.592157},
myColor {0.860533, 0.00127, 0.230606, 0.596078},
myColor {0.865006, 0.00127, 0.226055, 0.6},
myColor {0.869409, 0.00127, 0.221482, 0.603922},
myColor {0.873741, 0.00127, 0.216886, 0.607843},
myColor {0.878001, 0.00127, 0.212268, 0.611765},
myColor {0.882188, 0.00127, 0.207628, 0.615686},
myColor {0.886302, 0.00127, 0.202968, 0.619608},
myColor {0.890341, 0.00127, 0.198286, 0.623529},
myColor {0.894305, 0.00127, 0.193584, 0.627451},
myColor {0.898192, 0.00127, 0.18886, 0.631373},
myColor {0.902003, 0.00127, 0.184116, 0.635294},
myColor {0.905735, 0.00127, 0.17935, 0.639216},
myColor {0.90939, 0.00127, 0.174563, 0.643137},
myColor {0.912966, 0.00127, 0.169755, 0.647059},
myColor {0.916462, 0.00127, 0.164924, 0.65098},
myColor {0.919879, 0.00127, 0.16007, 0.654902},
myColor {0.923215, 0.00127, 0.155193, 0.658824},
myColor {0.92647, 0.00127, 0.150292, 0.662745},
myColor {0.929644, 0.00127, 0.145367, 0.666667},
myColor {0.932737, 0.00127, 0.140417, 0.670588},
myColor {0.935747, 0.00127, 0.13544, 0.67451},
myColor {0.938675, 0.00127, 0.130438, 0.678431},
myColor {0.941521, 0.00127, 0.125409, 0.682353},
myColor {0.944285, 0.00127, 0.120354, 0.686275},
myColor {0.946965, 0.00127, 0.115272, 0.690196},
myColor {0.949562, 0.00127, 0.110164, 0.694118},
myColor {0.952075, 0.00127, 0.105031, 0.698039},
myColor {0.954506, 0.00127, 0.099874, 0.701961},
myColor {0.956852, 0.00127, 0.094695, 0.705882},
myColor {0.959114, 0.00127, 0.089499, 0.709804},
myColor {0.961293, 0.00127, 0.084289, 0.713725},
myColor {0.963387, 0.00127, 0.079073, 0.717647},
myColor {0.965397, 0.00127, 0.073859, 0.721569},
myColor {0.967322, 0.00127, 0.068659, 0.72549},
myColor {0.969163, 0.00127, 0.063488, 0.729412},
myColor {0.970919, 0.00127, 0.058367, 0.733333},
myColor {0.97259, 0.00127, 0.053324, 0.737255},
myColor {0.974176, 0.00127, 0.048392, 0.741176},
myColor {0.975677, 0.00127, 0.043618, 0.745098},
myColor {0.977092, 0.00127, 0.03905, 0.74902},
myColor {0.978422, 0.00127, 0.034931, 0.752941},
myColor {0.979666, 0.00127, 0.031409, 0.756863},
myColor {0.980824, 0.00127, 0.028508, 0.760784},
myColor {0.981895, 0.00127, 0.02625, 0.764706},
myColor {0.982881, 0.00127, 0.024661, 0.768627},
myColor {0.983779, 0.00127, 0.02377, 0.772549},
myColor {0.984591, 0.00127, 0.023606, 0.776471},
myColor {0.985315, 0.00127, 0.024202, 0.780392},
myColor {0.985952, 0.00127, 0.025592, 0.784314},
myColor {0.986502, 0.00127, 0.027814, 0.788235},
myColor {0.986964, 0.00127, 0.030908, 0.792157},
myColor {0.987337, 0.00127, 0.034916, 0.796078},
myColor {0.987622, 0.00127, 0.039886, 0.8},
myColor {0.987819, 0.00127, 0.045581, 0.803922},
myColor {0.987926, 0.00127, 0.05175, 0.807843},
myColor {0.987945, 0.00127, 0.058329, 0.811765},
myColor {0.987874, 0.00127, 0.065257, 0.815686},
myColor {0.987714, 0.00127, 0.072489, 0.819608},
myColor {0.987464, 0.00127, 0.07999, 0.823529},
myColor {0.987124, 0.00127, 0.087731, 0.827451},
myColor {0.986694, 0.00127, 0.095694, 0.831373},
myColor {0.986175, 0.00127, 0.103863, 0.835294},
myColor {0.985566, 0.00127, 0.112229, 0.839216},
myColor {0.984865, 0.00127, 0.120785, 0.843137},
myColor {0.984075, 0.00127, 0.129527, 0.847059},
myColor {0.983196, 0.00127, 0.138453, 0.85098},
myColor {0.982228, 0.00127, 0.147565, 0.854902},
myColor {0.981173, 0.00127, 0.156863, 0.858824},
myColor {0.980032, 0.00127, 0.166353, 0.862745},
myColor {0.978806, 0.00127, 0.176037, 0.866667},
myColor {0.977497, 0.00127, 0.185923, 0.870588},
myColor {0.976108, 0.00127, 0.196018, 0.87451},
myColor {0.974638, 0.00127, 0.206332, 0.878431},
myColor {0.973088, 0.00127, 0.216877, 0.882353},
myColor {0.971468, 0.00127, 0.227658, 0.886275},
myColor {0.969783, 0.00127, 0.238686, 0.890196},
myColor {0.968041, 0.00127, 0.249972, 0.894118},
myColor {0.966243, 0.00127, 0.261534, 0.898039},
myColor {0.964394, 0.00127, 0.273391, 0.901961},
myColor {0.962517, 0.00127, 0.285546, 0.905882},
myColor {0.960626, 0.00127, 0.29801, 0.909804},
myColor {0.95872, 0.00127, 0.31082, 0.913725},
myColor {0.956834, 0.00127, 0.323974, 0.917647},
myColor {0.954997, 0.00127, 0.337475, 0.921569},
myColor {0.953215, 0.00127, 0.351369, 0.92549},
myColor {0.951546, 0.00127, 0.365627, 0.929412},
myColor {0.950018, 0.00127, 0.380271, 0.933333},
myColor {0.948683, 0.00127, 0.395289, 0.937255},
myColor {0.947594, 0.00127, 0.410665, 0.941176},
myColor {0.946809, 0.00127, 0.426373, 0.945098},
myColor {0.946392, 0.00127, 0.442367, 0.94902},
myColor {0.946403, 0.00127, 0.458592, 0.952941},
myColor {0.946903, 0.00127, 0.47497, 0.956863},
myColor {0.947937, 0.00127, 0.491426, 0.960784},
myColor {0.949545, 0.00127, 0.50786, 0.964706},
myColor {0.95174, 0.00127, 0.524203, 0.968627},
myColor {0.954529, 0.00127, 0.540361, 0.972549},
myColor {0.957896, 0.00127, 0.556275, 0.976471},
myColor {0.961812, 0.00127, 0.571925, 0.980392},
myColor {0.966249, 0.00127, 0.587206, 0.984314},
myColor {0.971162, 0.00127, 0.602154, 0.988235},
myColor {0.976511, 0.00127, 0.61676, 0.992157},
myColor {0.982257, 0.00127, 0.631017, 0.996078},
myColor {0.988362, 0.00127, 0.644924, 1}
};








#include "haloFitter.inl"

#endif // !LC_HALOFITTER_H

