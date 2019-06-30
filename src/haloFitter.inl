#include "haloFitter.h"
#include "parabola_fit.h"
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
#ifndef LC_HALOFITTER_INL
#define LC_HALOFITTER_INL

#include "inverseMatrix.h"


// ================================================================
//                        HALO FITTER
// ================================================================

// ================================
//            CONSTRUCTOR
// ================================

template<class LetVoid>
inline HaloFitter<LetVoid>::HaloFitter(std::string imgName):
	window(sf::VideoMode(800, 600), "halo fitter    @ Luca Ciucci"),
	maskColor(sf::Color::Cyan.r, sf::Color::Cyan.g, sf::Color::Cyan.b, 128), minLevel(20), chi2MaskSize(10), chi2MaskDivider(8)
{
	std::cout << "Starting..." << std::endl;
	
	// TODO sposta
	dr = 10.0;
	zoomFactor = 1.125;
	currZoom = 1.0;
	colorZoom = 1.0;
	showCircle = true;

	// load inìmage
	sf::Image background;
	if (!background.loadFromFile(imgName))
	{
		std::cout << "E: unable to load image!" << std::endl;
		exit(-1);
	}
	std::cout << imgName << " loaded!" << std::endl;

	wsX = background.getSize().x;
	wsY = background.getSize().y;
	window.setView(sf::View(sf::Vector2f(wsX * 0.5, wsY * 0.5), sf::Vector2f(wsX, wsY)));
	window.setSize(sf::Vector2u(wsX, wsY));
	window.setPosition(sf::Vector2i(50, 50));
	initialView = window.getView();

	// crea immagine
	if (!img_texture.create(wsX, wsY))
		exit(-1);
	if (!chi2_texture.create(chi2MaskSize * chi2MaskDivider, chi2MaskSize * chi2MaskDivider))
		exit(-1);
	img_pixels = new sf::Uint8[wsX * wsY * 4];
	tmp_pixels = new sf::Uint8[wsX * wsY * 4];
	mask_pixels = new sf::Uint8[wsX * wsY * 4]; memset(mask_pixels, 0, sizeof(sf::Uint8) * wsX * wsY * 4);
	chi2Mask_pixels = new sf::Uint8[sqr(chi2MaskSize * chi2MaskDivider) * 4]; memset(chi2Mask_pixels, 0, sizeof(sf::Uint8) * sqr(chi2MaskSize * chi2MaskDivider) * 4);
	memcpy(tmp_pixels, img_pixels, sizeof(sf::Uint8) * wsX * wsY * 4);
	img_sprite.setTexture(img_texture);
	img_sprite.setScale(sf::Vector2f(1.0, 1.0));
	img_sprite.setPosition(sf::Vector2f(0.0, 0.0));
	chi2_sprite.setTexture(chi2_texture);
	chi2_sprite.setScale(sf::Vector2f(1.0 / chi2MaskDivider, 1.0 / chi2MaskDivider));
	chi2_sprite.setPosition(sf::Vector2f(0.0, 0.0));
	chi2_sprite.setOrigin(sf::Vector2f(chi2MaskSize * chi2MaskDivider * 0.5, chi2MaskSize * chi2MaskDivider * 0.5));
	for (int y = 0; y < wsY; y++)
	{
		for (int x = 0; x < wsX; x++)
		{
			img_pixels[(y * wsX + x) * 4 + 0] = background.getPixel(x, y).r;
			img_pixels[(y * wsX + x) * 4 + 1] = background.getPixel(x, y).g;
			img_pixels[(y * wsX + x) * 4 + 2] = background.getPixel(x, y).b;
			img_pixels[(y * wsX + x) * 4 + 3] = background.getPixel(x, y).a;
		}
	}

	// crea finestra e mostra immagine
	while (window.pollEvent(event)) {/* clear events */ };
	plotComputerView();
	std::cout << "Started!" << std::endl;
	// TODO remove
	for (int i = 0; i < 0*200; i++)
	{
		while (window.pollEvent(event)) {/* clear events */ };
		lc::system::sleepSec(0.01);
	}
}


// ================================
//            FUNCTIONS
// ================================

template<class LetVoid>
inline void HaloFitter<LetVoid>::getUserPoints(void)
{
	memcpy(tmp_pixels, img_pixels, sizeof(sf::Uint8) * wsX * wsY * 4);
	while (window.isOpen())
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
			{
				window.close();
				std::cout << "Not supposed to do That! AHAH funny!" << std::endl;
				exit(-1);
			}
			if (event.type == sf::Event::MouseButtonPressed)
			{
				switch (event.mouseButton.button)
				{
				case sf::Mouse::Left:
				{
					sf::Vector2i pixelPos = sf::Mouse::getPosition(window);
					sf::Vector2f worldPos = window.mapPixelToCoords(pixelPos);
					userPts.push_back(Vector3<>(worldPos.x, worldPos.y, 0.0));
					updateUserPts();
					break;
				}
				case sf::Mouse::Middle:
				{
					// then pan
					sf::View currview = window.getView();
					sf::Vector2i startPos(sf::Mouse::getPosition());
					while (true)
					{
						sf::Vector2i mousePos = sf::Mouse::getPosition();
						sf::Vector2f realShift((double)(mousePos - startPos).x / wsX * window.getView().getSize().x, (double)(mousePos - startPos).y / wsY * window.getView().getSize().y);
						window.setView(sf::View(currview.getCenter() - realShift, currview.getSize()));
						plotUserView();
						while (window.pollEvent(event))
						{
							if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Middle)
								break;
						}
						if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Middle)
							break;
					}
					break;
				}
				default:
					break;
				}
			}
			if (event.type == sf::Event::KeyPressed)
			{
				switch (event.key.code)
				{
				case sf::Keyboard::Escape:
					// finito
					return;
					break;
				case sf::Keyboard::F1:
					colorZoom = 1.0;
					memcpy(tmp_pixels, img_pixels, sizeof(sf::Uint8) * wsX * wsY * 4);
					break;
				case sf::Keyboard::Add:
				{
					colorZoom += 0.1;
					for (int i = 0; i < wsX * wsY * 4; i++)
					{
						tmp_pixels[i] = clamp((double)colorZoom * (double)img_pixels[i], 0.0f, 255.0f);
					}
					break;
				}
				case sf::Keyboard::Subtract:
				{
					colorZoom -= 0.1;
					for (int i = 0; i < wsX * wsY * 4; i++)
					{
						tmp_pixels[i] = clamp((double)colorZoom *(double)img_pixels[i], 0.0f, 255.0f);
					}
					break;
				}
				case sf::Keyboard::F2:
				{
					showCircle = !showCircle;
					break;
				}
				default:
					break;
				}
			}
			if (event.type == sf::Event::MouseWheelScrolled)
			{
				sf::View currview = window.getView();
				sf::Vector2f currViewSize = currview.getSize();
				sf::Vector2f currViewCenter = currview.getCenter();
				sf::Vector2i mousePos = sf::Mouse::getPosition(window);

				double factor = pow(zoomFactor, -event.mouseWheelScroll.delta);
				currViewSize.x *= factor;
				currViewSize.y *= factor;
				currZoom *= factor;
				
				// center on mouse
				sf::Vector2i pixelPos = sf::Mouse::getPosition(window);
				sf::Vector2f worldPos = window.mapPixelToCoords(pixelPos);
				currViewCenter.x -= (worldPos - currViewCenter).x * (factor - 1.0);
				currViewCenter.y -= (worldPos - currViewCenter).y * (factor - 1.0);
				
				window.setView(sf::View(currViewCenter, currViewSize));
			}

			plotUserView();
		}
}

template<class LetVoid>
inline void HaloFitter<LetVoid>::plotUserView(void)
{
	window.clear();

	// plot image
	img_texture.update(tmp_pixels);
	img_texture.setSmooth(true);
	sf::Color showColor(sf::Color::Black);
	if (sf::Keyboard::isKeyPressed(sf::Keyboard::R))
		showColor += sf::Color::Red;
	if (sf::Keyboard::isKeyPressed(sf::Keyboard::G))
		showColor += sf::Color::Green;
	if (sf::Keyboard::isKeyPressed(sf::Keyboard::B))
		showColor += sf::Color::Blue;
	img_sprite.setColor(showColor == sf::Color::Black ? sf::Color::White : showColor);
	window.draw(img_sprite);
	
	// plot user points
	double r = 3.0 * currZoom;
	sf::CircleShape dot(r);
	dot.setFillColor(sf::Color::Red);
	dot.setOrigin(r, r);
	for (int i = 0; i < userPts.size(); i++)
	{
		dot.setPosition(sf::Vector2f(userPts[i].x, userPts[i].y));
		window.draw(dot);
	}

	if (showCircle)
	{
		// plot use circle
		sf::CircleShape center(r);
		center.setFillColor(sf::Color::Green);
		center.setOrigin(r, r);
		center.setPosition(sf::Vector2f(UsrFitCircle.x, UsrFitCircle.y));
		sf::CircleShape circle(UsrFitCircle.z, 1000);
		circle.setFillColor(sf::Color::Transparent);
		circle.setOutlineColor(sf::Color::Green);
		circle.setOutlineThickness(currZoom);
		circle.setOrigin(UsrFitCircle.z, UsrFitCircle.z);
		circle.setPosition(sf::Vector2f(UsrFitCircle.x, UsrFitCircle.y));
		window.draw(circle);
		window.draw(center);
	}

	window.display();
}

template<class LetVoid>
inline void HaloFitter<LetVoid>::plotComputerView(void)
{
	window.clear();

	// plot image
	img_texture.update(tmp_pixels);
	img_texture.setSmooth(false);// TODO false
	sf::Color showColor(sf::Color::Black);
	if (sf::Keyboard::isKeyPressed(sf::Keyboard::R))
		showColor += sf::Color::Red;
	if (sf::Keyboard::isKeyPressed(sf::Keyboard::G))
		showColor += sf::Color::Green;
	if (sf::Keyboard::isKeyPressed(sf::Keyboard::B))
		showColor += sf::Color::Blue;
	img_sprite.setColor(showColor == sf::Color::Black ? sf::Color::White : showColor);
	window.draw(img_sprite);

	// plot mask
	img_texture.update(mask_pixels);
	img_texture.setSmooth(false);
	window.draw(img_sprite);

	// plot computer points
	double r = 3.0 * currZoom;
	sf::CircleShape dot(r);
	dot.setFillColor(sf::Color::Red);
	dot.setOrigin(r, r);
	for (int i = 0; i < computerPts.size(); i++)
	{
		dot.setPosition(sf::Vector2f(computerPts[i].x, computerPts[i].y));
		window.draw(dot);
	}

	// chi2 mask
	if (computerPts.size() > 3)
	{
		chi2_texture.update(chi2Mask_pixels);
		chi2_texture.setSmooth(true);
		chi2_sprite.setPosition(sf::Vector2f(computerFitCircle.x, computerFitCircle.y));
		window.draw(chi2_sprite);
	}

	// plot use circle
	sf::CircleShape center(r);
	center.setFillColor(sf::Color::Green);
	center.setOrigin(r, r);
	center.setPosition(sf::Vector2f(computerFitCircle.x, computerFitCircle.y));
	sf::CircleShape circle(computerFitCircle.z, 1000);
	circle.setFillColor(sf::Color::Transparent);
	circle.setOutlineColor(sf::Color::Green);
	circle.setOutlineThickness(currZoom);
	circle.setOrigin(computerFitCircle.z, computerFitCircle.z);
	circle.setPosition(sf::Vector2f(computerFitCircle.x, computerFitCircle.y));
	window.draw(circle);
	window.draw(center);


	window.display();
}

template<class LetVoid>
inline void HaloFitter<LetVoid>::updateUserPts(void)
{
	UsrFitCircle = fastCircleFit(userPts);
}

template<class LetVoid>
inline void HaloFitter<LetVoid>::refine(int color, int nStep)
{
	window.setView(initialView); currZoom = 1.0;
	double dAlpha = acos(-1.0) / nStep; // note 1/2!
	double cosTolerance = cos(dAlpha);
	Vector3<> falseCenter(UsrFitCircle.x, UsrFitCircle.y, 0.0);
	double falseRadius(UsrFitCircle[2]);
	for (int t = 0; t < nStep; t++)
	{
		mask_pixels = new sf::Uint8[wsX * wsY * 4];
		memset(mask_pixels, 0, sizeof(sf::Uint8) * wsX * wsY * 4);
		selectedPxs.clear();
		double alpha = 2*acos(-1.0) * (double)t / nStep + acos(-1.0);
		Vector3<> p0(cos(alpha) * falseRadius, sin(alpha) * falseRadius, 0.0);
		//p0 += falseCenter;

		//std::cout << m_PI<> << std::endl;// TODO m_pi<> a volte non funzion?!?!?!?!?! non ha senso, bug di vs!!!!

		// find selected pixels
		for (int y = 0; y < wsY; y++)
		{
			for (int x = 0; x < wsX; x++)
			{
				sf::Vector2i pixelPos(x, y);
				sf::Vector2f worldPos = window.mapPixelToCoords(pixelPos);
				Vector3<> p(worldPos.x, worldPos.y, 0.0);
				p -= falseCenter;
				if (isIninterval(p.norm(), falseRadius - dr, falseRadius + dr) && (p0 * p) / (p0.norm() * p.norm()) > cosTolerance)
				//if ((p*p0) / (p.norm() * p0.norm()) > cosTolerance)
				{
					if (img_pixels[(y * wsX + x) * 4 + color] > minLevel)
						selectedPxs.push_back(Vector3<>(p.norm(), img_pixels[(y * wsX + x) * 4 + color], 0.0));
					mask_pixels[(y * wsX + x) * 4 + 0] = maskColor.r;
					mask_pixels[(y * wsX + x) * 4 + 1] = maskColor.g;
					mask_pixels[(y * wsX + x) * 4 + 2] = maskColor.b;
					mask_pixels[(y * wsX + x) * 4 + 3] = maskColor.a;
				}
			}
		}

		//
		plotComputerView();
		double r = parabolaFit(selectedPxs);
		if (r > 0)
		{
			computerPts.push_back(Vector3<>(cos(alpha) * r, sin(alpha) * r, 0.0) + falseCenter);
			computerFitCircle = fastCircleFit(computerPts);
			// TODO riattiva
			//if (computerPts.size() > 3)
				//computeChi2Map();
		}

		while (window.pollEvent(event)) { /* clear events*/ }
		lc::system::sleepSec(0.05);
	}
}

template<class LetVoid>
inline double& HaloFitter<LetVoid>::setDRadius(double new_dr)
{
	return dr = new_dr;
}

template<class LetVoid>
inline void HaloFitter<LetVoid>::showResult(void)
{
	/*hessianMatrix(computerFitCircle.x, computerFitCircle.y, computerFitCircle.z, 0.0001);
	hessianMatrix(computerFitCircle.x, computerFitCircle.y, computerFitCircle.z, 0.001);
	hessianMatrix(computerFitCircle.x, computerFitCircle.y, computerFitCircle.z, 0.01);
	hessianMatrix(computerFitCircle.x, computerFitCircle.y, computerFitCircle.z, 0.1);
	hessianMatrix(computerFitCircle.x, computerFitCircle.y, computerFitCircle.z, 1.0);*/
	computeCovMatrix();
	window.setView(initialView); currZoom = 1.0;
	memcpy(tmp_pixels, img_pixels, sizeof(sf::Uint8) * wsX * wsY * 4);
	memset(mask_pixels, 0, sizeof(sf::Uint8) * wsX * wsY * 4);
	computeChi2Map();
	while (window.isOpen())
	{
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
			if (event.type == sf::Event::MouseWheelScrolled)
			{
				sf::View currview = window.getView();
				sf::Vector2f currViewSize = currview.getSize();
				sf::Vector2f currViewCenter = currview.getCenter();
				sf::Vector2i mousePos = sf::Mouse::getPosition(window);

				double factor = pow(zoomFactor, -event.mouseWheelScroll.delta);
				currViewSize.x *= factor;
				currViewSize.y *= factor;
				currZoom *= factor;

				// center on mouse
				sf::Vector2i pixelPos = sf::Mouse::getPosition(window);
				sf::Vector2f worldPos = window.mapPixelToCoords(pixelPos);
				currViewCenter.x -= (worldPos - currViewCenter).x * (factor - 1.0);
				currViewCenter.y -= (worldPos - currViewCenter).y * (factor - 1.0);

				window.setView(sf::View(currViewCenter, currViewSize));
			}
			if (event.type == sf::Event::MouseButtonPressed)
			{
				switch (event.mouseButton.button)
				{
				case sf::Mouse::Middle:
				{
					// then pan
					sf::View currview = window.getView();
					sf::Vector2i startPos(sf::Mouse::getPosition());
					while (true)
					{
						sf::Vector2i mousePos = sf::Mouse::getPosition();
						sf::Vector2f realShift((double)(mousePos - startPos).x / wsX * window.getView().getSize().x, (double)(mousePos - startPos).y / wsY * window.getView().getSize().y);
						window.setView(sf::View(currview.getCenter() - realShift, currview.getSize()));
						plotComputerView();
						while (window.pollEvent(event))
						{
							if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Middle)
								break;
						}
						if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Middle)
							break;
					}
					break;
				}
				default:
					break;
				}
			}
			if (event.type == sf::Event::KeyPressed)
			{
				switch (event.key.code)
				{
				case sf::Keyboard::Escape:
					window.close();
					break;
				case sf::Keyboard::Add:
				{
					colorZoom += 0.1;
					for (int i = 0; i < wsX * wsY * 4; i++)
					{
						tmp_pixels[i] = clamp((double)colorZoom * (double)img_pixels[i], 0.0f, 255.0f);
					}
					break;
				}
				case sf::Keyboard::Subtract:
				{
					colorZoom -= 0.1;
					for (int i = 0; i < wsX * wsY * 4; i++)
					{
						tmp_pixels[i] = clamp((double)colorZoom * (double)img_pixels[i], 0.0f, 255.0f);
					}
					break;
				}
				default:
					break;
				}
			}
		}
		plotComputerView();
	}
	std::cout << "params: " << computerFitCircle << std::endl;
}

template<class LetVoid>
inline void HaloFitter<LetVoid>::plotRes(void)
{
	/*int minR, maxR;
	std::cout << "minR = " << std::endl;
	std::cin >> minR;
	std::cout << "MAXR = " << std::endl;
	std::cin >> maxR;
	int myWs = 500;
	sf::RenderWindow nWindow(sf::VideoMode(myWs, myWs));
	sf::Texture plot_textue;
	if (!plot_textue.create(myWs, myWs))
		exit(-1);
	sf::Uint8* img_pixels = new sf::Uint8[myWs * myWs * 4];
	sf::Sprite plot_sprite.setTexture(plot_texture);
	//plot_sprite.setScale(sf::Vector2f(1.0, 1.0));
	//plot_sprite.setPosition(sf::Vector2f(0.0, 0.0));

	for (int y = 0; y < myWs; y++)
		for (int x = 0; x < myWs; x++)
		{

		}*/
	for (int i = 0; i < computerPts.size(); i++)
	{
		Vector3<> center(computerFitCircle.x, computerFitCircle.y, 0.0);
		Vector3<> diff = computerPts[i] - center;
		double R = (diff).norm();
		double alpha = atan2(diff.y, diff.x);
		std::cout << alpha << " " << R << std::endl;
	}

	std::cout << "\n\n ora i punti in x-y:" << std::endl;
	for (int i = 0; i < computerPts.size(); i++)
	{
		Vector3<> center(computerFitCircle.x, computerFitCircle.y, 0.0);
		Vector3<> diff = computerPts[i] - center;
		double R = (diff).norm();
		double alpha = atan2(diff.y, diff.x);
		std::cout << computerPts[i] << std::endl;
	}
}

std::vector<Vector3<>> tmpPts;
Vector3<> tmpCircleCenter;
inline double chi2_on_r(double R)
{
	double sum_v = 0.0;
	for (int i = 0; i < tmpPts.size(); i++)
	{
		//sum_v += sqr((tmpPts[i] - tmpCircleCenter).norm2() - R * R);
		sum_v += sqr((tmpPts[i] - tmpCircleCenter).norm() - R);
	}
	return sum_v;
}

inline double chi2_on_all_parameters(double xc, double yc, double r)
{
	double sum_v = 0.0;
	for (int i = 0; i < tmpPts.size(); i++)
	{
		sum_v += sqr((tmpPts[i] - Vector3<>(xc, yc, 0.0)).norm() - r);
		//std::cout << "punti" << tmpPts[i] << std::endl;
		//std::cout << "centro" << Vector3<>(xc, yc, 0.0) << std::endl;
		//std::cout << "c2 " << sqr((tmpPts[i] - Vector3<>(xc, yc, 0.0)).norm2() - r * r) << std::endl;
		//std::cout << "c2 " << sqr((tmpPts[i] - Vector3<>(xc, yc, 0.0)).norm()) << std::endl;
	}
	return sum_v;
}


template<class LetVoid>
inline void HaloFitter<LetVoid>::hessianMatrix(double xc, double yc, double r, double dh)
{
	std::cout << "chi2 da analitica: " << chi2_on_all_parameters(xc, yc, r) << std::endl;
	VectorN<3> startParams, newParams; startParams.x = xc; startParams.y = yc; startParams.z = r;
	newParams = lc::math::OptimizeFunction(chi2_on_all_parameters, startParams, 0.1, 1000);
	xc = newParams.x;
	yc = newParams.y;
	r = newParams.z;
	std::cout << "chi2 da numerica: " << chi2_on_all_parameters(xc, yc, r) << std::endl;
	std::cout << "parametri:" << startParams << newParams << newParams - startParams << std::endl;
	computerFitCircle = Vector3<>(xc, yc, r);
	std::cout << xc << " " << yc << " " << r << std::endl;
	std::cout << chi2_on_all_parameters(xc, yc, r) << std::endl;
	std::cout << "derivata xc: " << (chi2_on_all_parameters(xc + dh, yc, r) - chi2_on_all_parameters(xc, yc, r)) / dh << std::endl;
	std::cout << "derivata yc: " << (chi2_on_all_parameters(xc, yc + dh, r) - chi2_on_all_parameters(xc, yc, r)) / dh << std::endl;
	std::cout << "derivata r : " << (chi2_on_all_parameters(xc, yc, r + dh) - chi2_on_all_parameters(xc, yc, r)) / dh << std::endl;


	// tmp
	double minChi2 = chi2_on_all_parameters(xc, yc, r);
	//double maxChi2 = minChi2 * ((double)computerPts.size() + sqrt(2.0 * computerPts.size())) / computerPts.size();
	double scaleFactor = minChi2 / (computerPts.size() - 3.0);
	double v0 = minChi2, v1, v2;

	lc::math::VectorN<3> p; p.x = xc; p.y = yc; p.z = r;
	lc::math::Mat<3, 3> mat = lc::math::getHessianMatrix(chi2_on_all_parameters, p);

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			hM[i][j] = (mat.M[i][j] /= scaleFactor);
	std::cout << mat << std::endl;
	return;
	// !tmp

	// rimetti
	//double minChi2 = chi2_on_all_parameters(xc, yc, r);
	//double maxChi2 = minChi2 * ((double)computerPts.size() + sqrt(2 * computerPts.size())) / computerPts.size();
	//double scaleFactor = maxChi2 - minChi2;
	//double v0 = minChi2, v1, v2;

	v1 = chi2_on_all_parameters(xc - dh, yc, r);
	v2 = chi2_on_all_parameters(xc + dh, yc, r);
	hM[0][0] = ((v2 - v0) - (v0 - v1)) / sqr(dh);
	std::cout << "derivata2 xc: " << hM[0][0] << std::endl;

	v1 = chi2_on_all_parameters(xc + dh, yc + dh, r) - chi2_on_all_parameters(xc + dh, yc - dh, r);
	v2 = chi2_on_all_parameters(xc - dh, yc + dh, r) - chi2_on_all_parameters(xc - dh, yc - dh, r);
	hM[1][0] = (v1 - v2) / sqr(2.0*dh);

	v1 = chi2_on_all_parameters(xc, yc - dh, r);
	v2 = chi2_on_all_parameters(xc, yc + dh, r);
	hM[1][1] = ((v1 - v0) - (v0 - v1)) / sqr(dh);
	std::cout << "derivata2 xc: " << hM[1][1] << std::endl;

	v1 = chi2_on_all_parameters(xc + dh, yc, r + dh) - chi2_on_all_parameters(xc + dh, yc, r - dh);
	v2 = chi2_on_all_parameters(xc - dh, yc, r + dh) - chi2_on_all_parameters(xc - dh, yc, r - dh);
	hM[2][0] = (v1 - v2) / sqr(2.0*dh);

	v1 = chi2_on_all_parameters(xc, yc + dh, r + dh) - chi2_on_all_parameters(xc, yc + dh, r - dh);
	v2 = chi2_on_all_parameters(xc, yc - dh, r + dh) - chi2_on_all_parameters(xc, yc - dh, r - dh);
	hM[2][1] = (v1 - v2) / sqr(2.0*dh);

	v1 = chi2_on_all_parameters(xc, yc, r - dh);
	v2 = chi2_on_all_parameters(xc, yc, r + dh);
	hM[2][2] = ((v2 - v0) - (v0 - v1)) / sqr(dh);
	std::cout << "derivata2 xc: " << hM[2][2] << std::endl;

	hM[0][1] = hM[1][0];
	hM[0][2] = hM[2][0];
	hM[1][2] = hM[2][1];

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			hM[i][j] /= scaleFactor;
	printf(hM, 3, 1/*mode*/);
}

template<class LetVoid>
inline void HaloFitter<LetVoid>::computeCovMatrix(void)
{
	tmpPts = computerPts;
	hessianMatrix(computerFitCircle.x, computerFitCircle.y, computerFitCircle.z);
	double deter = det(hM, 3);
	inverse(hM, covM, 3, deter);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			covM[i][j] *= 1.0;
	printf(covM, 3, 2/*mode*/);
}

template<class LetVoid>
inline Vector3<> HaloFitter<LetVoid>::fastCircleFit(std::vector<Vector3<>>& points)
{
	Vector3<> ptsBaricenter = 0.0;
	for (int i = 0; i < points.size(); i++)
	{
		ptsBaricenter += points[i];
	}
	ptsBaricenter /= points.size();

	std::vector<Vector3<>> Pts;
	double Suu = 0, Svv = 0, Suv = 0, Suuu = 0, Svvv = 0, Suuv = 0, Suvv = 0, uc, vc, R;
	for (int i = 0; i < points.size(); i++)
	{
		Pts.push_back(points[i] - ptsBaricenter);

		Suu = Suu + Pts[i].x * Pts[i].x;
		Svv = Svv + Pts[i].y * Pts[i].y;
		Suv = Suv + Pts[i].x * Pts[i].y;
		Suuu = Suuu + Pts[i].x * Pts[i].x * Pts[i].x;
		Svvv = Svvv + Pts[i].y * Pts[i].y * Pts[i].y;
		Suuv = Suuv + Pts[i].x * Pts[i].x * Pts[i].y;
		Suvv = Suvv + Pts[i].x * Pts[i].y * Pts[i].y;
	}

	uc = (Svv * (Suuu + Suvv) - Suv * (Svvv + Suuv)) / (2 * (Suu * Svv - Suv * Suv));
	vc = (Suv * (Suuu + Suvv) - Suu * (Svvv + Suuv)) / (2 * (Suv * Suv - Suu * Svv));
	R = sqrt(uc * uc + vc * vc + (Suu + Svv) / Pts.size());
	double xc = uc + ptsBaricenter.x;
	double yc = vc + ptsBaricenter.y;

	return Vector3<>(xc, yc, R);
}

template<class LetVoid>
inline void HaloFitter<LetVoid>::computeChi2Map(void)
{
	// con 1 sigma
	memset(chi2Mask_pixels, 0, sizeof(sf::Uint8) * sqr(chi2MaskSize * chi2MaskDivider) * 4);
	tmpPts = computerPts;
	Vector3<> computerCenter(computerFitCircle.x, computerFitCircle.y, 0.0);
	tmpCircleCenter = Vector3<>(computerFitCircle.x, computerFitCircle.y, 0.0);
	double minChi2 = chi2_on_r(computerFitCircle.z);
	double maxChi2 = minChi2 * ((double)computerPts.size() + sqrt(2.0 * computerPts.size())) / computerPts.size();
	//double maxChi2 = minChi2 * 10.0;
	for (int y = 0; y < chi2MaskSize * chi2MaskDivider; y++)
	{
		for (int x = 0; x < chi2MaskSize * chi2MaskDivider; x++)
		{
			tmpCircleCenter = computerCenter + (Vector3<>(x, y, 0.0) / (chi2MaskSize * chi2MaskDivider) - 0.5) * chi2MaskSize * 2.0;
			//FunctionOptimizer<> fo(chi2_on_r);
			//fo.params[0] = computerFitCircle.z; fo.oneStep(0.1, 50);
			//double chi2_v = chi2_on_r(fo.params[0]);
			VectorN<1> pOpt(computerFitCircle.z);
			pOpt = lc::math::OptimizeFunction(chi2_on_r, pOpt, 0.5, 50);
			double chi2_v = chi2_on_r(pOpt[0]);
			if (chi2_v < maxChi2)
			{
				//std::cout << "a" << std::endl;
				sf::Color color = toHeatmap(chi2_v, minChi2, maxChi2);
				//chi2Mask_pixels[(y * chi2MaskSize * chi2MaskDivider + x) * 4 + 0] = 255;
				chi2Mask_pixels[(y * chi2MaskSize * chi2MaskDivider + x) * 4 + 0] = color.r;
				chi2Mask_pixels[(y * chi2MaskSize * chi2MaskDivider + x) * 4 + 1] = color.g;
				chi2Mask_pixels[(y * chi2MaskSize * chi2MaskDivider + x) * 4 + 2] = color.b;
				chi2Mask_pixels[(y * chi2MaskSize * chi2MaskDivider + x) * 4 + 3] = 255;// color.a;
			}
		}
	}
}

template<class LetVoid>
inline sf::Color HaloFitter<LetVoid>::toHeatmap(double value, double min, double Max)
{
	int index = (1.0 - (value - min) / (Max - min))*255;

	if (index < 0)
		return sf::Color(0.0, 255, 0.0, 255);
	if (index > 254)
		return sf::Color(0.0, 255, 255, 255);

	return sf::Color(255 - index, 0.0, index, 255);
	
	return sf::Color(inferno[index].r * 255, inferno[index].g * 255, inferno[index].b * 255, inferno[index].a * 255);
}




// ================================
//             PRIVATE
// ================================


#endif // !LC_HALOFITTER_INL

