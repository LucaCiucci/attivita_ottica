#ifndef LC_PARABOLA_FIT_H
#define LC_PARABOLA_FIT_H


std::vector<Vector3<>> pts;
double chi2(double a, double x_c, double y_c)
{
	double chi2_v = 0;
	double min_v = INFINITY;
	//double maxVal = 20;
	for (int i = 0; i < pts.size(); i++)
	{
		Vector3<> tmp = pts[i];
		chi2_v += sqr((double)tmp.y - (a * sqr(tmp.x - x_c) + y_c));
		//chi2_v += sqr(atan(((double)tmp.y - (a * sqr(tmp.x - x_c) + y_c)) / (maxVal)) * maxVal);// TODO rimetti quello vero sopra
		//chi2_v += clamp(sqr((double)tmp.y - (a * sqr(tmp.x - x_c) + y_c)), 0, sqr(maxVal));// TODO rimetti quello vero sopra
	}
	return chi2_v / pts.size();
}

double parabolaFit(std::vector<Vector3<>> & points, int minOffset = 10)// ma con minOffset = 10 va ancora meglio
{
	double x_c, y_c, a;
	sf::RenderWindow window(sf::VideoMode(800, 600), "parabola fitter    @ Luca Ciucci");
	window.setPosition(sf::Vector2i(sf::VideoMode::getDesktopMode().width - 810, 10));
	Vector3<> vmin(INFINITY), vmax(-INFINITY);
	std::vector<int> histo(256);
	for (int i = 0; i < points.size(); i++)
	{
		histo[(int)points[i].y]++;
	}
	int maxHistoIndex = 0;
	{
		int maxIndex;
		for (maxIndex = 255; maxIndex >= 0; maxIndex--)// NOTA se parte da 128 ad esempio becca 255=normale
		{
			if (histo[maxIndex] > /*0*/4)
				break;
		}
		int max_v = 0;
		for (int i = maxIndex; i > 0 && i > maxIndex - minOffset; i--)
		{
			if (histo[i] > max_v)
			{
				maxHistoIndex = i;
				max_v = histo[i];
			}
		}
	}

	pts.clear();
	for (int i = 0; i < points.size(); i++)
	{
		if (points[i].y + minOffset > maxHistoIndex)
		{
			pts.push_back(points[i]);
			vmax.x = fmax(vmax.x, points[i].x);
			vmax.y = fmax(vmax.y, points[i].y);
			vmin.x = fmin(vmin.x, points[i].x);
			vmin.y = fmin(vmin.y, points[i].y);
		}
	}
	x_c = (vmax + vmin).x * 0.5;
	y_c = vmax.y;
	a = (vmin - vmax).y / sqr(x_c - vmin.x) * 0.5;

	FunctionOptimizer<> fo(chi2);
	fo.params[0] = a;
	fo.params[1] = x_c;
	fo.params[2] = y_c;
	//fo.oneStep(0.0005, 200001);

	sf::Event event;
	while (window.isOpen())
	{
		window.clear();


		double r = 2.0;
		sf::CircleShape dot(r);
		dot.setOrigin(sf::Vector2f(r, r));
		dot.setFillColor(sf::Color::Red);
		for (int i = 0; i < pts.size(); i++)
		{
			dot.setPosition(sf::Vector2f((pts[i] - vmin).x / (vmax - vmin).x * window.getSize().x, (1.0 - (pts[i] - vmin).y / (vmax - vmin).y) * window.getSize().y));
			window.draw(dot);
		}

		fo.oneStep(0.0005, 10001);

		int N = 1000;
		sf::VertexArray lines(sf::LinesStrip, N);
		double xStart = vmin.x;
		for (int i = 0; i < N; i++)
		{
			double x = xStart + (vmax - vmin).x * (double)i / (double)N;
			double y = fo.params[0] * sqr(x - fo.params[1]) + fo.params[2];
			lines[i].position.x = (x - vmin.x) / (vmax - vmin).x * window.getSize().x;
			lines[i].position.y = (1.0 - (y - vmin.y) / (vmax - vmin).y) * window.getSize().y;
			lines[i].color = sf::Color::Green;
		}
		window.draw(lines);

		window.display();

		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::KeyPressed)
			{
				switch (event.key.code)
				{
				case sf::Keyboard::Y:
					// finito
					return fo.params[1];
					break;
				case sf::Keyboard::N:
					return 0.0;
					break;
				default:
					break;
				}
			}
		}
	}
	return 0.0;
}

#endif // !LC_PARABOLA_FIT_H