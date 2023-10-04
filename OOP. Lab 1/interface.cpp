#include "interface.h"

void file_output_distribution(int n, vector<double>& x_s, HuberDistribution* HB)
{
	ofstream xs;
	ofstream fs_theoretical;
	ofstream fs_empirical;

	xs.open("xs.txt");
	fs_theoretical.open("fs_theoretical.txt");
	fs_empirical.open("fs_empirical.txt");

	sort(x_s.begin(), x_s.end());

	for (double& x : x_s)
	{
		xs << x << endl;

		double f_theoretical = Huber(x, HB);
		fs_theoretical << f_theoretical << endl;

		double f_empirical = empirical_huber(n, x, x_s);
		fs_empirical << f_empirical << endl;
	}

	xs.close();
	fs_theoretical.close();
	fs_empirical.close();
}

void file_output_mixture(int n, vector<double>& x_s, Mixture* M) 
{
	ofstream xs;
	ofstream fs_theoretical;
	ofstream fs_empirical;

	xs.open("xs.txt");
	fs_theoretical.open("fs_theoretical.txt");
	fs_empirical.open("fs_empirical.txt");

	sort(x_s.begin(), x_s.end());

	for (double& x : x_s)
	{
		xs << x << endl;

		double f_theoretical = mixture(x, M);
		fs_theoretical << f_theoretical << endl;

		double f_empirical = empirical_huber(n, x, x_s);
		fs_empirical << f_empirical << endl;
	}

	xs.close();
	fs_theoretical.close();
	fs_empirical.close();
}

void general_distribution()
{
	int n = 30000;
	int distribution_params_option;
	int file_option;
	double v, scale, shift, x;
	vector<double> x_s;
	HuberDistribution* HB;

	cout << "Параметры распределения" << endl << "1) Стандартное распределение" << endl << "2) Другие" << endl;
	cin >> distribution_params_option;

	if (distribution_params_option == 1)
	{
		cout << "Введите v: ";
		cin >> v;
		HB = init_huber_distribution(v, K(v));
		x = 0;
		
	}
	else
	{
		cout << "Введите v, lambda, mu, x: ";
		cin >> v >> scale >> shift >> x;
		HB = init_huber_distribution(v, K(v), scale, shift);

	}

	x_s = generate_sequence(n, HB);

	cout << "Для основного распределения:" << endl;
	cout << "Математическое ожидание: " << huber_expected_value(HB) << endl;
	cout << "Дисперсия: " << huber_variance(HB) << endl;
	cout << "Коэффициент асимметрии: " << huber_asymmetry(HB) << endl;
	cout << "Коэффициент эксцесса: " << huber_kurtosis(HB) << endl;
	cout << "Значение плотности в точке " << x << ": " << Huber(x, HB) << endl << endl;

	cout << "Для эмпирического распределения:" << endl;
	cout << "Математическое ожидание: " << empirical_expected_value(n, x_s) << endl;
	cout << "Дисперсия: " << empirical_variance(n, x_s) << endl;
	cout << "Коэффициент асимметрии: " << empirical_asymmetry(n, x_s) << endl;
	cout << "Коэффициент эксцесса: " << empirical_kurtosis(n, x_s) << endl;
	cout << "Значение плотности в точке " << x << ": " << empirical_huber(n, x, x_s) << endl;

	cout << "Вывести данные распределения в отдельный файл? 1 - Да, 0 - Нет: ";
	cin >> file_option;
	if (file_option)
	{
		file_output_distribution(n, x_s, HB);
	}
}

void mixture_distribution()
{
	int n = 1000;
	int file_option;
	double v1, scale1, shift1, v2, scale2, shift2, x, p;
	vector<double> x_s;
	Mixture* M;

	cout << "Введите v1, lambda1, mu1, v2, lambda2, mu2, p, x: ";
	cin >> v1 >> scale1 >> shift1 >> v2 >> scale2 >> shift2 >> p >> x;
	M = init_mixture(p, v1, v2, scale1, scale2, shift1, shift2);

	x_s = generate_sequence(n, M->HB1);

	cout << "Для смеси распределений: " << endl;
	cout << "Математическое ожидание: " << mixture_expected_value(M) << endl;
	cout << "Дисперсия: " << mixture_variance(M) << endl;
	cout << "Коэффициент асимметрии: " << mixture_asymmetry(M) << endl;
	cout << "Коэффициент эксцесса: " << mixture_kurtosis(M) << endl;
	cout << "Значение плотности в точке " << x << ": " << mixture(x, M) << endl << endl;

	cout << "Для эмпирического распределения:" << endl;
	cout << "Математическое ожидание: " << empirical_expected_value(n, x_s) << endl;
	cout << "Дисперсия: " << empirical_variance(n, x_s) << endl;
	cout << "Коэффициент асимметрии: " << empirical_asymmetry(n, x_s) << endl;
	cout << "Коэффициент эксцесса: " << empirical_kurtosis(n, x_s) << endl;
	cout << "Значение плотности в точке " << x << ": " << empirical_huber(n, x, x_s) << endl;

	cout << "Вывести данные распределения в отдельный файл? 1 - Да, 0 - Нет: ";
	cin >> file_option;
	if (file_option)
	{
		file_output_mixture(n, x_s, M);
	}
}


void interface()
{
	int distribution_option;
	int file_option;

	cout << "С каким распределением работаем?" << endl << "1) Основное распределение" << endl << "2) Смесь распределений" << endl << "3) Выйти" << endl;
	cin >> distribution_option;

	switch (distribution_option)
	{
		case 1: 
			general_distribution();
			break;
		case 2: 
			mixture_distribution();
			break;
		default: 
			return;
	}
}