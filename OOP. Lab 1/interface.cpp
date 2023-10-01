#include "interface.h"

void general_distribution()
{
	int n = 100;
	int distribution_params_option;
	double v, scale, shift, x;
	vector<double> x_s;
	HuberDistribution* HB;

	cout << "��������� �������������" << endl << "1) ����������� �������������" << endl << "2) ������" << endl;
	cin >> distribution_params_option;

	if (distribution_params_option == 1)
	{
		cout << "������� v: ";
		cin >> v;
		HB = init_huber_distribution(v, K(v));
		x = 0;
		
	}
	else
	{
		cout << "������� v, lambda, mu, x: ";
		cin >> v >> scale >> shift >> x;
		HB = init_huber_distribution(v, K(v), scale, shift);

	}

	x_s = generate_sequence(n, HB);

	cout << "��� ��������� �������������:" << endl;
	cout << "�������������� ��������: " << huber_expected_value(HB) << endl;
	cout << "���������: " << huber_variance(HB) << endl;
	cout << "����������� ����������: " << huber_asymmetry(HB) << endl;
	cout << "����������� ��������: " << huber_kurtosis(HB) << endl;
	cout << "�������� ��������� � ����� " << x << ": " << Huber(x, HB) << endl << endl;

	cout << "��� ������������� �������������:" << endl;
	cout << "�������������� ��������: " << empirical_expected_value(n, x_s) << endl;
	cout << "���������: " << empirical_variance(n, x_s) << endl;
	cout << "����������� ����������: " << empirical_asymmetry(n, x_s) << endl;
	cout << "����������� ��������: " << empirical_kurtosis(n, x_s) << endl;
	cout << "�������� ��������� � ����� " << x << ": " << empirical_huber(n, x, x_s) << endl;
}

void mixture_distribution()
{
	int n = 100;
	int mixture_params_option;
	double v1, scale1, shift1, v2, scale2, shift2, x, p;
	vector<double> x_s;
	Mixture* M;

	cout << "������� v1, lambda1, mu1, v2, lambda2, mu2, p, x: ";
	cin >> v1 >> scale1 >> shift1 >> v2 >> scale2 >> shift2 >> p >> x;
	M = init_mixture(p, v1, v2, scale1, scale2, shift1, shift2);

	//����� ������������� �����???
	x_s = generate_sequence(n, M->HB1);

	cout << "��� ����� �������������: " << endl;
	cout << "�������������� ��������: " << mixture_expected_value(M) << endl;
	cout << "���������: " << mixture_variance(M) << endl;
	cout << "����������� ����������: " << mixture_asymmetry(M) << endl;
	cout << "����������� ��������: " << mixture_kurtosis(M) << endl;
	cout << "�������� ��������� � ����� " << x << ": " << mixture(x, M) << endl << endl;

	cout << "��� ������������� �������������:" << endl;
	cout << "�������������� ��������: " << empirical_expected_value(n, x_s) << endl;
	cout << "���������: " << empirical_variance(n, x_s) << endl;
	cout << "����������� ����������: " << empirical_asymmetry(n, x_s) << endl;
	cout << "����������� ��������: " << empirical_kurtosis(n, x_s) << endl;
	cout << "�������� ��������� � ����� " << x << ": " << empirical_huber(n, x, x_s) << endl;
}


void interface()
{
	int distribution_option;

	cout << "� ����� �������������� ��������?" << endl << "1) �������� �������������" << endl << "2) ����� �������������" << endl << "3) �����" << endl;
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