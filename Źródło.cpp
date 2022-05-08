#include <iostream>
#include <vector>
#include <iomanip>

const double kt = 25;  //30 conductivity - przewodnosc cieplna

const double alpha = 300;	//wspolczynnik wymiany ciepla

const double ambient_temperature = 1200;	//temperatura otoczenia

const double c = 700;	//specific heat - cieplo wlasciwe

const double ro = 7800;	//gestosc

const double dt = 50;	//krok czasowy symulacji

const double t = 500;	//czas symulacji

using namespace std;

int rows = 4, cols = 4;

struct node {
	//kazdy wezel jest okreslony dwoma wspolrzednymi
	float x, y;
	int BC;	//boundary condition - okresla czy wezel znajduje sie na powierzchni
	double initial_t;	//temperatura poczatkowa
};

struct element {
	//kazdy wezel ma swoje id
	//a kazdy element jest reprezentowany przez cztery wezly
	int id[4];

	//pochodne czterech funkcji ksztaltu w czterech punktach calkowania
	double dn_dx[4][4];
	double dn_dy[4][4];

	//hesjan elementu - to suma hesjanow czastkowych,
	//czyli hesjanow obliczonych w kazdym punkcie calkowania
	vector<vector<double>> H = vector<vector<double>>(rows,vector<double>(cols));

	//macierz Hbc
	vector<vector<double>> Hbc = vector<vector<double>>(4, vector<double>(4));

	//wektor P
	vector<double> P = vector<double>(4);

	//macierz C
	vector<vector<double>> C = vector<vector<double>>(4, vector<double>(4));
};

struct GRID {
	float h;	//wysokosc siatki
	float b;	//szerokosc siatki
	int nh;		//liczba wezlow w pionie
	int nb;		//liczba wezlow w poziomie
	int nn;		//liczba wszystkich wezlow
	int ne;		//liczba wszystkich elementow
	vector<node> nodes;		//tablica wezlow siatki
	vector<element> elements;		//tablica elementow siatki
	vector<vector<double>> globalH;
	vector<double> globalP;
	vector<vector<double>> globalC;
	vector<vector<double>> solutionH;
	vector<double> solutionP;
};

struct points{
	double nodes2[3] = { -sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0) };
	double ratios2[2] = { 5.0 / 9, 8.0 / 9 };
	double nodes3[4] = { -0.861136,-0.339981,0.339981,0.861136 };
	double ratios3[4] = { 0.347855,0.652145, 0.652145, 0.347855 };
};

struct Element4_2D {
	//te dwie nazwy moga byc troche mylne,
	//bo eta_nodes przechowuje wartosci ksi
	//i na odwrot, ALE eta_nodes dotyczy pochodnych po eta
	//a ksi_nodes pochodnych po ksi
	//UWAGA KOLEJNOSC PUNKTOW - DO ZMIANY
	double eta_nodes[4] = { -1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3), -1 / sqrt(3) };
	double ksi_nodes[4] = { -1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3) };

	double eta[4][4];
	double ksi[4][4];

	double N_for_C[4][4];

	double ksi_nodes_fixed[4] = { -1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3), -1 / sqrt(3) };
	double eta_nodes_fixed[4] = { -1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3) };
	

	//wspolrzedne punktow calkowania po powierzchni
	//[bok][punkt]
	double pc_powierzchnia_eta[4][2] = { 
		{-1,-1},{-1 / sqrt(3),1 / sqrt(3)}, {1, 1},{1 / sqrt(3), -1 / sqrt(3)} };
	double pc_powierzchnia_ksi[4][2] = {
		{-1 / sqrt(3),1 / sqrt(3)}, {1, 1}, {1 / sqrt(3), -1 / sqrt(3)}, {-1, -1} };

	//funkcje ksztaltu
	//[bok][punkt][funkcja]
	//dla calkowania po powierzchni
	double N[4][2][4];

};

double f(double x) {
	return 5 * x * x + 3 * x + 6;
}

double f(double x, double y) {
	return 5 * x * x * y * y + 3 * x * y + 6;
}

double gauss1d(points p, int num_of_pts) {
	double solution = 0;

	if (num_of_pts == 2) {
		for (int i = 0; i < 3; i++) {
			solution += p.ratios2[i%2] * f(0 + p.nodes2[i]);
		}
	}
	else if (num_of_pts == 3) {
		for (int i = 0; i < 4; i++) {
			solution += p.ratios3[i] * f(0 + p.nodes3[i]);
		}
	}
	else
		cout << "Invalid number of points.";


	return solution;
}

double gauss2d(points p, int num_of_pts) {
	double solution = 0;

	if (num_of_pts == 2) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				solution += p.ratios2[i % 2] * p.ratios2[j % 2] * f(0 + p.nodes2[i], 0 + p.nodes2[j]);
			}
		}
	}
	else if (num_of_pts == 3) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				solution += p.ratios3[i] * p.ratios3[j] * f(0 + p.nodes3[i], 0 + p.nodes3[j]);
			}
		}
	}
	else
		cout << "Invalid number of points.";

	return solution;
}

void fillTheGrid(GRID &g) {

	int index;

	//making nodes
	for (int i = 0; i < g.nb; i++) {
		for (int j = 0; j < g.nh; j++) {
			index = i * g.nh + j;
			g.nodes.push_back(node());
			g.nodes[index].x = i * g.b / (g.nb - 1);
			g.nodes[index].y = j * g.h / (g.nh - 1);
			g.nodes[index].initial_t = 100;
			if (g.nodes[index].y == 0 || g.nodes[index].y == g.h ||
				g.nodes[index].x == 0 || g.nodes[index].x == g.b) {
				g.nodes[index].BC = 1;
			}
			//else
				//g.nodes[index].BC = 0;
		}
	}

	//making elements
	int k = -1;

	for (int i = 0; i < g.ne; i++) {
		g.elements.push_back(element());

		if (i % (g.nh - 1) == 0)
			k++;

		g.elements[i].id[0] = i + 1 + k;
		g.elements[i].id[1] = i + g.nh + 1 + k;
		g.elements[i].id[2] = i + g.nh + 2 + k;
		g.elements[i].id[3] = i + 2 + k;
		
	}
}

//wypisuje siatke:
//wspolrzedne kolejnych wezlow
//oraz elementy w postaci id wezlow z ktorych sie skladaja
void printOutGrid(GRID g) {
	cout << "nodes:";
	for (int i = 0; i < g.nn; i++) {
		cout << "\nnode " << i + 1 << ": " << "\n\tx = " << g.nodes[i].x << "\n\ty = " << g.nodes[i].y;
		cout << "\n\tBC = " << g.nodes[i].BC;
	}
	cout << "\n\nelements:";
	for (int i = 0; i < g.ne; i++) {
		cout << "\nelement " << i + 1 << ": " << g.elements[i].id[0] << " ";
		cout << g.elements[i].id[1] << " " << g.elements[i].id[2] << " " << g.elements[i].id[3];
	}
}

//wypelnia tablice 4-wezlowego elementu 2d
//pochodnymi kolejnych funkcji ksztaltu w kazdym z czterech punktow calkowania 
void fillArrays(Element4_2D& e) {

	//[punkt calkowania][funkcja]
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			switch (j) {
			case 0:		//N1
				//pochodne funkcji ksztaltu w elemencie
				e.eta[i][0] = -(1 - e.eta_nodes[i]) / 4;
				e.ksi[i][0] = -(1 - e.ksi_nodes[i]) / 4;
				//funkcje ksztaltu na powierzchni
				e.N[i][0][0] = (1 - e.pc_powierzchnia_ksi[i][0])*(1 - e.pc_powierzchnia_eta[i][0])/4;
				e.N[i][1][0] = (1 - e.pc_powierzchnia_ksi[i][1])*(1 - e.pc_powierzchnia_eta[i][1])/4;

				e.N_for_C[i][j]= (1 - e.ksi_nodes[(4 - i) % 4]) * (1 - e.eta_nodes[(4 - i) % 4]) / 4;
				break;
			case 1:		//N2
				e.eta[i][1] = -(1 + e.eta_nodes[i]) / 4;
				e.ksi[i][1] = (1 - e.ksi_nodes[i]) / 4;

				e.N[i][0][1] = (1 + e.pc_powierzchnia_ksi[i][0]) * (1 - e.pc_powierzchnia_eta[i][0]) / 4;
				e.N[i][1][1] = (1 + e.pc_powierzchnia_ksi[i][1]) * (1 - e.pc_powierzchnia_eta[i][1]) / 4;

				e.N_for_C[i][j] = (1 + e.ksi_nodes[(4 - i) % 4]) * (1 - e.eta_nodes[(4 - i) % 4]) / 4;
				break;
			case 2:		//N3
				e.eta[i][2] = (1 + e.eta_nodes[i]) / 4;
				e.ksi[i][2] = (1 + e.ksi_nodes[i]) / 4;

				e.N[i][0][2] = (1 + e.pc_powierzchnia_ksi[i][0]) * (1 + e.pc_powierzchnia_eta[i][0]) / 4;
				e.N[i][1][2] = (1 + e.pc_powierzchnia_ksi[i][1]) * (1 + e.pc_powierzchnia_eta[i][1]) / 4;

				e.N_for_C[i][j] = (1 + e.ksi_nodes[(4 - i) % 4]) * (1 + e.eta_nodes[(4 - i) % 4]) / 4;
				break;
			case 3:		//N4
				e.eta[i][3] = (1 - e.eta_nodes[i]) / 4;
				e.ksi[i][3] = -(1 + e.ksi_nodes[i]) / 4;

				e.N[i][0][3] = (1 - e.pc_powierzchnia_ksi[i][0]) * (1 + e.pc_powierzchnia_eta[i][0]) / 4;
				e.N[i][1][3] = (1 - e.pc_powierzchnia_ksi[i][1]) * (1 + e.pc_powierzchnia_eta[i][1]) / 4;

				e.N_for_C[i][j] = (1 - e.ksi_nodes[(4 - i) % 4]) * (1 + e.eta_nodes[(4 - i) % 4]) / 4;
				break;
			}
		}
	}
}

//wyznacznik macierzy drugiego stopnia
double det2(vector < vector < double >> j) {
	return j[0][0] * j[1][1] - j[0][1] * j[1][0];
}

//jakobian przeksztalcenia
//jakobian odwrotny
void jakobian_odwrotny(vector<vector<double>> j, vector<vector<double>>& j_inv) {
	double detJ = det2(j);

	//jakobian odwrotny
	j_inv[0][0] = 1/detJ * j[1][1];
	j_inv[0][1] = -1/detJ * j[0][1];
	j_inv[1][0] = -1/detJ * j[1][0];
	j_inv[1][1] = 1/detJ * j[0][0];
}

//en-numer elementu
//pcn-numer punktu calkowania
void jakobian(int en, int pcn, vector<vector <double>> &j, vector<vector<double>> &j_inv, GRID g, Element4_2D e) {

	//interpolacja pochodnych x,y po ksi,eta
	for (int p = 0; p < 4; p++) {
		//jakobian
		// dx/dksi
		j[0][0] += e.ksi[pcn][p] * g.nodes[g.elements[en].id[p]-1].x;
		// dy/dksi
		j[0][1] += e.ksi[pcn][p] * g.nodes[g.elements[en].id[p]-1].y;
		// dx/deta
		j[1][0] += e.eta[pcn][p] * g.nodes[g.elements[en].id[p]-1].x;
		// dy/deta
		j[1][1] += e.eta[pcn][p] * g.nodes[g.elements[en].id[p]-1].y;
	}

	jakobian_odwrotny(j, j_inv);
}

void printJakobian(vector<vector<double>> j) {
	
	for (int i = 0; i < 2; i++) {
		for (int k = 0; k < 2; k++) {
			cout << j[i][k] << "\t";
		}
		cout << endl;
	}
}

void printDerivatives(Element4_2D e) {
	cout << "\npochodne dn/dksi:\n";
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << e.ksi[i][j] << "\t";
		}
		cout << "\n";
	}

	cout << "\npochodne dn/deta:\n";
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << e.eta[i][j] << "\t";
		}
		cout << "\n";
	}
}

void twoMatricesSum(vector<vector<double>> &m, vector<vector<double>> m1, vector<vector<double>> m2, int r, int c) {
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			m[i][j] = m1[i][j] + m2[i][j];
		}
	}
}

//h-macierz wyjsciowa, dn_dxy - pochodne przechowywane w elemencie 
void solveH(vector<vector<double>> &h, double* dn_dx, double* dn_dy, int r, int c, double kt, double dv) {
	vector<vector<double>> dn_dx_pom(r, vector<double>(c));
	vector<vector<double>> dn_dy_pom(r, vector<double>(c));
	//vector<vector<double>> h_pom;

	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			dn_dx_pom[i][j] = dn_dx[j] * dn_dx[i];
			dn_dy_pom[i][j] = dn_dy[j] * dn_dy[i];
		}
	}

	twoMatricesSum(h, dn_dx_pom, dn_dy_pom, 4, 4);
	
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			h[i][j] *= kt * dv;
		}
	}
}

void dn_dxy(GRID &g, Element4_2D e, vector<vector<double>> &jk, vector<vector<double>> &j_inv) {

	vector<vector<double>> Hpom(4,vector<double>(4));

	//dla kazdego elementu siatki
	for (int i = 0; i < g.ne; i++) {
		//dla kazdego punktu calkowania
		for (int j = 0; j < 4; j++) {
			//wyzerowanie jakobianu
			for (int x = 0; x < 2; x++) {
				for (int y = 0; y < 2; y++) {
					jk[x][y] = 0;
					j_inv[x][y] = 0;
				}
			}
			jakobian(i, j, jk, j_inv, g, e);
			//pochodne czterech funkcji ksztaltu
			for (int k = 0; k < 4; k++) {
				//jakobian(i, j, jk, j_inv, g, e);
				g.elements[i].dn_dx[j][k] = j_inv[0][0] * e.ksi[j][k] + j_inv[0][1] * e.eta[j][k];
				g.elements[i].dn_dy[j][k] = j_inv[1][0] * e.ksi[j][k] + j_inv[1][1] * e.eta[j][k];
			}

			solveH(Hpom,g.elements[i].dn_dx[j],g.elements[i].dn_dy[j],4,4,kt,det2(jk));
			//H = Hstare + Hnowe (H1+H2+H3+H4)
			twoMatricesSum(g.elements[i].H,g.elements[i].H,Hpom,4,4);
		}
	}
}

void solveHbc(GRID& g, Element4_2D e) {

	vector<vector<double>> Hbc_pom = vector<vector<double>>(4, vector<double>(4));
	vector<vector<double>> Npc1_pom = vector<vector<double>>(4, vector<double>(4));
	vector<vector<double>> Npc2_pom = vector<vector<double>>(4, vector<double>(4));
	double detJ;
	double L;
	double w1 = 1, w2 = 1;

	//dla kazdego elementu siatki
	for (int i = 0; i < g.ne; i++) {
		//dla kazdego wezla/boku
		for (int j = 0; j < 4; j++) {
			//sprawdzany jest warunek tworzenia przez dwa kolejne wezly
			//sciany bedacej na powierzchni
			if (g.nodes[g.elements[i].id[j] - 1].BC == 1 &&
				g.nodes[g.elements[i].id[(j + 1) % 4] - 1].BC == 1) {
				//jesli spelniony to liczona jest macierz Hbc
				//dla jednej sciany
				//hbc = (w1 * alpha * ({n} * {n}) + w2 * alpha * ({n} * {n})) *  L/2
				L = sqrt(pow(g.nodes[g.elements[i].id[j] - 1].x -
					g.nodes[g.elements[i].id[(j + 1) % 4] - 1].x, 2) +
					pow(g.nodes[g.elements[i].id[j] - 1].y -
						g.nodes[g.elements[i].id[(j + 1) % 4] - 1].y, 2));
				detJ = L / 2;
				//calkowanie (powierzchnia)
				for (int k = 0; k < 4; k++) {
					for (int l = 0; l < 4; l++) {
						Npc1_pom[k][l] = e.N[j][0][k] * e.N[j][0][l] * alpha * w1 * detJ;
						Npc2_pom[k][l] = e.N[j][1][k] * e.N[j][1][l] * alpha * w2 * detJ;
					}
				}
				//Hbc = Hbc_pc1 + Hbc_pc2
				twoMatricesSum(Hbc_pom, Npc1_pom, Npc2_pom, 4, 4);
				twoMatricesSum(g.elements[i].Hbc, g.elements[i].Hbc, Hbc_pom,4,4);

			}
		}
	}
}

void solveP(GRID& g, Element4_2D e) {

	vector<double> p1_pom = vector<double>(4);
	vector<double> p2_pom = vector<double>(4);

	double w1 = 1, w2 = 1;

	double detJ;
	double L;

	//dla kazdego elementu siatki
	for (int i = 0; i < g.ne; i++) {
		//dla kazdego wezla/boku
		for (int j = 0; j < 4; j++) {
			//sprawdzany jest warunek tworzenia przez dwa kolejne wezly
			//sciany bedacej na powierzchni
			if (g.nodes[g.elements[i].id[j] - 1].BC == 1 &&
				g.nodes[g.elements[i].id[(j + 1) % 4] - 1].BC == 1) {

				L = sqrt(pow(g.nodes[g.elements[i].id[j] - 1].x -
					g.nodes[g.elements[i].id[(j + 1) % 4] - 1].x, 2) +
					pow(g.nodes[g.elements[i].id[j] - 1].y -
						g.nodes[g.elements[i].id[(j + 1) % 4] - 1].y, 2));

				detJ = L / 2;

				for (int k = 0; k < 4; k++) {
					p1_pom[k] = e.N[j][0][k] * w1 * ambient_temperature * alpha * detJ;
					p2_pom[k] = e.N[j][1][k] * w2 * ambient_temperature * alpha *detJ;
					g.elements[i].P[k] += p1_pom[k] + p2_pom[k];
				}
			}
		}
	}
}

void solveGlobalH(GRID& g) {
	vector<vector<double>> globalH = vector<vector<double>>(g.nn, vector<double>(g.nn));
	for (int i = 0; i < g.ne; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				g.globalH[g.elements[i].id[j] - 1][g.elements[i].id[k] - 1] +=
					g.elements[i].H[j][k] + g.elements[i].Hbc[j][k];
			}
		}
	}
}

void solveGlobalP(GRID& g) {
	
	for (int i = 0; i < g.ne; i++) {
		for (int j = 0; j < 4; j++) {
			g.globalP[g.elements[i].id[j] - 1] += g.elements[i].P[j];
		}
	}
}

void solutionH(GRID &g, Element4_2D e) {

	for (int i = 0; i < g.nn; i++) {
		for (int j = 0; j < g.nn; j++) {
			g.solutionH[i][j] = g.globalH[i][j] + g.globalC[i][j] / dt;
		}
	}
}

void solutionP(GRID& g, Element4_2D e) {
	vector<double> pom1(g.nn);
	double pom2=0;
	for (int i = 0; i < g.nn; i++) {
		for (int j = 0; j < g.nn; j++) {
			pom2 += g.globalC[i][j] / dt * g.nodes[j].initial_t;
		}
		//pom1[i] = pom2;
		g.solutionP[i] = g.globalP[i] + pom2;
		pom2 = 0;
	}


}

void solveC(GRID& g, Element4_2D e) {

	vector<vector<double>> jk = vector<vector<double>>(2,vector<double>(2));
	vector<vector<double>> j_inv = vector<vector<double>>(2, vector<double>(2));
	vector<vector<double>> C_pom = vector<vector<double>>(4, vector<double>(4));

	double detJ;
	double w1 = 1, w2 = 1;

	//dla kazdego elementu siatki
	for (int i = 0; i < g.ne; i++) {

		//dla kazdego punktu calkowania
		for (int j = 0; j < 4; j++) {
			//wyzerowanie jakobianu
			for (int x = 0; x < 2; x++) {
				for (int y = 0; y < 2; y++) {
					jk[x][y] = 0;
					j_inv[x][y] = 0;
				}
			}

			//policzenie jakobianu
			jakobian(i, j, jk, j_inv, g, e);
			detJ = det2(jk);

			//wymnozenie 
			for (int k = 0; k < 4; k++) {
				for (int m = 0; m < 4; m++) {
					//tu nastapila zmiana w drugim N_for_C[j][k]
					C_pom[k][m] = ro * c * w1 * w2 * e.N_for_C[j][k] * e.N_for_C[j][m] * detJ;
				}
			}

			twoMatricesSum(g.elements[i].C, g.elements[i].C, C_pom, 4, 4);
		}
	}
}

void solveGlobalC(GRID& g, Element4_2D e) {

	for (int i = 0; i < g.ne; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				g.globalC[g.elements[i].id[j] - 1][g.elements[i].id[k] - 1] += g.elements[i].C[j][k];
			}
		}
	}
}

void upperTriangularMatrix(vector<vector<double>> &matrix, vector<double>& constantTerms) {
	unsigned int n = matrix.size();
	double ratio;

	for (int i = 0; i < n - 1; i++) {
		for (int j = i + 1; j < n; j++) {
			ratio = matrix[j][i] / matrix[i][i];
			for (int k = 0; k < n; k++) {
				matrix[j][k] -= matrix[i][k] * ratio;
			}
			constantTerms[j] -= constantTerms[i] * ratio;
		}
	}
}

vector<double> gaussEliminationSolve(vector<vector<double>>& matrix, vector<double>& constantTerms) {
	upperTriangularMatrix(matrix, constantTerms);

	double difference;
	unsigned int n = matrix.size();
	vector<double> solution(n);

	for (int i = n - 1; i >= 0; i--) {
		difference = 0;
		if (i == n - 1) {
			solution[i] = constantTerms[i] / matrix[i][i];
		}
		else {
			for (int j = i + 1; j < n; j++) {
				difference += matrix[i][j] * solution[j];
			}
			solution[i] = (constantTerms[i] - difference) / matrix[i][i];
		}
	}

	return solution;
}

double findMin(vector<double> v) {
	double min = v[0];
	for (int i = 1; i < v.size(); i++) {
		if (v[i] < min)
			min = v[i];
	}
	return min;
}

double findMax(vector<double> v) {
	double max = v[0];
	for (int i = 1; i < v.size(); i++) {
		if (v[i] > max)
			max = v[i];
	}
	return max;
}

int main() {

	GRID g;
	g.h = 0.1;  //0.1
	g.b = 0.1;	//0.1
	g.nh = 4;	//5 //31
	g.nb = 4;	//5 //31
	g.nn = g.nh * g.nb;
	g.ne = (g.nh - 1) * (g.nb - 1);
	g.globalH = vector<vector<double>>(g.nn, vector<double>(g.nn));
	g.globalC = vector<vector<double>>(g.nn, vector<double>(g.nn));
	g.solutionH = vector<vector<double>>(g.nn, vector<double>(g.nn));
	g.solutionP = vector<double>(g.nn);
	g.globalP = vector<double>(g.nn);

	fillTheGrid(g);
	printOutGrid(g);

	points p;

	/*double result2 = gauss1d(p, 2);
	double result3 = gauss1d(p, 3);

	cout << setprecision(14) << result2 << endl;
	cout << setprecision(14) << result3 << endl;

	result2 = gauss2d(p, 2);
	result3 = gauss2d(p, 3);

	cout << setprecision(14) << result2 << endl;
	cout << setprecision(14) << result3 << endl;*/

	Element4_2D e;

	fillArrays(e);

	//printDerivatives(e);

	vector<vector<double>> j(2, vector<double>(2, 0));
	vector<vector<double>> j_inv(2, vector<double>(2, 0));

	//jakobian(0, 0, j, j_inv, g, e);

	//cout << "\njakobian: " << endl;
	//printJakobian(j);

	//cout << "\n\njakobian odwrotny: " << endl;
	//printJakobian(j_inv);

	cout << "\n\npochodne dn/dx elementu 1: " << endl;

	dn_dxy(g, e, j, j_inv);

	/*for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << g.elements[0].dn_dx[i][j] << "\t";
		}
		cout << endl;
	}*/

	cout << "\n\nmacierz H elementu 1: " << endl;

	/*for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << g.elements[0].H[i][j] << "\t";
		}
		cout << endl;
	}*/

	printOutGrid(g);

	cout << "funkcje ksztaltu w punktach calkowania na powierzchni (element 1): \n";


	solveHbc(g, e);
	cout << "\n\nmacierz Hbc elementu 1: " << endl;
	/*for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << g.elements[0].Hbc[i][j] << "\t";
		}
		cout << endl;
	}*/
	cout << "\n\nmacierz Hbc elementow: " << endl;
	/*for (int k = 0; k < g.ne; k++) {
		cout << "\nelement " << k + 1 << ": \n";
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				cout << g.elements[k].Hbc[i][j] << "\t";
			}
			cout << endl;
		}
	}*/
	cout << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << e.N[i][0][j] << endl;
		}
	}

	solveP(g, e);
	cout << "\nwektor p:\n";
	for (int i = 0; i < 4; i++) {
		cout << g.elements[0].P[i] << endl;
	}
	
	solveGlobalH(g);
	cout << "\nglobalna macierz H (agregacja):\n";
	/*for (int i = 0; i < g.nn; i++) {
		for (int j = 0; j < g.nn; j++) {
			cout << g.globalH[i][j] << " ";
		}
		cout << endl;
	}*/

	solveGlobalP(g);
	/*cout << "\nwektor P (agregacja):\n";
	for (int i = 0; i < g.nn; i++) {
		cout << g.globalP[i] << " ";
	}*/

	solveC(g, e);
	cout << "\nmacierz C elementu 1:" << endl;
	/*for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << g.elements[0].C[i][j] << " ";
		}
		cout << endl;
	}*/

	cout << "\nFunkcje ksztaltu:" << endl;
	/*for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << e.N_for_C[i][j] << " ";
		}
		cout << endl;
	}*/

	solveGlobalC(g, e);
	cout << "\nglobalna macierz C (agregacja):\n";
	/*for (int i = 0; i < g.nn; i++) {
		for (int j = 0; j < g.nn; j++) {
			cout << g.globalC[i][j] << " ";
		}
		cout << endl;
	}*/

	/*solutionH(g, e);
	cout << "\nMacierz H (rownanie H=H+C/dt):\n";
	for (int i = 0; i < g.nn; i++) {
		for (int j = 0; j < g.nn; j++) {
			cout << g.solutionH[i][j] << " ";
		}
		cout << endl;
	}*/

	/*solutionP(g, e);
	cout << "\nWektor P (rownanie P=P+(C/dt)*T0):\n"; 
	for (int i = 0; i < g.nn; i++) {
		cout << g.solutionP[i] << " ";
	}*/

	vector<double> solutionT;// = gaussEliminationSolve(g.solutionH, g.solutionP);
	double min;
	double max;

	for (int i = 0; i < t / dt; i++) {
		if (i == 0) {
			cout << "\nSymulacja\n";
			cout << "Time \tMinTemp \tMaxTemp\n";
		}
		solutionH(g, e);
		solutionP(g, e);
		solutionT = gaussEliminationSolve(g.solutionH, g.solutionP);
		for (int j = 0; j < g.nn; j++) {
			g.nodes[j].initial_t = solutionT[j];
		}
		min = findMin(solutionT);
		max = findMax(solutionT);
		cout << (i + 1) * dt << " \t" << min << " \t" << max << endl;
	}

	return 0;
}