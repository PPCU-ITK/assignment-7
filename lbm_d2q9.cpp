#include <iostream>
#include <fstream>
#include <cstdlib>

static inline int mod(int v, int m) {
	int val = v%m;
	if (val<0) val = m+val;
	return val;
}

int main(int argc, char ** argv) {
	const int NX = 128;
	const int NY = 128;
	const double OMEGA = 1.0;
	const double rho0 = 1.0;
	const double deltaUX=10e-6;

	const double W[] = {4.0/9.0,1.0/9.0,1.0/36.0,1.0/9.0,1.0/36.0,1.0/9.0,1.0/36.0,1.0/9.0,1.0/36.0};
	const int cx[] = {0,0,1,1, 1, 0,-1,-1,-1};
	const int cy[] = {0,1,1,0,-1,-1,-1, 0, 1};

	const int opposite[] = {0,5,6,7,8,1,2,3,4};

	//Generate random obstacles
	srand(0);
	int * __restrict__ SOLID = new int[NX*NY];
	for (int j = 0; j < NY; j++) {
		for (int i = 0; i < NX; i++) {
			if (rand()%10 >= 7)
				SOLID[j*NX+i] = 1;
			else
				SOLID[j*NX+i] = 0;
		}
	}

	//Initial values
	double * __restrict__ N = new double[NX*NY*9];
	for (int j = 0; j < NY; j++) {
		for (int i = 0; i < NX; i++) {
			for (int f = 0; f < 9; f++) {
				N[(j*NX+i)*9 + f] = rho0 * W[f];
			}
		}
	}

	//Work arrays
	double * __restrict__ workArray = new double[NX*NY*9];
	double * __restrict__ N_SOLID = new double[NX*NY*9];
	double * __restrict__ rho = new double[NX*NY];
	double * __restrict__ ux = new double[NX*NY];
	double * __restrict__ uy = new double[NX*NY];

	//Main time loop
	for (int t = 0; t < 4000; t++) {

		//Backup values
		for (int j = 0; j < NY; j++) {
			for (int i = 0; i < NX; i++) {
				for (int f = 0; f < 9; f++) {
					workArray[(j*NX+i)*9 + f] = N[(j*NX+i)*9 + f];
				}
			}
		}

		//Gather neighbour values
		for (int j = 0; j < NY; j++) {
			for (int i = 0; i < NX; i++) {
				for (int f = 1; f < 9; f++) {
					N[(j*NX+i)*9 + f] = workArray[(mod(j-cy[f],NY)*NX+mod(i-cx[f],NX))*9 + f];
				}
			}
		}

		//Bounce back from solids, no collision
		for (int j = 0; j < NY; j++) {
			for (int i = 0; i < NX; i++) {
				if (SOLID[j*NX+i]==1) {
					for (int f = 0; f < 9; f++) {
						N_SOLID[(j*NX+i)*9 + opposite[f]] = N[(j*NX+i)*9 + f];
					}
				}
			}
		}

		for (int j = 0; j < NY; j++) {
			for (int i = 0; i < NX; i++) {
				rho[j*NX+i] = 0;
				for (int f = 0; f < 9; f++) {
					rho[j*NX+i] += N[(j*NX+i)*9 + f];
				}
			}
		}

		for (int j = 0; j < NY; j++) {
			for (int i = 0; i < NX; i++) {
				ux[j*NX+i] = 0;
				for (int f = 0; f < 9; f++) {
					ux[j*NX+i] += N[(j*NX+i)*9 + f] * cx[f];
				}
				ux[j*NX+i] = ux[j*NX+i] / rho[j*NX+i] + deltaUX;
			}
		}

		for (int j = 0; j < NY; j++) {
			for (int i = 0; i < NX; i++) {
				uy[j*NX+i] = 0;
				for (int f = 0; f < 9; f++) {
					uy[j*NX+i] += N[(j*NX+i)*9 + f] * cy[f];
				}
				uy[j*NX+i] = uy[j*NX+i] / rho[j*NX+i];
			}
		}

		for (int j = 0; j < NY; j++) {
			for (int i = 0; i < NX; i++) {
				for (int f = 0; f < 9; f++) {
					workArray[(j*NX+i)*9 + f] = ux[j*NX+i]*cx[f] + uy[j*NX+i]*cy[f];
				}
			}
		}

		for (int j = 0; j < NY; j++) {
			for (int i = 0; i < NX; i++) {
				for (int f = 0; f < 9; f++) {
					workArray[(j*NX+i)*9 + f] = (3+4.5*workArray[(j*NX+i)*9 + f])*workArray[(j*NX+i)*9 + f];
				}
			}
		}

		for (int j = 0; j < NY; j++) {
			for (int i = 0; i < NX; i++) {
				for (int f = 0; f < 9; f++) {
					workArray[(j*NX+i)*9 + f] = workArray[(j*NX+i)*9 + f] - 1.5 * (ux[j*NX+i]*ux[j*NX+i] + uy[j*NX+i]*uy[j*NX+i]);
				}
			}
		}

		for (int j = 0; j < NY; j++) {
			for (int i = 0; i < NX; i++) {
				for (int f = 0; f < 9; f++) {
					workArray[(j*NX+i)*9 + f] = (1+workArray[(j*NX+i)*9 + f]) * W[f] * rho[j*NX+i];
				}
			}
		}

		for (int j = 0; j < NY; j++) {
			for (int i = 0; i < NX; i++) {
				for (int f = 0; f < 9; f++) {
					N[(j*NX+i)*9 + f] += (workArray[(j*NX+i)*9 + f] - N[(j*NX+i)*9 + f]) * OMEGA;
				}
			}
		}

		for (int j = 0; j < NY; j++) {
			for (int i = 0; i < NX; i++) {
				if (SOLID[j*NX+i]==1) {
					for (int f = 0; f < 9; f++) {
						N[(j*NX+i)*9 + f] = N_SOLID[(j*NX+i)*9 + f];
					}
				}
			}
		}

		//Calculate kinetic energy
		double energy = 0;
		for (int j = 0; j < NY; j++) {
			for (int i = 0; i < NX; i++) {
				energy += ux[j*NX+i]*ux[j*NX+i]+uy[j*NX+i]*uy[j*NX+i];
			}
		}
		if (t%100==0) std::cout << energy << std::endl;
	}

	if (true) {
		std::ofstream myfile;
		myfile.open ("output_velocity.txt");
		for (int j = 0; j < NY; j++) {
			for (int i = 0; i < NX; i++) {
				myfile << SOLID[j*NX+i] << " " << ux[j*NX+i] << " " << uy[j*NX+i] << std::endl;
			}
		}
		myfile.close();
	}
 }
