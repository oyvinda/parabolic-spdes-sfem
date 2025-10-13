#define transformed_icosasphere_mesh
//#define dome_mesh
//#define square_mesh
//#define noflow_boundary

// FEM functionality
#include "fem.h"

// Functions to simulate u(1) and save it as a .vtk file
#include "simulate-2d-example.h"

// Use: g++ -O3 -fopenmp -march=native main.cpp
// and: compute ./a.out
// to run.
// Use: export OMP_NUM_THREADS=28 to set the number of cores used



int main()
{
	int Nt = 10,
		Nsim = 1;
	double gamma = 0.9;

	int NN_reference = 6;

	simulate_example((1<<Nt), Nsim, NN_reference, gamma);

#ifdef rates
	int Nt = (1 << 12), 
		Nsim = 1;
	double gamma = 1.;

	int NN_reference = 6;

	
	std::vector<int> NN_approxs = {}; 
	for (int i = 1; i < NN_reference; i++)
		NN_approxs.push_back(i);


	// Create matrices that will be reused. 
	for (int nn : NN_approxs)
	{
		create_mesh_matrix(NN_reference, nn, false); 
	}

	std::cout << "Nt = " << Nt << std::endl << "gamma = " << gamma << std::endl;




	std::vector<double> vls = test_example_convergence(Nt, Nsim, NN_reference, NN_approxs, gamma),
						eoc;

	std::cout << "NN_reference: " << NN_reference << '\n' << "gamma: " << gamma << '\n' << "Nsim: " << Nsim <<  std::endl;

	for (int i = 0; i < NN_approxs.size(); i++)
	{
		double x = vls[i];
		int n = NN_approxs[i];
		std::cout << x << ',' << std::endl;

		if (i > 0)
			eoc.push_back(std::log(vls[i-1]/vls[i])/std::log(2.));
	}
	std::cout << "Theoretical convergence rate: " << 2. * gamma + 1 - 0.5 << std::endl;
	for (double x : eoc)
		std::cout << "EOC: " << x << std::endl;
#endif

	return 0;
}




