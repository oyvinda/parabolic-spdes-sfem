#define line_mesh

// FEM functionality
#include "fem.h"

// Include functions to compute the L^2(\Omega; H) and pathwise error
#include "compute-error.h"


// Use: g++ -O3 -fopenmp -march=native main.cpp
// and: compute ./a.out
// to run.
// Use: export OMP_NUM_THREADS=28 to set the number of cores used

#define space_convergence

int main()
{

#ifdef time_convergence
	int Nt = 20, // 20
		Nsim = 1;
	double gamma = 0.75;

	int NN_reference = 11; //11;

	std::vector<int> NN_approxs = {};
	for (int i = 1; i < Nt; i++)
		NN_approxs.push_back(1 << i);

	std::cout << "Nt = " << (1 << Nt) << std::endl << "gamma = " << gamma << std::endl;

	std::vector<double> vls = test_time_convergence(NN_reference, (1 << Nt), Nsim, NN_approxs, gamma),
		eoc;
#endif

#ifdef space_convergence
	int Nt = 20,
		Nsim = 1;
	double gamma = 0.75;

	int NN_reference = 11;

	std::vector<int> NN_approxs = {};
	for (int i = 1; i < NN_reference; i++)
		NN_approxs.push_back(i);

	// Create matrices that will be reused.
	for (int nn : NN_approxs)
	{
		create_mesh_matrix(NN_reference, nn, false);
	}

	std::cout << "Nt = " << (1 << Nt) << std::endl << "gamma = " << gamma << std::endl;

	std::vector<double> vls = test_space_convergence(NN_reference, (1 << Nt), Nsim, NN_approxs, gamma),
		eoc;
#endif


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


	return 0;
}




