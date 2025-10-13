#pragma once 

// System right hand side
std::vector<Matrix<double, Dynamic,Dynamic> > F;

// Solution
Matrix<double, Dynamic,Dynamic> u;

// Matrix to map between different meshes
SparseMatrix<double> mesh_matrix;

// Name of write files
const std::string vtk_animation_path = "figures/motivation-transformed-sphere-";

// Include additional files
#include "write.h"
#include "helpful.h"


void simulate_and_save_noise(int Nt, int Nsim, int NN_reference)
{
	// Simulates (W, \varphi_i)_i for i = 1, ..., n, Nt x Nsim times

	// Create reference mesh
	create_plane_mesh(NN_reference);
	reorder_vertices();
	assemble_system();

	std::cout << "Created mesh and assembled system matrices." << std::endl;

	// Sample new noise

	int N = interior.size();

	// Compute square root of mass matrix
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int> > cholesky(massm);

	M_sqrt.resize(N,N);
	M_sqrt = cholesky.matrixL();

	std::cout << "Computed Cholesky factor of mass matrix." << std::endl;

	// Insert values into F

	F = {};

	Matrix<double, Dynamic,Dynamic> F_;
	F_.resize(N,Nt);
	F_.setZero();

	for (int k = 0; k < Nsim; k++)
	{
		for (int i = 0; i < Nt; i++)
		{

			// Print to see progress. 
			if (i % 2000==0)
				std::cout << i/(double)Nt*100. << std::endl; 

			for (int j = 0; j < N; j++)
			{
				// Add the new term
				F_(j,i) = distribution(generator) / std::sqrt(Nt);
				if (i > 0)
					F_(j,i) += F_(j,i-1);
			}
		}

		// Compute (W, \varphi_i) for i = 1, ..., N
		F_ = M_sqrt * F_;

		// Save
		F.push_back(F_);

		// Set F_ to zero
		F_.setZero();
	}

	std::cout << "Simulated noise." << std::endl;

}


void simulate(int NN, int Nsim, double gamma, double tt = 1.)
{
	// Compute h
	double h = std::sqrt(2.)*std::pow(0.5,NN);

	// Compute k
	double Qk = 0.5;

	// Number of time steps is fixed
	int Nt = F[0].cols();

	// Compute time step size 
	double delta_t = tt / Nt;

	// Quadrature constants
	int Kplus = std::ceil(PI*PI/(2.*gamma*Qk*Qk)),
		Kminus = 0; 
	if (gamma < 1.-1.e-9)
		Kminus = -std::ceil(PI*PI/(2.*(1.-gamma)*Qk*Qk));


	// Dimension of V_h
	int N = interior.size();

	// Quadrature parameters
	std::cout << "h: " << h << '\n' << "delta_t: " << delta_t << '\n' <<  "Qk: " << Qk << '\n' << "Kminus: " << Kminus << '\n' << "Kplus: "<< Kplus << std::endl;

	u.resize(vertex_value.size(), Nsim);
	u.setZero(); 

	// Compute Backward Euler system matrix

	// Compute cholesky decomposition of A
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > A_fac(massm + delta_t*vvtm);

	// to hold value of fractional term
	Matrix<double,Dynamic,Dynamic> delta_Q_W(interior.size(),1);
	delta_Q_W.setZero();

	// gamma \in (0,1)
	if (gamma < 1.-1.e-9)
	{
		for (int k = 0; k < Nsim; k++)
		{
		
			auto F_ = F[k];

			// Iterate to time t with timestep size delta_t
			for (int i = 1; i <= Nt; i++)
			{
				// Print to see progress. 
				if (i % 2000==0)
					std::cout << i/(double)Nt*100. << std::endl; 

				// Compute increment of Wiener process
				
				int delta_N = F_.cols() / Nt;
				assert(delta_N == 1);
				
				auto delta_W = F_.col(i*delta_N - 1); // F[k] does not work...
				if (i > 1)
					delta_W -= F_.col((i-1)*delta_N - 1);
				
				// Perform quadrature for the fractional part
				delta_Q_W.setZero();
				#pragma omp parallel for
				for (int ell = Kminus; ell <= Kplus; ell++)
				{
					double yell = ell*Qk;

					// Compute Cholesky decomposion for each term in quadrature
					Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > B_fac(massm + std::exp(-yell)*(massm + stiffnm)); 

					delta_Q_W += std::exp((-gamma)*yell) * B_fac.solve( delta_W );
				}
				u.block(0,k,interior.size(),1) = A_fac.solve( massm * (u.block(0,k,interior.size(),1) + Qk * (std::sin(PI * gamma) / PI) * delta_Q_W) );

			}

		}

	} else { // gamma = 1

		for (int k = 0; k < Nsim; k++)
		{
			auto F_ = F[k];

			// Iterate to time t with timestep size delta_t
			for (int i = 1; i <= Nt; i++)
			{
				// Print to see progress. 
				if (i % 2000==0)
					std::cout << i/(double)Nt*100. << std::endl; 

				// Compute increment of Wiener process
				int delta_N = F_.cols() / Nt;
				assert(delta_N == 1);
				
				auto delta_W = F_.col(i*delta_N - 1); // F[k] does not work...
				if (i > 1)
					delta_W -= F_.col((i-1)*delta_N - 1);
				
				Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > B_fac(massm + stiffnm);
				delta_Q_W = B_fac.solve(delta_W);

				u.block(0,k,interior.size(),1) = A_fac.solve( massm * (u.block(0,k,interior.size(),1) + delta_Q_W)); 

				delta_Q_W.setZero();
			}

		}

	}

}

void simulate_example(int Nt, int Nsim, int NN_reference,
											double gamma)
{
	// Compute reference solution and noise, in addition
	// to creating the reference mesh. 

	// Set seed
	srand(1);

	// Nt \times Nsim 
	simulate_and_save_noise(Nt, Nsim, NN_reference); 

	// Save reference solution
	simulate(NN_reference, Nsim, gamma);
	for (int i = 0; i < u.rows(); i++)
		u(i,0) = -u(i,0);

	// Write reference solution to see
	std::string path_vtk = vtk_animation_path+std::to_string((int)(gamma*100.+1.e-1))+"-"+std::to_string(NN_reference)+".vtk";
	write_vtk_file(path_vtk);

	std::cout << "Computed and saved reference solution" << std::endl;

	return;

}






