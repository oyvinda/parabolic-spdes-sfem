#pragma once 

// System right hand side
std::vector<Matrix<double, Dynamic,Dynamic> > F;

// b process for paper 4
std::vector<Matrix<double, Dynamic,1> > b_process;

// Solution
Matrix<double, Dynamic,Dynamic> u;

// Matrix to map between different meshes
SparseMatrix<double> mesh_matrix;

// Name of write files
const std::string vtk_animation_path = "figures/animation-test-04-";

// Include additional files
#include "write.h"
#include "helpful.h"


void simulate_and_save_noise(int NN_reference, int Nt, int Nsim)
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


	// Simulate b process

	// Assemble 1d mass and stiffness matrices
	SparseMatrix<double> massm_, stiffnm_;
	massm_.resize(Nt,Nt);
	stiffnm_.resize(Nt,Nt);
	std::vector<Triplet<double>> coeff_massm, 
								 coeff_stiffnm;

	for (int i = 0; i < Nt; i++)
	{
		Triplet<double> massm_aij = Triplet<double>(i,i, 2./(3.*Nt));
		coeff_massm.push_back(massm_aij);

		Triplet<double> stiffn_aij = Triplet<double>(i,i, 2.*Nt);
		coeff_stiffnm.push_back(stiffn_aij);

		if (i > 0)
		{
			massm_aij = Triplet<double>(i,i-1, 1./(6.*Nt));
			coeff_massm.push_back(massm_aij);

			stiffn_aij = Triplet<double>(i,i-1, -Nt);
			coeff_stiffnm.push_back(stiffn_aij);
		}
		if (i < Nt - 1)
		{
			massm_aij = Triplet<double>(i,i+1, 1./(6.*Nt));
			coeff_massm.push_back(massm_aij);

			stiffn_aij = Triplet<double>(i,i+1, -Nt);
			coeff_stiffnm.push_back(stiffn_aij);
		}
	}
	massm_.setFromTriplets(coeff_massm.begin(), coeff_massm.end());
	stiffnm_.setFromTriplets(coeff_stiffnm.begin(), coeff_stiffnm.end());

	// Cholesky of 1d massm
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int> > cholesky_(massm_);

	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > B_fac(massm_ + stiffnm_);

	for (int k = 0; k < Nsim; k++)
	{
		Matrix<double, Dynamic,1> bt;
		bt.resize(Nt,1);

		for (int i = 0; i < Nt; i++)
			bt(i,0) = distribution(generator);

		bt = B_fac.solve(cholesky_.matrixL() * bt);

		for (int i = 0; i < Nt; i++)
			bt(i,0) = std::exp(1.*bt(i,0)*bt(i,0));

		b_process.push_back(bt);
	}

}

std::vector<Matrix<double,Dynamic,Dynamic> > compute_load_vectors_from_noise(const std::vector<Matrix<double,Dynamic,Dynamic> >& noise)
{
	// Assumes mesh_matrix is already created

	// Return the transformed noise
	std::vector<Matrix<double,Dynamic,Dynamic> > res;
	for (auto& mat : noise)
	{
		res.push_back(mesh_matrix*mat);
	}

	return res;
}


void simulate(int NN, int Nt, int Nsim, double gamma, double tt = 1.)
{
	// Compute h
	double h = std::sqrt(2.)*std::pow(0.5,NN);

	// Compute k
	double Qk = 0.5;

	// Compute time step size 
	double delta_t = tt / Nt;

	// Quadrature constants
	int Kplus = 0,
		Kminus = 0;
	if (gamma > 1.e-5)
		Kplus = std::ceil(PI*PI/(2.*gamma*Qk*Qk)); 
	if (gamma < 1.-1.e-5)
		Kminus = -std::ceil(PI*PI/(2.*(1.-gamma)*Qk*Qk));


	// Dimension of V_h
	int N = interior.size();

	// Quadrature parameters
	std::cout << "h: " << h << '\n' << "delta_t: " << delta_t << '\n' <<  "Qk: " << Qk << '\n' << "Kminus: " << Kminus << '\n' << "Kplus: "<< Kplus << std::endl;

	u.resize(vertex_value.size(), Nsim);
	u.setZero(); 

	// Compute Backward Euler system matrix

	// Compute cholesky decomposition of A
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > A_fac(massm + delta_t*(stiffnm));

	// to hold value of fractional term
	Matrix<double,Dynamic,Dynamic> delta_Q_W(interior.size(),1);
	delta_Q_W.setZero();

	// gamma \in [0,1)
	if (gamma < 1.-1.e-5)
	{
		for (int k = 0; k < Nsim; k++)
		{
		
			Matrix<double,Dynamic,Dynamic> F_ = F[k];

			// Iterate to time t with timestep size delta_t
			for (int i = 1; i <= Nt; i++)
			{
				// Print to see progress. 
				if (i % 2000==0)
					std::cout << i/(double)Nt*100. << std::endl; 

				// Compute increment of Wiener process
				
				int delta_N = F_.cols() / Nt;
				//assert(delta_N == 1);
				
				Matrix<double,Dynamic,1> delta_W = F_.col(i*delta_N - 1); // F[k] does not work...
				if (i > 1)
					delta_W -= F_.col((i-1)*delta_N - 1);
				/*
				// Perform quadrature for the fractional part
				//#pragma omp parallel for
				delta_Q_W.setZero();
				for (int ell = Kminus; ell <= Kplus; ell++)
				{
					double yell = ell*Qk;

					// Compute Cholesky decomposion for each term in quadrature
					Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > B_fac(massm + std::exp(-yell)*(massm + stiffnm)); 

					delta_Q_W += std::exp((-gamma)*yell) * B_fac.solve( delta_W );
				}
				u.block(0,k,interior.size(),1) = A_fac.solve( massm * (u.block(0,k,interior.size(),1) + Qk * (std::sin(PI * gamma) / PI) * delta_Q_W) );
				delta_Q_W.setZero();
				*/
				
				u.block(0,k,interior.size(),1) = A_fac.solve(massm * u.block(0,k,interior.size(),1) + b_process[k]((i-1)*delta_N,0)*delta_W);

			}

			// If gamma = 0
			if (gamma < 1.e-5)
				return;

			// Perform quadrautre for the fractional part when matrices commute
			//#pragma omp parallel for
			delta_Q_W.setZero();
			for (int ell = Kminus; ell <= Kplus; ell++)
			{
				double yell = ell*Qk;

				// Compute Cholesky decomposion for each term in quadrature
				Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > B_fac(massm + std::exp(-yell)*(massm + stiffnm));

				delta_Q_W += std::exp((-gamma)*yell) * B_fac.solve( massm * u.block(0,k,interior.size(),1) ); // Remember massm
			}
			u.block(0,k,interior.size(),1) = Qk * (std::sin(PI * gamma) / PI) * delta_Q_W;

		}

	} else { // gamma = 1

		for (int k = 0; k < Nsim; k++)
		{
			Matrix<double,Dynamic,Dynamic> F_ = F[k];

			// Iterate to time t with timestep size delta_t
			for (int i = 1; i <= Nt; i++)
			{
				// Print to see progress. 
				if (i % 2000==0)
					std::cout << i/(double)Nt*100. << std::endl; 

				// Compute increment of Wiener process
				int delta_N = F_.cols() / Nt;
				//assert(delta_N == 1);
				
				Matrix<double,Dynamic,1> delta_W = F_.col(i*delta_N - 1); // F[k] does not work...
				if (i > 1)
					delta_W -= F_.col((i-1)*delta_N - 1);
				
				//Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > B_fac(massm + stiffnm);
				//delta_Q_W = B_fac.solve(delta_W);

				//u.block(0,k,interior.size(),1) = A_fac.solve( massm * (u.block(0,k,interior.size(),1) + delta_Q_W)); 

				//delta_Q_W.setZero();

				u.block(0,k,interior.size(),1) = A_fac.solve( massm * u.block(0,k,interior.size(),1) + b_process[k]((i-1)*delta_N,0)*delta_W);
			}

			// When matrices commute
			Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > B_fac(massm + stiffnm);
			u.block(0,k,interior.size(),1) = B_fac.solve( massm * u.block(0,k,interior.size(),1) ); // Remember massm

		}

	}

}


double approximate_L2_error(const Matrix<double,Dynamic,Dynamic>& u_ref, const Matrix<double,Dynamic,Dynamic>& u_approx)
{
	// Approximates error on finer mesh.
	// Assumes mesh_matrix has been created,
	// and that massm is of the finer mesh.

	Matrix<double,Dynamic,Dynamic> u_approx_projected_up(u_ref.rows(),u_ref.cols());
	u_approx_projected_up.setZero();

	u_approx_projected_up.block(0,0,mesh_matrix.cols(),u_ref.cols()) = mesh_matrix.transpose() * u_approx.block(0,0,mesh_matrix.rows(),u_ref.cols());

	double res = 0.;
	for (int k = 0; k < u_ref.cols(); k++)
	{
		// Sum up for every sample
		auto d = (u_approx_projected_up-u_ref).block(0,k,mesh_matrix.cols(),1);
		res += (d.transpose() * massm * d)(0,0) / u_ref.cols();
	}
	res = std::sqrt(res);

	return res;
}
double approximate_L2_norm(const Matrix<double,Dynamic,Dynamic>& u_ref)
{
	// Approximates error on finer mesh.
	// Assumes mesh_matrix has been created,
	// and that massm is of the finer mesh.

	double res = 0.;
	for (int k = 0; k < u_ref.cols(); k++)
	{
		// Sum up for every sample
		auto d = u_ref.block(0,k, massm.rows(), 1);
		res += (d.transpose() * massm * d)(0,0) / u_ref.cols();
	}
	res = std::sqrt(res);

	return res;
}

std::vector<double> test_space_convergence(int NN_reference, int Nt, int Nsim, const std::vector<int>& NN_approxs,
											double gamma)
{
	// Compute reference solution and noise, in addition
	// to creating the reference mesh. 

	// Set seed
	srand(1);

	// Nt \times Nsim 
	simulate_and_save_noise(NN_reference, Nt, Nsim); 

	// Save reference solution
	simulate(NN_reference, Nt, Nsim, gamma);

	Matrix<double, Dynamic,Dynamic> u_reference;
	u_reference.resize(vertex_value.size(), Nsim);
	u_reference.setZero();
	u_reference = u;

	// Write reference solution to see
	std::string path_vtk = vtk_animation_path+std::to_string(NN_reference)+".vtk";
	write_vtk_file(path_vtk);

	std::cout << "Computed and saved reference solution" << std::endl;

	// Save mass matrix for later use
	auto massm_reference = massm;

	// Save noise for later use
	auto F_reference = F;

	
	std::vector<double> res;

	for (int NN_approx : NN_approxs)
	{
		// Create approximation mesh
		create_plane_mesh(NN_approx);
		reorder_vertices();
		assemble_system();

		// Compute load vectors for the approximate models
		create_mesh_matrix(NN_reference, NN_approx, true);

		// Compute the noise
		F = compute_load_vectors_from_noise(F_reference);

		// Simulate
		simulate(NN_approx, Nt, Nsim, gamma);

		//Matrix<double, Dynamic, Dynamic> u_approx;
		//u_approx.resize(vertex_value.size(), Nsim);
		//u_approx.setZero();
		//u_approx = u;

		// Write vtk file to see
		std::string path_vtk = vtk_animation_path+std::to_string(NN_approx)+".vtk";
		write_vtk_file(path_vtk);

		// Set mass matrix equal to that of the reference mesh 
		massm.resize(massm_reference.rows(), massm_reference.cols());
		massm.setZero();
		massm = massm_reference;

		// Approximate the (relative) L2 error
		double e = approximate_L2_error(u_reference, u) / approximate_L2_norm(u_reference);
		res.push_back(e); 
	}

	return res;
}

std::vector<double> test_time_convergence(int NN_reference, int Nt, int Nsim, const std::vector<int>& NN_approxs,
											double gamma)
{
	// Compute reference solution and noise, in addition
	// to creating the reference mesh. 

	// Set seed
	srand(1);

	// Nt \times Nsim 
	simulate_and_save_noise(NN_reference, Nt, Nsim); 

	// Save reference solution
	simulate(NN_reference, Nt, Nsim, gamma);

	Matrix<double, Dynamic,Dynamic> u_reference;
	u_reference.resize(vertex_value.size(), Nsim);
	u_reference.setZero();
	u_reference = u;

	// Write reference solution to see
	std::string path_vtk = vtk_animation_path+std::to_string(NN_reference)+".vtk";
	write_vtk_file(path_vtk);

	std::cout << "Computed and saved reference solution" << std::endl;


	// Make mesh matrix the identity
	mesh_matrix.resize(vertex_value.size(),vertex_value.size());
	mesh_matrix.setIdentity();

	std::vector<double> res;

	// This is now the number of time steps (different powers of 2)
	for (int NN_approx : NN_approxs)
	{
		// Simulate
		simulate(NN_reference, NN_approx, Nsim, gamma);

		// Write vtk file to see
		std::string path_vtk = vtk_animation_path+std::to_string(NN_approx)+".vtk";
		write_vtk_file(path_vtk);

		// Approximate the (relative) L2 error
		double e = approximate_L2_error(u_reference,u) / approximate_L2_norm(u_reference);
		res.push_back(e); 
	}

	return res;
}





