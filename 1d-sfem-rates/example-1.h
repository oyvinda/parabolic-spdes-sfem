#pragma once

// System matrix (sparse)
SparseMatrix<double> massm, stiffnm, vvtm, advm;

// Square root of mass matrix
SparseMatrix<double> M_sqrt;

std::default_random_engine generator;
std::normal_distribution<double> distribution(0.,1.);


// Bilinear form with varying diffusion 
struct vvt_term
{
	vvt_term() { }

	double operator()(const Matrix<double,1,2>& x, const Matrix<int,2,1>& T, int i, int j) const
	{

		// Define A matrix with varying diffusion
		Matrix<double,2,2> vvt;
		vvt.setIdentity();

		return Dphii(x,T,i).dot(vvt*Dphii(x,T,j).transpose());
	}

};

// Bilinear form for advection term
struct adv_term
{
	adv_term() { }
	double operator()(const Matrix<double, 1,2>& x, const Matrix<int,2,1>& T, int i, int j) const
	{
		Matrix<double,2,1> adv;
		adv.setZero();

		return phii(x,T,i) * Dphii(x,T,j).dot(adv);
	}
};


void assemble_system()
{

	// Assemble sparse matrices

	// Number of vertices
	int N = interior.size();

	massm.resize(N,N);
	stiffnm.resize(N,N);
	vvtm.resize(N,N);

	// Non zero matrix coefficients
	std::vector<Triplet<double>> coeff_massm, 
						coeff_stiffnm,
						coeff_vvt;

  	// Functions for mass and stiffness matrices
	phiiphij massf;
	DphiiDphij stiffnf;
	vvt_term vvtf;
	adv_term advf;

	// Assemble matrices
	for (int i = 0; i < N; i++)
	{

		// Compute coefficients in the matrices
		Triplet<double> massm_aij,
			stiffn_aij,
			vvt_aij;

		massm_aij = Triplet<double>(i,i, integrate_ij(massf,i,i));
		coeff_massm.push_back(massm_aij);

		stiffn_aij = Triplet<double>(i,i, integrate_ij(stiffnf,i,i));
		coeff_stiffnm.push_back(stiffn_aij);

		vvt_aij = Triplet<double>(i,i, integrate_ij(vvtf,i,i) + integrate_ij(advf,i,i));
		coeff_vvt.push_back(vvt_aij);


		for (int j : adj_n[i])
		{
			// Omit if at Dirichlet boundary
			if (j >= interior.size())
				continue;

			massm_aij = Triplet<double>(i,j, integrate_ij(massf,i,j));
			coeff_massm.push_back(massm_aij);

			stiffn_aij = Triplet<double>(i,j, integrate_ij(stiffnf,i,j));
			coeff_stiffnm.push_back(stiffn_aij);

			vvt_aij = Triplet<double>(i,j, integrate_ij(vvtf,i,j) + integrate_ij(advf,i,j));
			coeff_vvt.push_back(vvt_aij);

		}
	}

	// Insert into matrices
	massm.setFromTriplets(coeff_massm.begin(), coeff_massm.end());
	stiffnm.setFromTriplets(coeff_stiffnm.begin(), coeff_stiffnm.end());
	vvtm.setFromTriplets(coeff_vvt.begin(), coeff_vvt.end());

}

