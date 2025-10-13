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

	double operator()(const Matrix<double,1,3>& x, const Matrix<int,3,1>& T, int i, int j) const
	{

		// Define A matrix with varying diffusion
		Matrix<double,3,3> vvt;
		vvt.setIdentity();

		// Lifted x
		double x1 = x(0,0),
			x2 = x(0,1),
			x3 = x(0,2);

		double gamma_x3 = 1.*std::cos(PI * x3 / 2.)*std::cos(PI * x3 / 2.);

		Matrix<double,1,3> x0;
		x0 << 1., 0., 0.;

		double eta = 2.e-1 + 0.*std::exp(-10.*(x-x0).norm());

		vvt.block(0,0,1,3) << eta+gamma_x3*x2*x2, -gamma_x3*x1*x2, 0.;
		vvt.block(1,0,1,3) << -gamma_x3*x1*x2, eta+gamma_x3*x1*x1, 0.;
		vvt.block(2,0,1,3) << 0., 0., eta;

		return Dphii(x,T,i).dot(vvt*Dphii(x,T,j).transpose());

	}

};

// Bilinear form for advection term
struct adv_term
{
	adv_term() { }
	double operator()(const Matrix<double, 1,3>& x, const Matrix<int,3,1>& T, int i, int j) const
	{
		Matrix<double,3,1> adv;
		adv.setZero();

		// Lifted x
		double x1 = x(0,0),
			x2 = x(0,1),
			x3 = x(0,2);

		double gamma_x3 = 1. * std::cos(PI * x3 / 2.)*std::cos(PI * x3 / 2.);
		
		adv << gamma_x3*x2, -gamma_x3*x1, 0.;

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

