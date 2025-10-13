#pragma once

inline double area_T(const Matrix<int, 2,1>& T)
{
	// Computes the area of triange T
	Matrix<double,1,2> v0 = vertex_value[T(0,0)],
					   v1 = vertex_value[T(1,0)];

	double res = (v1-v0).norm();
	return res;
}

inline double phii(const Matrix<double, 1,2>& x, const Matrix<int,2,1>& T, int i)
{
	// Computes the value of \varphi_i at x \in T

	int x0=T(0,0),
		x1=T(1,0);

	Matrix<double, 1,2> v0=vertex_value[x0],
		v1 = vertex_value[x1],
		vx;

	if (i==x0)
	{
		vx = x-v0;
	}
	if (i==x1)
	{
		vx = x-v1;
	}

	double res = 1. - vx.norm() / (v1-v0).norm();

	return res;
}

inline Matrix<double, 1,2> Dphii(const Matrix<double, 1,2>& x, const Matrix<int,2,1>& T, int i)
{
	// Computes the value of \nabla varphi_i at x \in T 

	int x0=T(0,0),
		x1=T(1,0);

	Matrix<double, 1,2> v0=vertex_value[x0],
		v1 = vertex_value[x1],
		v;

	Matrix<double,1,2> res;
	res.setZero();

	if (i==x0)
	{
		v = v1-v0;
	}
	if (i==x1)
	{
		v = v0-v1;
	}

	res = -v / v.norm() * 1. / v.norm();

	return res;
}

struct phiiphij
{
	// x \mapsto \varphi_i(x) \varphi_j(x)
	// Function to be integrated over triangle
	// to create mass matrix

	phiiphij() { }

	double operator()(const Matrix<double, 1,2>& x, const Matrix<int,2,1>& T, int i, int j) const
	{
		return phii(x,T,i)*phii(x,T,j);
	}
};

struct DphiiDphij
{
	// x \mapsto \nabla \varphi_i(x) \cdot \nabla \varphi_j(x).
	// Function to be integrated over triangle
	// to create stiffness matrix

	DphiiDphij() { }

	double operator()(const Matrix<double, 1,2>& x, const Matrix<int,2,1>& T, int i, int j) const
	{
		return Dphii(x,T,i).dot(Dphii(x,T,j));
	}
};


template<typename bilinear_form> 
double integrate_T(const bilinear_form& f, const Matrix<int,2,1>& T, int i, int j)
{
	// Integrates f over T using quadrature.
	// f takes arguments x, T, i, j.
	
	int x0=T(0,0),
		x1=T(1,0);

	Matrix<double,1,2> v0 = vertex_value[x0],
					   v1 = vertex_value[x1];

	Matrix<double,1,2> p1 = v0 + (std::sqrt(3./5.) * 0.5 + 0.5) * (v1-v0),
		p2 = v0 + (-std::sqrt(3./5.) * 0.5 + 0.5) * (v1-v0),
		p3 = v0 + 0.5 * (v1 - v0);
	double w01 = f(p1,T,i,j) * 5./18.,
		w02 = f(p2,T,i,j) * 5./18.,
		w03 = f(p3,T,i,j) * 8./18.;

	double res = area_T(T)*(w01 + w02 + w03);

	return res;
}

template<typename bilinear_from>
double integrate_ij(const bilinear_from& a, int i, int j)
{
	// Computes a(\varphi_i, \varphi_j) for some bilinear form a

	if (i==j)
	{
		double res = 0.;

		for (int idx : adj[i])
		{
			auto T = cells[idx];
			res = res + integrate_T(a,T,i,i);
		}

		return res;
	}

	if ( std::count(adj_n[i].begin(),adj_n[i].end(), j))
	{
		// i is adjacent to j
		double res = 0.;

		// they share trianges 
		for (int Tiidx : adj[i])
		{
			for (int Tjidx : adj[j])
			{
				if (Tiidx==Tjidx)
				{
					auto Ti = cells[Tiidx];
					res = res + integrate_T(a,Ti,i,j);
					break;
				}
			}
		}

		return res;
	}

	return 0.;
}
