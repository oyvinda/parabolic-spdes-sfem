#pragma once

inline double area_T(const Matrix<int, 3,1>& T)
{
	// Computes the area of triange T
	Matrix<double,1,3> v0 = vertex_value[T(0,0)],
		v1 = vertex_value[T(1,0)],
		v2 = vertex_value[T(2,0)];

	double res = 0.5 * (v1-v0).cross(v2-v0).norm();
	return res;
}

inline double phii(const Matrix<double, 1,3>& x, const Matrix<int,3,1>& T, int i)
{
	// Computes the value of \varphi_i at x \in T

	int x0=T(0,0),
		x1=T(1,0),
		x2=T(2,0);

	Matrix<double, 1,3> v0=vertex_value[x0],
		v1 = vertex_value[x1],
		v2 = vertex_value[x2],
		vs,vt,vx;

	if (i==x0)
	{
		vs = v1-v0;
		vt = v2-v0;
		vx = x-v0;
	}
	if (i==x1)
	{
		vs=v0-v1;
		vt=v2-v1;
		vx=x-v1;
	}
	if (i==x2)
	{
		vs=v0-v2;
		vt=v1-v2;
		vx=x-v2;
	}

	// Compute x in terms of s and t

	Matrix<double,2,2> M;
	M.setZero();
	M(0,0)=vs.dot(vs);
	M(1,1)=vt.dot(vt);
	M(0,1)=vs.dot(vt);
	M(1,0)=vt.dot(vs);

	Matrix<double,2,1> y;
	y.setZero();
	y(0,0)=vx.dot(vs);
	y(1,0)=vx.dot(vt);

	Matrix<double,2,1> s_t = M.inverse()*y;

	// phi(s,t) = 1. - s - t
	double res = 1. - s_t(0,0) - s_t(1,0);

	return res;
}

inline Matrix<double, 1,3> Dphii(const Matrix<double, 1,3>& x, const Matrix<int,3,1>& T, int i)
{
	// Computes the value of \nabla varphi_i at x \in T 

	int x0=T(0,0),
		x1=T(1,0),
		x2=T(2,0);

	Matrix<double, 1,3> v0=vertex_value[x0],
		v1 = vertex_value[x1],
		v2 = vertex_value[x2],
		v, vstar;

	Matrix<double,1,3> res;
	res.setZero();

	if (i==x0)
	{
		v = v2-v1;
		vstar=v1;
	}
	if (i==x1)
	{
		v=v2-v0;
		vstar=v0;
	}
	if (i==x2)
	{
		v=v1-v0;
		vstar=v0;
	}

	// vstar + alpha * v - x0, v = 0
	double alpha = (vertex_value[i] - vstar).dot(v) / v.dot(v);

	res = vertex_value[i] - vstar - alpha * v;
	res = res / res.dot(res);

	return res;
}

struct phiiphij
{
	// x \mapsto \varphi_i(x) \varphi_j(x)
	// Function to be integrated over triangle
	// to create mass matrix

	phiiphij() { }

	double operator()(const Matrix<double, 1,3>& x, const Matrix<int,3,1>& T, int i, int j) const
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

	double operator()(const Matrix<double, 1,3>& x, const Matrix<int,3,1>& T, int i, int j) const
	{
		return Dphii(x,T,i).dot(Dphii(x,T,j));
	}
};


template<typename bilinear_form> 
double integrate_T(const bilinear_form& f, const Matrix<int,3,1>& T, int i, int j)
{
	// Integrates f over T using quadrature.
	// f takes arguments x, T, i, j.
	
	int x0=T(0,0),
		x1=T(1,0),
		x2=T(2,0);

	Matrix<double, 1,3> v0=vertex_value[x0],
		v1 = vertex_value[x1],
		v2 = vertex_value[x2];

	/*
	double w01 = f(0.5*(v0+v1),T,i,j),
		w02 = f(0.5*(v0+v2),T,i,j),
		w12 = f(0.5*(v1+v2),T,i,j);

	double res = area_T(T)*(w01+w02+w12)/3.;
	*/
	double w1 = f((v0+v1+v2)/3.,T,i,j) * (-9./16.);

	double w2 = f((3.*v0+v1+v2)/5.,T,i,j) * (25./48.);
	double w3 = f((v0+3.*v1+v2)/5.,T,i,j) * (25./48.);
	double w4 = f((v0+v1+3.*v2)/5.,T,i,j) * (25./48.);

	double res = area_T(T) * (w1+w2+w3+w4);

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
