#pragma once

// Object to save current mesh and the nodal basis
struct mesh
{
	// Vertex number n -> value in R^3
	std::vector< Matrix<double, 1,3>> vertex_value_;

	// Triangle number n -> mesh triangle in (vertex idx)^3 
	std::vector< Matrix<int,3,1>> cells_;

	 // Vertex number n -> of vector idx of adjacent triangles
	std::vector< std::vector<int>> adj_;

	// Vertex number to adjacent vertex numbers
	std::vector< std::vector<int>> adj_n_;

	// Interior vertices
	std::vector<int> interior_;

	// Indices (interior first, then boundary)
	std::vector<int> indices_;

	mesh() : vertex_value_(vertex_value), cells_(cells), adj_(adj), adj_n_(adj_n), interior_(interior), indices_(indices) { }

	inline double phii_(const Matrix<double, 1,3>& x, const Matrix<int,3,1>& T, int i) const
	{
		// Computes the value of \varphi_i at x \in T

		int x0=T(0,0),
			x1=T(1,0),
			x2=T(2,0);

		Matrix<double, 1,3> v0=vertex_value_[x0],
			v1 = vertex_value_[x1],
			v2 = vertex_value_[x2],
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

};

void create_mesh_matrix(int from, int to, bool read = true)
{
	// Creates the mesh matrix
	assert(to < from);

	std::string path = "some-data/from-"+std::to_string(from)+"-to-"+std::to_string(to)+".txt";

	std::cout << "Started making mesh matrix:" << std::endl;
	
	// Finer mesh
	create_plane_mesh(from);
	reorder_vertices();
	mesh reference_mesh;

	// Coarser mesh
	create_plane_mesh(to);
	reorder_vertices();


	mesh_matrix.resize(interior.size(), reference_mesh.interior_.size());
	mesh_matrix.setZero();

	std::vector<Triplet<double>> coeff;

	if (read)
	{
		std::ifstream fin(path);

		int i, j;
		double val;
		while (fin >> i)
		{
			fin >> j;
			fin >> val;

			coeff.push_back(Triplet<double>(i,j,val));
		}

	} else {

		std::ofstream fout(path);
		fout << std::fixed << std::setprecision(20);

		for (int j = 0; j < reference_mesh.interior_.size(); j++)
		{
			// Print to see progress. 
			if (j % 2000==0)
				std::cout << j/(double)reference_mesh.interior_.size()*100. << std::endl; 

			auto x = reference_mesh.vertex_value_[j];

			// Iterate over triangles
			for (auto T : cells)
			{
					// Check if x is inside T
					int x0=T(0,0),
						x1=T(1,0),
						x2=T(2,0);

					Matrix<double, 1,3> v0=vertex_value[x0],
						v1 = vertex_value[x1],
						v2 = vertex_value[x2],
						vs,vt,vx;

					vs = v1-v0;
					vt = v2-v0;
					vx = x-v0;

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


					// Only nescessary for sphere mesh
					double EPS = 0.;
#ifdef sphere_mesh
					EPS = 0.1;
#endif
#ifdef icosasphere_mesh
					EPS = 0.1;
#endif

					// x is out of T
					if (s_t(0,0) < 0.-EPS || s_t(1,0) < 0.-EPS || s_t(0,0)+s_t(1,0) > 1.+EPS)
						continue;

#ifdef sphere_mesh
					// Only nescessary for sphere mesh
					if (x.dot(v0) < 0.-EPS || x.dot(v1) < 0.-EPS || x.dot(v2) < 0.-EPS) 
						continue;
#endif

#ifdef icosasphere_mesh
					// Only nescessary for sphere mesh
					if (x.dot(v0) < 0.-EPS || x.dot(v1) < 0.-EPS || x.dot(v2) < 0.-EPS) 
						continue;
#endif

					// x is inside. Compute coefficient

					if (x0 < interior.size())
					{
						double val = phii(x,T,x0);
						coeff.push_back(Triplet<double>(x0,j,val));
						// Write the data
						fout << x0 << ' ' << j << ' ' << val << '\n';
					}

					if (x1 < interior.size())
					{
						double val = phii(x,T,x1);
						coeff.push_back(Triplet<double>(x1,j,val));
						// Write the data
						fout << x1 << ' ' << j << ' ' << val << '\n';
					}

					if (x2 < interior.size())
					{
						double val = phii(x,T,x2);
						coeff.push_back(Triplet<double>(x2,j,val));
						// Write the data
						fout << x2 << ' ' << j << ' ' << val << '\n';
					}

					break;

			}
		}
	}

	// Fill obsmat
	mesh_matrix.setFromTriplets(coeff.begin(), coeff.end());

	std::cout << "Finished making mesh matrix." << std::endl;
}
