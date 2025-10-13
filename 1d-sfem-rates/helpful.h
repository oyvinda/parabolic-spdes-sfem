#pragma once

// Object to save current mesh and the nodal basis
struct mesh
{
	// Vertex number n -> value in R^2
	std::vector< Matrix<double, 1,2>> vertex_value_;

	// Triangle number n -> mesh triangle in (vertex idx)^2
	std::vector< Matrix<int,2,1>> cells_;

	// Vertex number n -> of vector idx of adjacent triangles
	std::vector< std::vector<int>> adj_;

	// Vertex number to adjacent vertex numbers
	std::vector< std::vector<int>> adj_n_;

	// Interior vertices
	std::vector<int> interior_;

	// Indices (interior first, then boundary)
	std::vector<int> indices_;

	mesh() : vertex_value_(vertex_value), cells_(cells), adj_(adj), adj_n_(adj_n), interior_(interior), indices_(indices) { }

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
						x1=T(1,0);

					Matrix<double, 1,2> v0 = vertex_value[x0],
										v1 = vertex_value[x1];


					// Only nescessary for sphere mesh
					double EPS = 0.;
#ifdef sphere_mesh
					EPS = 1.e-12;
#endif

					// x is out of T
					double varphix = (v1-v0).dot(x-v0) / ((v1-v0).norm()*(v1-v0).norm());
					if (varphix < 0.-EPS || varphix > 1.+EPS)
						continue;

#ifdef sphere_mesh
					// Only nescessary for sphere mesh
					if (x.dot(v0) < 0. || x.dot(v1) < 0.) 
						continue;
#endif

					// x is inside. Compute coefficient
					// Adjust x to make sure that it is on the line 
					auto x_ = v0 + (x-v0).dot(v1-v0)*(v1-v0) / (v1-v0).dot(v1-v0); 
					//assert((x-x_).norm() < 1.e-10);

					if (x0 < interior.size())
					{
						double val = phii(x_,T,x0);
						coeff.push_back(Triplet<double>(x0,j,val));
						// Write the data
						fout << x0 << ' ' << j << ' ' << val << '\n';
					}

					if (x1 < interior.size())
					{
						double val = phii(x_,T,x1);
						coeff.push_back(Triplet<double>(x1,j,val));
						// Write the data
						fout << x1 << ' ' << j << ' ' << val << '\n';
					}

					break;

			}
		}
	}

	// Fill obsmat
	mesh_matrix.setFromTriplets(coeff.begin(), coeff.end());

	std::cout << "Finished making mesh matrix." << std::endl;
}
