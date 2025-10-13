#pragma once

// Vertex number n -> value in R^3
std::vector< Matrix<double, 1,3>> vertex_value;

// Triangle number n -> mesh triangle in (vertex idx)^3 
std::vector< Matrix<int,3,1>> cells;

 // Vertex number n -> of vector idx of adjacent triangles
std::vector< std::vector<int>> adj;

// Vertex number to adjacent vertex numbers
std::vector< std::vector<int>> adj_n;

// Interior vertices
std::vector<int> interior;

// Indices order after ordering (interior first, then boundary)
std::vector<int> indices;

#ifdef square_mesh

void create_plane_mesh(int n)
{
	// Reset variable values
	vertex_value = {};
	cells = {};
	adj = {};
	adj_n = {};

	const double Z = 0.;

	static double vdata[4][3] = {
	   {0., 1., Z}, {1., 1., Z}, {0., 0., Z}, {1., 0., Z}
	};

	static int tindices[2][3] = {
	   {0,1,2}, {2,1,3} };

	std::vector< Matrix<int, 3,1>> triangles,
		triangles_;

	// Add initial vertices
	int idx = 0;
	for (int i = 0; i < 4; i++)
	{ 
		vertex_value.push_back({vdata[i][0], vdata[i][1], vdata[i][2]});
		idx++;
	}

	// Initialize triangles
	for (int i = 0; i < 2; i++)
	{
		triangles.push_back({tindices[i][0], tindices[i][1], tindices[i][2]});
	}

	// To check if we have already added a vertex
	std::map<std::pair<int,int> , int> middle_vertice;
	
	// Compute icosasphere mesh iteratively
	for (int i = 0; i < n; i++)
	{

		for (int j = 0; j < triangles.size(); j++)
		{
			// Index of triangle vertices 
			int x0 = triangles[j](0,0),
				x1 = triangles[j](1,0),
				x2 = triangles[j](2,0),
				x01,x02,x12;


			// Values of new vertices
			if (!middle_vertice[{x0,x1}]) {
				vertex_value.push_back((vertex_value[x0] + vertex_value[x1])/2.);
				middle_vertice[{x0,x1}] = idx;
				middle_vertice[{x1,x0}] = idx;
				x01 = idx++;
			} else {
				x01 = middle_vertice[{x0,x1}];
			}

			if (!middle_vertice[{x0,x2}]) {
				vertex_value.push_back((vertex_value[x0] + vertex_value[x2])/2.);
				middle_vertice[{x0,x2}] = idx;
				middle_vertice[{x2,x0}] = idx;
				x02 = idx++;
			} else {
				x02 = middle_vertice[{x0,x2}];
			}

			if (!middle_vertice[{x1,x2}]) {
				vertex_value.push_back((vertex_value[x1] + vertex_value[x2])/2.);
				middle_vertice[{x1,x2}] = idx;
				middle_vertice[{x2,x1}] = idx;
				x12 = idx++;
			} else {
				x12 = middle_vertice[{x1,x2}];
			}


			// Values of new traingles
			triangles_.push_back({x0, x01, x02});
			triangles_.push_back({x01, x12, x1});
			triangles_.push_back({x02, x12, x2});
			triangles_.push_back({x02, x01, x12});


		}

		triangles = triangles_;
		triangles_ = {};

	}


	cells = triangles;

	// Make adjacency lists
	int N = vertex_value.size();

	adj = std::vector<std::vector<int>>(N);
	adj_n = std::vector<std::vector<int>>(N);

	std::map<std::pair<int,int>, int> has_edge;


	for (int i = 0; i < cells.size(); i++)
	{
		int x0 = cells[i](0,0),
			x1 = cells[i](1,0),
			x2 = cells[i](2,0);

		adj[x0].push_back(i);
		adj[x1].push_back(i);
		adj[x2].push_back(i);

		if (!has_edge[{x0,x1}]) {
			adj_n[x0].push_back(x1);
			adj_n[x1].push_back(x0);
			has_edge[{x0,x1}] = 1;
			has_edge[{x1,x0}] = 1;
		}
		if (!has_edge[{x0,x2}]) {
			adj_n[x0].push_back(x2);
			adj_n[x2].push_back(x0);
			has_edge[{x0,x2}] = 1;
			has_edge[{x2,x0}] = 1;
		}
		if (!has_edge[{x1,x2}]) {
			adj_n[x1].push_back(x2);
			adj_n[x2].push_back(x1);
			has_edge[{x1,x2}] = 1;
			has_edge[{x2,x1}] = 1;
		}

	}


}


// Reorder vertices 
void reorder_vertices()
{

	// Collect interior and boundary vertices
	interior = {};
	std::vector<int> boundary = {};

#define noflow_boundary
#ifdef noflow_boundary
	for (int i = 0; i < vertex_value.size(); i++)
		interior.push_back(i);
	return;
#endif

	for (int i = 0; i < vertex_value.size(); i++)
	{
		Matrix<double,1,3> x = vertex_value[i];
		double x1 = x(0,0),
			x2 = x(0,1);
		if (x1 < 1.e-12 || x1 > 1.-1.e-12 || x2 < 1.e-12 || x2 > 1.-1.e-12)
		{
			boundary.push_back(i);
			continue;
		}
		interior.push_back(i);
	}

	// Reorder vertex_value so that first indices are interior

	indices = interior;
	for (int i : boundary)
		indices.push_back(i);

	std::map<int,int> reordering;
	for (int i = 0; i < indices.size(); i++)
		reordering[indices[i]] = i;


	// Reorder vertex values
	auto vertex_value_ = vertex_value;
	for (int i = 0; i < vertex_value.size(); i++)
		vertex_value[i] = vertex_value_[indices[i]];

	// Reorder cells
	for (auto& cell : cells)
	{
		cell(0,0) = reordering[cell(0,0)];
		cell(1,0) = reordering[cell(1,0)];
		cell(2,0) = reordering[cell(2,0)];
	}

	// Reorder adj
	auto adj_ = adj;
	for (int i = 0; i < adj.size(); i++)
		adj[i] = adj_[indices[i]];

	// Reorder adj_n
	auto adj_n_ = adj_n;
	for (int i = 0; i < adj_n.size(); i++)
		adj_n[i] = adj_n_[indices[i]];
	for (auto& ind : adj_n)
	{
		for (int& k : ind)
			k = reordering[k];
	}


}

#endif


#ifdef sphere_mesh

void create_plane_mesh(int n)
{
	// Reset variable values
	vertex_value = {};
	cells = {};
	adj = {};
	adj_n = {};

	static double vdata[6][3] = {
		{1.,0.,0.}, {-1.,0.,0.}, 
		{0.,1.,0.}, {0.,-1.,0.},
		{0.,0.,1.}, {0.,0.,-1.}
	};

	static int tindices[8][3] = {
		{0,4,3}, {0,4,2}, {1,4,3}, {1,4,2},
		{0,5,3}, {0,5,2}, {1,5,3}, {1,5,2}
	};

	std::vector< Matrix<int, 3,1>> triangles,
		triangles_;

	// Add initial vertices
	int idx = 0;
	for (int i = 0; i < 6; i++)
	{ 
		vertex_value.push_back({vdata[i][0], vdata[i][1], vdata[i][2]});
		idx++;
	}

	// Initialize triangles
	for (int i = 0; i < 8; i++)
	{
		triangles.push_back({tindices[i][0], tindices[i][1], tindices[i][2]});
	}

	// To check if we have already added a vertex
	std::map<std::pair<int,int> , int> middle_vertice;
	
	// Compute icosasphere mesh iteratively
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < triangles.size(); j++)
		{
			// Index of triangle vertices 
			int x0 = triangles[j](0,0),
				x1 = triangles[j](1,0),
				x2 = triangles[j](2,0),
				x01,x02,x12;


			// Values of new vertices
			if (!middle_vertice[{x0,x1}]) {
				vertex_value.push_back((vertex_value[x0] + vertex_value[x1]).normalized());
				middle_vertice[{x0,x1}] = idx;
				middle_vertice[{x1,x0}] = idx;
				x01 = idx++;
			} else {
				x01 = middle_vertice[{x0,x1}];
			}

			if (!middle_vertice[{x0,x2}]) {
				vertex_value.push_back((vertex_value[x0] + vertex_value[x2]).normalized());
				middle_vertice[{x0,x2}] = idx;
				middle_vertice[{x2,x0}] = idx;
				x02 = idx++;
			} else {
				x02 = middle_vertice[{x0,x2}];
			}

			if (!middle_vertice[{x1,x2}]) {
				vertex_value.push_back((vertex_value[x1] + vertex_value[x2]).normalized());
				middle_vertice[{x1,x2}] = idx;
				middle_vertice[{x2,x1}] = idx;
				x12 = idx++;
			} else {
				x12 = middle_vertice[{x1,x2}];
			}

			// Values of new traingles
			triangles_.push_back({x0, x01, x02});
			triangles_.push_back({x01, x12, x1});
			triangles_.push_back({x02, x12, x2});
			triangles_.push_back({x02, x01, x12});

		}

		triangles = triangles_;
		triangles_ = {};
	}

	cells = triangles;

	// Make adjacency lists
	int N = vertex_value.size();

	adj = std::vector<std::vector<int>>(N);
	adj_n = std::vector<std::vector<int>>(N);

	std::map<std::pair<int,int>, int> has_edge;

	for (int i = 0; i < cells.size(); i++)
	{
		int x0 = cells[i](0,0),
			x1 = cells[i](1,0),
			x2 = cells[i](2,0);

		adj[x0].push_back(i);
		adj[x1].push_back(i);
		adj[x2].push_back(i);

		if (!has_edge[{x0,x1}]) {
			adj_n[x0].push_back(x1);
			adj_n[x1].push_back(x0);
			has_edge[{x0,x1}] = 1;
			has_edge[{x1,x0}] = 1;
		}
		if (!has_edge[{x0,x2}]) {
			adj_n[x0].push_back(x2);
			adj_n[x2].push_back(x0);
			has_edge[{x0,x2}] = 1;
			has_edge[{x2,x0}] = 1;
		}
		if (!has_edge[{x1,x2}]) {
			adj_n[x1].push_back(x2);
			adj_n[x2].push_back(x1);
			has_edge[{x1,x2}] = 1;
			has_edge[{x2,x1}] = 1;
		}
	}
}


// Reorder vertices 
void reorder_vertices()
{
	// Collect interior and boundary vertices
	interior = {};
	std::vector<int> boundary = {};

	for (int i = 0; i < vertex_value.size(); i++)
	{
		interior.push_back(i);
	}
	
}

#endif


#ifdef icosasphere_mesh

void create_plane_mesh(int n)
{
	
	// Reset variable values
	vertex_value = {};
	cells = {};
	adj = {};
	adj_n = {};

	const double X = 0.525731112119133606;
	const double Z = 0.850650808352039932;

	static double vdata[12][3] = {
	   {-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z}, {X, 0.0, -Z},
	   {0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},
	   {Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0}
	};

	static int tindices[20][3] = {
	   {0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},
	   {8,10,1}, {8,3,10}, {5,3,8}, {5,2,3}, {2,7,3},
	   {7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6},
		   {6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11} };

	std::vector< Matrix<int, 3,1>> triangles,
		triangles_;

	// Add initial vertices
	int idx = 0;
	for (int i = 0; i < 12; i++)
	{ 
		vertex_value.push_back({vdata[i][0], vdata[i][1], vdata[i][2]});
		idx++;
	}

	// Initialize triangles
	for (int i = 0; i < 20; i++)
	{
		triangles.push_back({tindices[i][0], tindices[i][1], tindices[i][2]});
	}

	// To check if we have already added a vertex
	std::map<std::pair<int,int> , int> middle_vertice;
	
	// Compute icosasphere mesh iteratively
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < triangles.size(); j++)
		{
			// Index of triangle vertices 
			int x0 = triangles[j](0,0),
				x1 = triangles[j](1,0),
				x2 = triangles[j](2,0),
				x01,x02,x12;


			// Values of new vertices
			if (!middle_vertice[{x0,x1}]) {
				vertex_value.push_back((vertex_value[x0] + vertex_value[x1]).normalized());
				middle_vertice[{x0,x1}] = idx;
				middle_vertice[{x1,x0}] = idx;
				x01 = idx++;
			} else {
				x01 = middle_vertice[{x0,x1}];
			}

			if (!middle_vertice[{x0,x2}]) {
				vertex_value.push_back((vertex_value[x0] + vertex_value[x2]).normalized());
				middle_vertice[{x0,x2}] = idx;
				middle_vertice[{x2,x0}] = idx;
				x02 = idx++;
			} else {
				x02 = middle_vertice[{x0,x2}];
			}

			if (!middle_vertice[{x1,x2}]) {
				vertex_value.push_back((vertex_value[x1] + vertex_value[x2]).normalized());
				middle_vertice[{x1,x2}] = idx;
				middle_vertice[{x2,x1}] = idx;
				x12 = idx++;
			} else {
				x12 = middle_vertice[{x1,x2}];
			}

			// Values of new traingles
			triangles_.push_back({x0, x01, x02});
			triangles_.push_back({x01, x12, x1});
			triangles_.push_back({x02, x12, x2});
			triangles_.push_back({x02, x01, x12});

		}

		triangles = triangles_;
		triangles_ = {};
	}

	cells = triangles;

	// Make adjacency lists
	int N = vertex_value.size();

	adj = std::vector<std::vector<int>>(N);
	adj_n = std::vector<std::vector<int>>(N);

	std::map<std::pair<int,int>, int> has_edge;

	for (int i = 0; i < cells.size(); i++)
	{
		int x0 = cells[i](0,0),
			x1 = cells[i](1,0),
			x2 = cells[i](2,0);

		adj[x0].push_back(i);
		adj[x1].push_back(i);
		adj[x2].push_back(i);

		if (!has_edge[{x0,x1}]) {
			adj_n[x0].push_back(x1);
			adj_n[x1].push_back(x0);
			has_edge[{x0,x1}] = 1;
			has_edge[{x1,x0}] = 1;
		}
		if (!has_edge[{x0,x2}]) {
			adj_n[x0].push_back(x2);
			adj_n[x2].push_back(x0);
			has_edge[{x0,x2}] = 1;
			has_edge[{x2,x0}] = 1;
		}
		if (!has_edge[{x1,x2}]) {
			adj_n[x1].push_back(x2);
			adj_n[x2].push_back(x1);
			has_edge[{x1,x2}] = 1;
			has_edge[{x2,x1}] = 1;
		}
	}
}

// Reorder vertices 
void reorder_vertices()
{
	// Collect interior and boundary vertices
	interior = {};
	std::vector<int> boundary = {};

	for (int i = 0; i < vertex_value.size(); i++)
	{
		interior.push_back(i);
	}
	
}

#endif


#ifdef dome_mesh

void create_plane_mesh(int n)
{
	// Reset variable values
	vertex_value = {};
	cells = {};
	adj = {};
	adj_n = {};

	static double vdata[5][3] = {
		{1.,0.,0.}, {-1.,0.,0.}, 
		{0.,1.,0.}, {0.,-1.,0.},
		{0.,0.,1.}
	};

	static int tindices[4][3] = {
		{0,4,3}, {0,4,2}, {1,4,3}, {1,4,2}
		//{0,5,3}, {0,5,2}, {1,5,3}, {1,5,2}
	};

	std::vector< Matrix<int, 3,1>> triangles,
		triangles_;

	// Add initial vertices
	int idx = 0;
	for (int i = 0; i < 5; i++)
	{ 
		vertex_value.push_back({vdata[i][0], vdata[i][1], vdata[i][2]});
		idx++;
	}

	// Initialize triangles
	for (int i = 0; i < 4; i++)
	{
		triangles.push_back({tindices[i][0], tindices[i][1], tindices[i][2]});
	}

	// To check if we have already added a vertex
	std::map<std::pair<int,int> , int> middle_vertice;
	
	// Compute icosasphere mesh iteratively
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < triangles.size(); j++)
		{
			// Index of triangle vertices 
			int x0 = triangles[j](0,0),
				x1 = triangles[j](1,0),
				x2 = triangles[j](2,0),
				x01,x02,x12;


			// Values of new vertices
			if (!middle_vertice[{x0,x1}]) {
				vertex_value.push_back((vertex_value[x0] + vertex_value[x1]).normalized());
				middle_vertice[{x0,x1}] = idx;
				middle_vertice[{x1,x0}] = idx;
				x01 = idx++;
			} else {
				x01 = middle_vertice[{x0,x1}];
			}

			if (!middle_vertice[{x0,x2}]) {
				vertex_value.push_back((vertex_value[x0] + vertex_value[x2]).normalized());
				middle_vertice[{x0,x2}] = idx;
				middle_vertice[{x2,x0}] = idx;
				x02 = idx++;
			} else {
				x02 = middle_vertice[{x0,x2}];
			}

			if (!middle_vertice[{x1,x2}]) {
				vertex_value.push_back((vertex_value[x1] + vertex_value[x2]).normalized());
				middle_vertice[{x1,x2}] = idx;
				middle_vertice[{x2,x1}] = idx;
				x12 = idx++;
			} else {
				x12 = middle_vertice[{x1,x2}];
			}

			// Values of new traingles
			triangles_.push_back({x0, x01, x02});
			triangles_.push_back({x01, x12, x1});
			triangles_.push_back({x02, x12, x2});
			triangles_.push_back({x02, x01, x12});

		}

		triangles = triangles_;
		triangles_ = {};
	}

	cells = triangles;

	// Make adjacency lists
	int N = vertex_value.size();

	adj = std::vector<std::vector<int>>(N);
	adj_n = std::vector<std::vector<int>>(N);

	std::map<std::pair<int,int>, int> has_edge;

	for (int i = 0; i < cells.size(); i++)
	{
		int x0 = cells[i](0,0),
			x1 = cells[i](1,0),
			x2 = cells[i](2,0);

		adj[x0].push_back(i);
		adj[x1].push_back(i);
		adj[x2].push_back(i);

		if (!has_edge[{x0,x1}]) {
			adj_n[x0].push_back(x1);
			adj_n[x1].push_back(x0);
			has_edge[{x0,x1}] = 1;
			has_edge[{x1,x0}] = 1;
		}
		if (!has_edge[{x0,x2}]) {
			adj_n[x0].push_back(x2);
			adj_n[x2].push_back(x0);
			has_edge[{x0,x2}] = 1;
			has_edge[{x2,x0}] = 1;
		}
		if (!has_edge[{x1,x2}]) {
			adj_n[x1].push_back(x2);
			adj_n[x2].push_back(x1);
			has_edge[{x1,x2}] = 1;
			has_edge[{x2,x1}] = 1;
		}
	}
}


// Reorder vertices 
void reorder_vertices()
{
	// Collect interior and boundary vertices
	interior = {};
	std::vector<int> boundary = {};

	for (int i = 0; i < vertex_value.size(); i++)
	{
		interior.push_back(i);
	}
	
}

#endif


#ifdef dome_w_dirichlet_mesh

void create_plane_mesh(int n)
{
	// Reset variable values
	vertex_value = {};
	cells = {};
	adj = {};
	adj_n = {};

	static double vdata[5][3] = {
		{1.,0.,0.}, {-1.,0.,0.}, 
		{0.,1.,0.}, {0.,-1.,0.},
		{0.,0.,1.}
	};

	static int tindices[4][3] = {
		{0,4,3}, {0,4,2}, {1,4,3}, {1,4,2}
		//{0,5,3}, {0,5,2}, {1,5,3}, {1,5,2}
	};

	std::vector< Matrix<int, 3,1>> triangles,
		triangles_;

	// Add initial vertices
	int idx = 0;
	for (int i = 0; i < 5; i++)
	{ 
		vertex_value.push_back({vdata[i][0], vdata[i][1], vdata[i][2]});
		idx++;
	}

	// Initialize triangles
	for (int i = 0; i < 4; i++)
	{
		triangles.push_back({tindices[i][0], tindices[i][1], tindices[i][2]});
	}

	// To check if we have already added a vertex
	std::map<std::pair<int,int> , int> middle_vertice;
	
	// Compute icosasphere mesh iteratively
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < triangles.size(); j++)
		{
			// Index of triangle vertices 
			int x0 = triangles[j](0,0),
				x1 = triangles[j](1,0),
				x2 = triangles[j](2,0),
				x01,x02,x12;


			// Values of new vertices
			if (!middle_vertice[{x0,x1}]) {
				vertex_value.push_back((vertex_value[x0] + vertex_value[x1]).normalized());
				middle_vertice[{x0,x1}] = idx;
				middle_vertice[{x1,x0}] = idx;
				x01 = idx++;
			} else {
				x01 = middle_vertice[{x0,x1}];
			}

			if (!middle_vertice[{x0,x2}]) {
				vertex_value.push_back((vertex_value[x0] + vertex_value[x2]).normalized());
				middle_vertice[{x0,x2}] = idx;
				middle_vertice[{x2,x0}] = idx;
				x02 = idx++;
			} else {
				x02 = middle_vertice[{x0,x2}];
			}

			if (!middle_vertice[{x1,x2}]) {
				vertex_value.push_back((vertex_value[x1] + vertex_value[x2]).normalized());
				middle_vertice[{x1,x2}] = idx;
				middle_vertice[{x2,x1}] = idx;
				x12 = idx++;
			} else {
				x12 = middle_vertice[{x1,x2}];
			}

			// Values of new traingles
			triangles_.push_back({x0, x01, x02});
			triangles_.push_back({x01, x12, x1});
			triangles_.push_back({x02, x12, x2});
			triangles_.push_back({x02, x01, x12});

		}

		triangles = triangles_;
		triangles_ = {};
	}

	cells = triangles;

	// Make adjacency lists
	int N = vertex_value.size();

	adj = std::vector<std::vector<int>>(N);
	adj_n = std::vector<std::vector<int>>(N);

	std::map<std::pair<int,int>, int> has_edge;

	for (int i = 0; i < cells.size(); i++)
	{
		int x0 = cells[i](0,0),
			x1 = cells[i](1,0),
			x2 = cells[i](2,0);

		adj[x0].push_back(i);
		adj[x1].push_back(i);
		adj[x2].push_back(i);

		if (!has_edge[{x0,x1}]) {
			adj_n[x0].push_back(x1);
			adj_n[x1].push_back(x0);
			has_edge[{x0,x1}] = 1;
			has_edge[{x1,x0}] = 1;
		}
		if (!has_edge[{x0,x2}]) {
			adj_n[x0].push_back(x2);
			adj_n[x2].push_back(x0);
			has_edge[{x0,x2}] = 1;
			has_edge[{x2,x0}] = 1;
		}
		if (!has_edge[{x1,x2}]) {
			adj_n[x1].push_back(x2);
			adj_n[x2].push_back(x1);
			has_edge[{x1,x2}] = 1;
			has_edge[{x2,x1}] = 1;
		}
	}
}


// Reorder vertices 
void reorder_vertices()
{

	// Collect interior and boundary vertices
	interior = {};
	std::vector<int> boundary = {};

	for (int i = 0; i < vertex_value.size(); i++)
	{
		Matrix<double,1,3> x = vertex_value[i];
		double x1 = x(0,0),
			x2 = x(0,1),
			x3 = x(0,2);
		if (x3 < 1.e-8)
		{
			boundary.push_back(i);
			continue;
		}
		interior.push_back(i);
	}

	// Reorder vertex_value so that first indices are interior

	indices = interior;
	for (int i : boundary)
		indices.push_back(i);

	std::map<int,int> reordering;
	for (int i = 0; i < indices.size(); i++)
		reordering[indices[i]] = i;


	// Reorder vertex values
	auto vertex_value_ = vertex_value;
	for (int i = 0; i < vertex_value.size(); i++)
		vertex_value[i] = vertex_value_[indices[i]];

	// Reorder cells
	for (auto& cell : cells)
	{
		cell(0,0) = reordering[cell(0,0)];
		cell(1,0) = reordering[cell(1,0)];
		cell(2,0) = reordering[cell(2,0)];
	}

	// Reorder adj
	auto adj_ = adj;
	for (int i = 0; i < adj.size(); i++)
		adj[i] = adj_[indices[i]];

	// Reorder adj_n
	auto adj_n_ = adj_n;
	for (int i = 0; i < adj_n.size(); i++)
		adj_n[i] = adj_n_[indices[i]];
	for (auto& ind : adj_n)
	{
		for (int& k : ind)
			k = reordering[k];
	}


}

#endif