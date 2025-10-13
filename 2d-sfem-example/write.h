#pragma once

void write_vtk_file(const std::string& path)
{

	std::ofstream fout(path);

	// Set precision
	fout << std::fixed << std::setprecision(5);

	// Write header
	fout << "# vtk DataFile Version 2.0\nUnstructured Grid Example\nASCII\nDATASET UNSTRUCTURED_GRID\n\n";


	// Write points
	fout << "POINTS " << vertex_value.size() << " float\n";
	for (int i = 0; i < vertex_value.size(); i++)
	{
		auto x = vertex_value[i];
		fout << x(0,0) << ' ' << x(0,1) << ' ' << x(0,2) << '\n';
	}
	fout << '\n';

	// Write cells
	fout << "CELLS " << cells.size() << ' ' << 4 * cells.size() << '\n';
	for (int i = 0; i < cells.size(); i++)
	{
		auto T = cells[i];
		fout << 3 << ' ' << T(0,0) << ' ' << T(1,0) << ' ' << T(2,0) << '\n'; 
	}
	fout << '\n';

	// Write cell type
	fout << "CELL_TYPES " << cells.size() << '\n';
	for (int i = 0; i < cells.size(); i++)
	{
		fout << 5 << '\n';
	}
	fout << '\n';

	// Write function value
	fout << "POINT_DATA " << vertex_value.size() << '\n';
	fout << "SCALARS scalars float 1\nLOOKUP_TABLE default\n";
	for (int i = 0; i < u.rows(); i++)
	{
		fout << -u(i,0) << '\n';
	}

}










