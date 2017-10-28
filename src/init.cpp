#include "./../include/init.h"



Init::Init()
{

}

Init::~Init()
{
    //dtor
}

vector<int> Init::get_mols_each_dir()
{
	int mols_in_directions[] = MOLECULES_IN_EACH_DIRECTION;
	int max_mols_each_direction_size = SIZE_OF_ARRAY(mols_in_directions);
	if (max_mols_each_direction_size != DIMENSIONS)
	{
		cout << "dimensions of MOLECULES_IN_EACH_DIRECTION is " << max_mols_each_direction_size << " and it doesn't match with DIMENSIONS which is " << DIMENSIONS << "\n";
		exit(EXIT_FAILURE);
	}
	std::vector<int> molecules_in_each_directions(mols_in_directions, mols_in_directions + max_mols_each_direction_size);
	return molecules_in_each_directions;
}

vector<double> Init::get_sys_sizes()
{
	double sys_size[] = SYSTEM_SIZES;
	int max_sys_size = SIZE_OF_ARRAY(sys_size);
	if (max_sys_size != DIMENSIONS)
	{
		cout << "dimensions of SYSTEM_SIZE is " << max_sys_size << " and it doesn't match with DIMENSIONS which is " << DIMENSIONS << "\n";
		exit(EXIT_FAILURE);
	}
	vector<double> sys_sizes(sys_size, sys_size + max_sys_size);
	return sys_sizes;
}

vector< vector<double> > Init::get_col_location()
{
	vector< vector<double> > colloid_molecules;
#ifndef IGNORE_COLLOIDS
	int coll_mols[][DIMENSIONS] = COLLOID_MOLS;
	int max_coll_mols = SIZE_OF_ARRAY(coll_mols);
	for (int i = 0; i < max_coll_mols; i++)
	{
		vector<double> coll_location(coll_mols[i], coll_mols[i] + DIMENSIONS);
		colloid_molecules.push_back(coll_location);
	}
#endif // IGNORE_COLLOIDS
	return colloid_molecules;
}

vector<double> Init::get_temp_range()
{
	double temp_range[] = TEMPERATURE_RANGE;
	int temp_size = SIZE_OF_ARRAY(temp_range);
	vector<double> temperature_range(temp_range, temp_range + temp_size);
	return temperature_range;
}

vector<double> Init::get_init_spin()
{
	double init_spin[] = INIT_SPIN;
	int init_spin_size = SIZE_OF_ARRAY(init_spin);
	vector<double> initial_spin(init_spin, init_spin + init_spin_size);
	//check dimensions:
	if (init_spin_size != DIMENSIONS)
	{
		cout << "dimensions of INIT_SPIN is " << init_spin_size << " and it doesn't match with DIMENSIONS which is " << DIMENSIONS << "\n";
		exit(EXIT_FAILURE);
	}
	return initial_spin;
}

vector<Molecule> Init::get_molecules(const vector<int>& molecules_in_each_directions, const vector<double> & initial_spin)
{
	int max_mol = 1;
	for (unsigned int i = 0; i < DIMENSIONS; i++)
		max_mol *= molecules_in_each_directions[i];

	std::vector<Molecule> molecules(max_mol); //default initialize the vectors --with lc-mols

#if DIMENSIONS == 2
	for (int i = 0; i < molecules_in_each_directions[0]; i++)
	{
		for (int j = 0; j < molecules_in_each_directions[1]; j++)
		{
			int molecule_to_manipulate = molecules_in_each_directions[1] * i + j;
			molecules[molecule_to_manipulate].m_location[0] = i * INIT_SPACING;
			molecules[molecule_to_manipulate].m_location[1] = j * INIT_SPACING;
			molecules[molecule_to_manipulate].m_spin[0] = initial_spin[0];
			molecules[molecule_to_manipulate].m_spin[1] = initial_spin[1];
		}
	}
#elif DIMENSIONS == 3
	for (int i = 0; i < molecules_in_each_directions[0]; i++)
	{
		for (int j = 0; j < molecules_in_each_directions[1]; j++)
		{
			for (int k = 0; k < molecules_in_each_directions[2]; k++)
			{
				int molecule_to_manipulate = molecules_in_each_directions[2] * molecules_in_each_directions[1] * i + molecules_in_each_directions[2] * j + k;
				molecules[molecule_to_manipulate].m_location[0] = i * INIT_SPACING;
				molecules[molecule_to_manipulate].m_location[1] = j * INIT_SPACING;
				molecules[molecule_to_manipulate].m_location[2] = k * INIT_SPACING;
				molecules[molecule_to_manipulate].m_spin[0] = initial_spin[0];
				molecules[molecule_to_manipulate].m_spin[1] = initial_spin[1];
				molecules[molecule_to_manipulate].m_spin[2] = initial_spin[2];
			}
		}
	}
#endif //DIMENSIONS
	return molecules;
}

vector<Molecule>& Init::add_colloids(vector<Molecule>& molecules, vector<vector<double>> colloid_location) {
	Molecule col_mol;
	vector<double> loc;
	for (vector<vector<double>>::iterator it = colloid_location.begin(); it != colloid_location.end(); it++) {
		col_mol.m_mol_type = COLLOID;
		loc = (*it);
#if DIMENSIONS == 2
		col_mol.m_location[0] = loc[0];
		col_mol.m_location[1] = loc[1];
		col_mol.m_spin[0] = 1;
		col_mol.m_spin[1] = 0;
#elif DIMENSIONS == 3
		col_mol.m_location[0] = loc[0];
		col_mol.m_location[1] = loc[1];
		col_mol.m_location[2] = loc[2];
		col_mol.m_spin[0] = 1;
		col_mol.m_spin[1] = 0;
		col_mol.m_spin[2] = 0;
#endif //DIMENSIONS
		molecules.push_back(col_mol);

		return molecules;
	}
}

BoundaryType Init::get_boundary_condition() {
	return static_cast<BoundaryType>(BOUNDARY);
}

int Init::get_range() {
	return RANGE;
}

void Init::add_randomization(vector<Molecule>& molecules, const vector<int> & molecules_in_each_directions, 
	const vector<double> & sys_sizes)
{
	///initiate the random generators:
	srand((unsigned int)time(0));
	std::default_random_engine loc_gen((unsigned int)time(0));
	std::normal_distribution<double> loc_dist(0.0, INIT_SPACING_STD);

	std::default_random_engine spin_gen((unsigned int)time(0));
	std::normal_distribution<double> spin_dist(0.0, INIT_SPIN_STD);

	for (unsigned int i = 0; i < molecules.size(); i++)
	{
#ifdef DONT_MOVE_COLS
		if (molecules[i].m_mol_type == COLLOID)
			return; //return; colloid are packed in the end of molecules vector
#endif //DONT_MOVE_COLS
		//change the location with std of location:
		for (int j = 0; j < DIMENSIONS; j++)
		{
			double suggested_init_loc;
			do
			{
				suggested_init_loc = molecules[i].m_location[j] + loc_dist(loc_gen);
			} while ((suggested_init_loc + 0.5 < 0.0001) || (sys_sizes[j] - 0.5 - suggested_init_loc < 0.0001));
			molecules[i].m_location[j] = suggested_init_loc;
		}

		//change spin with the std of the spin and normalize it:
		double norm_spin = 0;
		for (int j = 0; j < DIMENSIONS; j++)
		{
			molecules[i].m_spin[j] = molecules[i].m_spin[j] + spin_dist(spin_gen);
			norm_spin += (molecules[i].m_spin[j] * molecules[i].m_spin[j]);
		}
		norm_spin = sqrt(norm_spin);
		for (int j = 0; j < DIMENSIONS; j++)
		{
			molecules[i].m_spin[j] = molecules[i].m_spin[j] / norm_spin;
		}
	}
}
