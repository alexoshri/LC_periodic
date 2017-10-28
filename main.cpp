#include "./include/mol_sys.h"
#include "./include/grid.h"
#include "./include/molecule.h"
#include "./include/defined.h"
#include "./include/init.h"

using namespace std;

int main(int argc, char* argv[])
{
	//initializing all the vectors:
	vector<double> sys_sizes = Init::get_sys_sizes();
	vector<int> molecules_in_each_directions = Init::get_mols_each_dir();
	vector<vector<double> > colloid_location = Init::get_col_location();
	vector<double> temperature_range = Init::get_temp_range();
	vector<double> initial_spin = Init::get_init_spin();
	vector<Molecule> molecules;
#ifdef SARTING_SYSTEM
	molecules = Init::get_molecules_from_file(molecules_in_each_directions);
#else
	molecules = Init::get_molecules(molecules_in_each_directions, initial_spin);
#endif // DEBUG

	
	int num_lc = molecules.size();
	molecules = Init::add_colloids(molecules, colloid_location);
	BoundaryType bc = Init::get_boundary_condition();
	int range = Init::get_range();

	//add randomization to the initial location and orientation of the vectors:
	Init::add_randomization(molecules, molecules_in_each_directions, sys_sizes);
	cout << "finished initialization, start colloing:" << endl;	
	//here we can call to some function that will modify the system if user want to.
	Mol_Sys * lc_system = new Mol_Sys(sys_sizes, molecules, num_lc, temperature_range,bc,range);
	lc_system->start_cooling();
	delete lc_system;

	cout << "finished successfully" << endl;
}

