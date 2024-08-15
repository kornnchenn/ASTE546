#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include <chrono>
#include <memory>
#include <sstream>
#include <iomanip>
#include <cstring>

#include "PotentialSolver.h"
#include "Species.h" 
#include "Output.h"
#include "Source.h" //fine as is from Week 9?
#include "World.h" //


using namespace std;
using namespace Const;
using dvector = vector<double>;

int main(int argc, char *args[]) {
	
	//Initialize Simulation Domain
	World world(100,100,100); 

	const int ni = 100;	//number of nodes
	const double x0 = 0;	//origin
	const double xd = 5.35e-4;	//10 Debye Lengths in Meters
	double dx = (xd-x0)/(ni-1);	//node spacing

	/*3d
	World world(100,100,100); //at least 10 nodes per Debye Length
    world.setExtents({0,0,0},{5.35e-4,5.35e-4,5.35e-4});
	*/

	world.setTime(7.94e-13,10); //doing 10 time steps instead of 1000 for now	

	//define particle species
	vector<Species> species; // create a vector of type Species named species
	species.reserve(4); //define vector size
	species.push_back(Species("Deuterium+", 2*AMU, QE, 500, world));
	species.push_back(Species("Tritium+", 2*AMU, QE, 500, world));
	species.push_back(Species("Alpha", 2*AMU, QE, 500, world));
	species.push_back(Species("Electrons?")); //do I need to simulate electrons?
	//create pointers to the sim species
	Species &deuterium = species[0];
	Species &tritium = species[1];
	Species &alpha = species[2];

	//No Material Interaction due to plasma not touching walls in a Tokamak ?

	//Setup Injection Sources
	//Warm source?
	double den; //injection density <---- need to define

	sources.emplace_back(new WarmBeamSource(deuterium,world, 0.2, den ,100e3 )); //100keV temp
	sources.emplace_back(new WarmBeamSource(tritium,world, 0.2, den , 100e3)); //0.2 m/s drift velocity


	//Initialize Potential Solve and solve initial potential


	//Obtain initial EF


	//Simulation Main Loop
	while( world.advanceTime()) {


	}

}