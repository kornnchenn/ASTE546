#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <memory>
#include "World.h"
#include "PotentialSolver.h"
#include "Species.h"
#include "Output.h"
#include "Source.h"
#include "Collisions.h"

using namespace std;		//to avoid having to write std::cout
using namespace Const;		//to avoid having to write Const::ME

/*program execution starts here*/
int main(int argc, char *args[])
{
	/*initialize domain*/
    World world(31,81,71);
    world.setExtents({0,-0.4,0},{0.3,0.4,0.7});
    world.setTime(2e-7,10000);

	/*set objects*/
	double phi_object = -100;		//set default
	if (argc>1)
		phi_object = atof(args[1]);	//convert argument to float
	cout<<"Object potential: "<<phi_object<<" V"<<endl;
    world.addBox({-0.1,-0.1,0},{0.1,0.1,0.1},phi_object);
    //world.addInlet(); is this going to be the thruster inlet??

	/*set up particle species*/
    vector<Species> species;
	Species Xenon("Xenon", 131.3*AMU, 0, 2e9, world);
	Species XenonIon("Xenon+", 131.3*AMU, QE, 5e8, world);
    species.push_back(Xenon);
	species.push_back(XenonIon);
	Species &neutrals = species[0];
	Species &ions = species[1];
	
	/*setup injection sources*/
	const double nda = 2e20;			//neutral density
	const double ndi = 1e17;			//mean ion density --- changed in step 5
	vector<unique_ptr<Source>> sources;
	sources.emplace_back(new WarmBeamSource(neutrals,world,{0,0,0.1}, 0.02, 0.1e-6 ,500, 300));  //neutral source
	sources.emplace_back(new WarmBeamSource(ions,world,{0,0,0.1}, 0.02, 0.8e-6 , 15000, 1e4));	 //ions source
	
	/*setup material interactions*/
	vector<unique_ptr<Interaction>> interactions;
	interactions.emplace_back(new MCC_CEX(ions,neutrals,world));
	

	/*initialize potential solver and solve initial potential*/
    PotentialSolver solver(world,SolverType::QN,1000,1e-4);
    solver.setReferenceValues(30,2.5,ndi); //30 V, 2.5eV
    solver.solve();

    /*obtain initial electric field*/
    solver.computeEF();

    /* main loop*/
	while(world.advanceTime())
    {
		/*inject particles*/
    	for (auto &source:sources)
    		source->sample();

    	/*perform material interactions*/
    	for (auto &interaction:interactions)  interaction->apply(world.getDt());

		/*move particles*/
		for (Species &sp:species)
		{
			sp.advance();
			sp.computeNumberDensity();
			sp.sampleMoments();
			sp.computeMPC();
		}

		/*compute charge density*/
		world.computeChargeDensity(species);

        /*update potential*/
       // solver.solve();

        /*obtain electric field*/
        solver.computeEF();

        /*update averages at steady state*/
        if (world.steadyState(species)) {
        	for (Species &sp:species)
        		sp.updateAverages();
        }

		/*screen and file output*/
        Output::screenOutput(world,species);
        Output::diagOutput(world,species);

		/*periodically write out results*/
        if (world.getTs()%50==0 || world.isLastTimeStep()) {
			Output::fields(world, species);
		//	Output::particles(world, species,10000);
        }

    }
	
	/* grab starting time*/
	cout<<"Simulation took "<<world.getWallTime()<<" seconds"<<endl;
	return 0;		//indicate normal exit
}
