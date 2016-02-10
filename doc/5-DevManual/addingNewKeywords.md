
Adding a keyword to the input file {#newKeywords}
=================================================


Follow the instructions to add a new keyword: 

 * Add the new item in "include/SimulationParams.h"
 * Initialise and destroy the new item in SimulationParams::SimulationParams()
      and SimulationParams::~SimulationParams()
 * Read in the value in SimulationParams::assign
 * Validate the value in SimulationParams::validate
 * Go to World::init and manage your new item to be used to initialise your new stuff.
 
