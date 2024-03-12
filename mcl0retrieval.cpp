/**
 * (c) 2023
 * the p-one collaboration
 * 
 * @date Nov 29 2023
 * @author cweaver with contributions from vparrish
*/


// This script is mc only retrieval of triggered dom information 
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

//#include "cl_options.h"
#include "l1-mockup.h"


#include <icetray/I3Tray.h>
#include <dataclasses/ModuleKey.h>
#include <dataclasses/geometry/I3Geometry.h>
#include <dataclasses/physics/I3RecoPulse.h>
#include <dataclasses/physics/I3Particle.h>

//bleh lazy
#include <nuSQuIDS/marray.h>
// struct ModuleTrigger{
// 	ModuleKey module;
// 	unsigned int multiplicity;
// 	double time;
// };


//find the multiplicity of each module, given the time window
//at the moment this only reports the highest multiplicity trigger for each module, not all triggers
std::vector<ModuleTrigger> findModuleMultiplicities(std::map<ModuleKey, std::map<unsigned int,PulseQueue>> modules, double timeWindow){
	std::vector<ModuleTrigger> triggers;
	for(auto& [mk, pmts] : modules){
		//std::cout << "Module " << mk << " has pulses on " << pmts.size() << " PMTs\n";
		unsigned int maxMultFound=0;
		
		//skip pulses too small to trigger
		for(auto& [pmt, pulses] : pmts){
			if(pulses.next().GetCharge()<0.25){
				pulses.advance();
				if(pulses.empty())
					pmts.erase(pmt);
			}
		}
		
		ModuleTrigger trigger{mk,0,0};
		while(!pmts.empty()){
			//figure out the earliest pulse not yet considered as a trigger start	
			double startTime=std::numeric_limits<double>::infinity();
			unsigned int leadTube;
			for(const auto& [pmt, pulses] : pmts){
				if(pulses.next().GetTime()<startTime){
					leadTube=pmt;
					startTime=pulses.next().GetTime();
				}
			}
			//std::cout << "  Start time is now " << startTime << " with PMT " << leadTube << " as leader\n";
			//check all PMTs to count those with a pulse in the window
			//this will always be at least 1 for the leadTube
			unsigned int mult=0;
			for(const auto& [pmt, pulses] : pmts){
				if(pulses.next().GetTime()<startTime+timeWindow){
					//std::cout << "      Pulse at " << pulses.next().GetTime() << " on PMT " << pmt << " is in\n";
					mult++;
				}
				if(pmt==leadTube){
					assert(pulses.next().GetTime()<startTime+timeWindow);
					assert(mult>0);
				}
			}
			//std::cout << "    Multiplicity in window: " << mult << '\n';
			if(mult>trigger.multiplicity){
				trigger.multiplicity=mult;
				trigger.time=startTime;
			}
			pmts[leadTube].advance(); //move past the start pulse used in this round
			if(pmts[leadTube].empty()) //if pulse series consumed, remove
				pmts.erase(leadTube);
		}
		if(trigger.multiplicity)
			triggers.push_back(trigger);
	}
	return triggers;
}


// Filtering function
std::vector<ModuleTrigger> levelZeroTriggers(const std::vector<ModuleTrigger>& triggers, unsigned int minMultiplicity) {
    std::vector<ModuleTrigger> levelZeroTriggers;

    for (const auto& trigger : triggers) {
        if (trigger.multiplicity > minMultiplicity) {
            // Add the trigger to the levelZeroTriggers vector
            levelZeroTriggers.push_back(trigger);
        }
    }

    return levelZeroTriggers;
}


int level0(int argc, char* argv[]){
	// l0 for pulses on different channels on same DOM for coincidence
	double timeWindow=10; //ns
	unsigned int moduleRequirement=1;
	
	std::string moduleTrigMultFile="/dev/stdout";
	std::string moduleTrigEffFile="/dev/stdout";

	OptionParser op;
	op.addOption({"w","window"},timeWindow,"Length of the coincidence time window in nanoseconds");
	op.addOption({"doms","modules"},moduleRequirement,"Minimum number of modules which must have PEs for an event to be considered");
	op.addOption("modTrigMultOutput",moduleTrigMultFile,"Path to the output file for module trigger efficiency data as function of multiplicity");
	op.addOption("modTrigEffOutput",moduleTrigEffFile,"Path to the output file for effective module trigger efficiency data");
	std::vector<std::string> positionals = op.parseArgs(argc, argv);
	if(op.didPrintUsage())
		return 0;
	
	if(positionals.size()<2){
		std::cerr << "No input files" << std::endl;geome
		return 1;
	}
	std::vector<std::string> inputFiles(positionals.begin()+1,positionals.end());
	
	unsigned int viableEvents=0;
	unsigned int lostEvents=0;
	double totalReadoutTime=0;
	unsigned int minMultiplicity=1;
	
	I3::init_icetray_lib();
	I3Tray tray;
	tray.AddModule("I3Reader","Reader")("FileNameList",inputFiles);
	auto getTriggerMultiplicity=[&](boost::shared_ptr<I3Frame> frame)->void{
		const std::string muonName="MCMuon";
		const std::string pulsesName="PMTResponse_nonoise";
		if(!frame->Has(muonName) || !frame->Has(pulsesName))
			return;
		const I3Map<OMKey, std::vector<I3RecoPulse>>& pulses=frame->Get<I3Map<OMKey, std::vector<I3RecoPulse>>>(pulsesName);
		const I3Particle& muon=frame->Get<I3Particle>(muonName);
		if(!geometry){
			if(frame->Has("I3Geometry"))
				geometry=frame->Get<boost::shared_ptr<const I3Geometry>>("I3Geometry");
			else
				return;
		}

		
		std::map<ModuleKey, std::map<unsigned int,PulseQueue>> modules;
		
		//group all data by module
		for(const auto& p : pulses)
			modules[ModuleKey(p.first.GetString(), p.first.GetOM())].insert(
				std::make_pair(p.first.GetPMT(), PulseQueue{p.second.begin(), p.second.end()}));
		
		if(modules.size()<moduleRequirement)
			return;
		viableEvents++;
		
		// modules that have PMT hit
		std::vector<ModuleTrigger> moduleTriggers=findModuleMultiplicities(std::move(modules), timeWindow);
		
		// std::sort(moduleTriggers.begin(),moduleTriggers.end(),multiplicityDescending);

		// modules that have passed l0 trigger
		std::vector<ModuleTrigger> levelZeroTriggers = levelZeroTriggers(triggers, minMultiplicity);
		
		// unsigned int maxMult=0;
		// unsigned int nTriggers=0;
		// for(const auto& trigger : moduleTriggers){
		// 	if(trigger.multiplicity>maxMult)
		// 		maxMult=trigger.multiplicity;
		// 	if(trigger.multiplicity>1)
		// 		nTriggers++;
		// }
				
		
		totalReadoutTime+=totalReadoutTimeMcElroy(moduleTriggers);
	};
	tray.AddModule("I3NullSplitter","NullSplit");
	tray.AddModule(getTriggerMultiplicity,"GetTMult");
	tray.Execute();



}