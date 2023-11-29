#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

#include "cl_options.h"

#include <icetray/I3Tray.h>
#include <dataclasses/ModuleKey.h>
#include <dataclasses/geometry/I3Geometry.h>
#include <dataclasses/physics/I3RecoPulse.h>
#include <dataclasses/physics/I3Particle.h>

//bleh lazy
#include <nuSQuIDS/marray.h>

struct PulseQueue{
	using PulseIterator=std::vector<I3RecoPulse>::const_iterator;
	PulseIterator it, end;
	
	bool empty() const{ return it==end; }
	const I3RecoPulse& next() const{ return *it; }
	void advance(){
		while(it!=end){
			it++;
			//stop if we find a pulse large enough to trigger
			if(it!=end && it->GetCharge()>=0.25)
				break;
		}
	}
};

struct ModuleTrigger{
	ModuleKey module;
	unsigned int multiplicity;
	double time;
};

bool multiplicityDescending(const ModuleTrigger& t1, const ModuleTrigger& t2){
	return t1.multiplicity>t2.multiplicity;
}

bool timeAscending(const ModuleTrigger& t1, const ModuleTrigger& t2){
	return t1.time<t2.time;
}

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

//fraction of the time when each multiplicity threshold is expected to be in use
std::vector<double> moduleMultiplicityDutyCycles={0,0,.409,.365,.129,.041,0,0,0,0,0,0,0,0,0,0};
//fraction of the time when the multiplicity threshold is smaller than or equal to a count
std::vector<double> moduleMultiplicityThresholdProbabilities(moduleMultiplicityDutyCycles.size(),0);

namespace{
struct threshMultInit{
	threshMultInit(){
		for(std::size_t i=1; i<moduleMultiplicityDutyCycles.size(); i++)
			moduleMultiplicityThresholdProbabilities[i]=moduleMultiplicityThresholdProbabilities[i-1]+moduleMultiplicityDutyCycles[i];
	}
} tmInitializer;
}

auto triggerSuccessProbability(unsigned int multiplicity)->double{
	if(multiplicity>=moduleMultiplicityThresholdProbabilities.size())
		return moduleMultiplicityThresholdProbabilities.back();
	return moduleMultiplicityThresholdProbabilities[multiplicity];
}

//compute the periods of time over which all module triggers occur
double moduleTriggerTimeSpan(const std::vector<ModuleTrigger>& moduleTriggers, unsigned int minimumThreshold){
	double minTime=std::numeric_limits<double>::infinity();
	double maxTime=-std::numeric_limits<double>::infinity();
	for(const auto& trigger : moduleTriggers){
		if(trigger.multiplicity<minimumThreshold)
			continue;
		if(trigger.time<minTime)
			minTime=trigger.time;
		if(trigger.time>maxTime)
			maxTime=trigger.time;
	}
	if(maxTime<minTime)
		return 0;
	return maxTime-minTime;
}

double windowscale=1.2;
double attenlength=35; //meters
double maxphotonatten=2.3; //attenuation lengths => 10% survival
double MaxTracktoReferenceDOM=50; //meters

const double c=0.299792458; //m/ns
const double n=1.34744898746;
const double ngroup=1.376924748647469;
const double theta_c=acos(1/n); //rad

std::pair<double,double> tMcElroyWindow(double dis){
	if(dis==0)
		return(std::make_pair(0.,300.));
	double dis_atten=attenlength*maxphotonatten; //distance photons can travel before considered 'fully' attenuated
	double dis_maxphoton=std::min(dis,dis_atten);
	double dis_part=0; //distance to particle? From reference/seed point?
	if(dis>dis_atten)
		dis_part=sin(theta_c-asin(sin(theta_c)*(dis_atten/dis)))*dis/sin(theta_c);
	double windowMin=std::max(0.,(dis/c)-MaxTracktoReferenceDOM*ngroup/c)/windowscale;
	double windowMax=windowscale*(dis_part/c + dis_maxphoton*ngroup/c);
	return std::make_pair(windowMin,windowMax);
}

boost::shared_ptr<const I3Geometry> geometry;

//assumes moduleTriggers is sorted with multiplicityDescending
double totalReadoutTimeMcElroy(const std::vector<ModuleTrigger>& moduleTriggers){
	if(!geometry)
		throw std::runtime_error("No geometry!");
	if(moduleTriggers.empty())
		return 0;
	auto modulePosition=[&](ModuleKey m){
		auto geoIt=geometry->omgeo.find(OMKey(m.GetString(),m.GetOM(),1));
		if(geoIt==geometry->omgeo.end())
			log_fatal_stream("Module " << m << " not in geometry?");
// 			throw std::runtime_error("Module not in geometry?");
		return geoIt->second.position;
	};
	
	auto seedModule=moduleTriggers.front().module;
	I3Position seedPos=modulePosition(seedModule);
	double seedTime=moduleTriggers.front().time;
	
	double readoutTime=0;
	std::set<ModuleKey> modulesSeen;
	for(auto& [omkey, omgeo] : geometry->omgeo){
		ModuleKey mk(omkey.GetString(),omkey.GetOM());
		if(modulesSeen.count(mk))
			continue;
		modulesSeen.insert(mk);
		double distance=(seedPos-omgeo.position).Magnitude();
		auto timeRange=tMcElroyWindow(distance);
		readoutTime+=2*(timeRange.second-timeRange.first);
	}
	return readoutTime;
}

//computes the probability that this event would be able to trigger the detector globally, assuming
//that the noise rate is fully correlated in all modules, such that they always share the same
//threshold for module triggering, and two module triggers are required for a global trigger
//assumes moduleTriggers is sorted with multiplicityDescending
double globalTriggerProbabilityFullCorrelation(const std::vector<ModuleTrigger>& moduleTriggers){
	double totalProb=0;
	for(unsigned int multThresh=2; multThresh<moduleMultiplicityDutyCycles.size() && moduleMultiplicityDutyCycles[multThresh];multThresh++){
		unsigned int triggersAboveThresh=0;
		for(; triggersAboveThresh<moduleTriggers.size(); triggersAboveThresh++){
			if(moduleTriggers[triggersAboveThresh].multiplicity<multThresh)
				break;
		}
		bool triggers=triggersAboveThresh>=2;
		if(triggers)
			totalProb+=moduleMultiplicityDutyCycles[multThresh];
	}
	return totalProb;
	
	/*for(const auto& t1 : moduleTriggers){
		for(const auto& t2 : moduleTriggers){
			if(t1.module==t2.module)
				continue; //no self-coincidence
			if(t1.module.GetString()==t1.module.GetString() && 
			   std::abs((int)t1.module.GetOM()-(int)t2.module.GetOM())==1 &&
			   (t1.multiplicity>1 || t2.multiplicity>1)){
				if(t1.multiplicity==1)
					return triggerSuccessProbability(t2.multiplicity);
				else
					return triggerSuccessProbability(t1.multiplicity);
			}
		}
	}
	return 0;*/
}

int main(int argc, char* argv[]){
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
		std::cerr << "No input files" << std::endl;
		return 1;
	}
	std::vector<std::string> inputFiles(positionals.begin()+1,positionals.end());
	
	
// 	auto triggerSuccessProbability=[&](unsigned int multiplicity)->double{
// 		if(multiplicity>=moduleMultiplicityThresholdProbabilities.size())
// 			return moduleMultiplicityThresholdProbabilities.back();
// 		return moduleMultiplicityThresholdProbabilities[multiplicity];
// 	};
	
	const double emin=100, espacing=10;
	nusquids::marray<double,2> hist({40,16});
	auto histInsert=[&](double energy, unsigned int multiplicity){
		if(energy<emin)
			return;
		unsigned int eIdx=(log10(energy)-log10(emin))*espacing;
		if(eIdx>hist.extent(0))
			return;
		if(multiplicity>=hist.extent(1))
			multiplicity=hist.extent(1)-1; //treat as an overflow bin
		hist[eIdx][multiplicity]+=1;
	};
	nusquids::marray<double,1> globalTriggerCorrelated({40});
	auto insertGlobalTriggerCorrelated=[&](double energy, double probability){
		if(energy<emin)
			return;
		unsigned int eIdx=(log10(energy)-log10(emin))*espacing;
		if(eIdx>globalTriggerCorrelated.extent(0))
			return;
		globalTriggerCorrelated[eIdx]+=probability;
	};
	unsigned int viableEvents=0;
	unsigned int lostEvents=0;
	double totalReadoutTime=0;
	
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
// 		if(muon.GetEnergy()<10000 || muon.GetEnergy()>12589.3)
// 			return;
		
		std::map<ModuleKey, std::map<unsigned int,PulseQueue>> modules;
		
		//group all data by module
		for(const auto& p : pulses)
			modules[ModuleKey(p.first.GetString(), p.first.GetOM())].insert(
				std::make_pair(p.first.GetPMT(), PulseQueue{p.second.begin(), p.second.end()}));
			//modules[ModuleKey(p.first.GetString(), 0)].insert(
			//	std::make_pair(p.first.GetOM()*20+p.first.GetPMT(), PulseQueue{p.second.begin(), p.second.end()}));
		
		if(modules.size()<moduleRequirement)
			return;
		viableEvents++;
		
		std::vector<ModuleTrigger> moduleTriggers=findModuleMultiplicities(std::move(modules), timeWindow);
		
		std::sort(moduleTriggers.begin(),moduleTriggers.end(),multiplicityDescending);
		
		unsigned int maxMult=0;
		unsigned int nTriggers=0;
		for(const auto& trigger : moduleTriggers){
			if(trigger.multiplicity>maxMult)
				maxMult=trigger.multiplicity;
			if(trigger.multiplicity>1)
				nTriggers++;
		}
		
// 		//look for multiplicity 1 triggers adjacent to higher multiplicities
//		//this is not very useful, as it only contributes to the calculation of lostEvents, 
//		//not the global trigger probability
// 		for(const auto& trigger : moduleTriggers){
// 			if(trigger.multiplicity>1)
// 				continue;
// 			bool nearby=false;
// 			if(trigger.module.GetOM()>0){
// 				ModuleKey below(trigger.module.GetString(), trigger.module.GetOM()-1);
// 				for(const auto& otherTrigger : moduleTriggers){
// 					if(otherTrigger.module==below && otherTrigger.multiplicity>1){
// 						nearby=true;
// 						break;
// 					}
// 				}
// 			}
// 			if(!nearby && trigger.module.GetOM()<20){
// 				ModuleKey above(trigger.module.GetString(), trigger.module.GetOM()+1);
// 				for(const auto& otherTrigger : moduleTriggers){
// 					if(otherTrigger.module==above && otherTrigger.multiplicity>1){
// 						nearby=true;
// 						break;
// 					}
// 				}
// 			}
// 			if(nearby){
// 				std::cout << "Nearby\n";
// 				nTriggers++;
// 				break;
// 			}
// 		}
				
// 		std::cout << "Energy: " << muon.GetEnergy() 
// 		<< " Maximum coincidence multiplicity: " << maxMult 
// 		<< " Number of LC module triggers: " << nTriggers
// 		<< '\n';
		if(maxMult>1 && nTriggers<2){
			lostEvents++;
		}
		if(maxMult)
			histInsert(muon.GetEnergy(),maxMult-1);
// 		if(maxMultString)
// 			histInsert(muon.GetEnergy(),maxMultString-1);
		//std::cout << "----\n";
		
		double triggerSpan=moduleTriggerTimeSpan(moduleTriggers,2);
		double globalProbFullCorr=globalTriggerProbabilityFullCorrelation(moduleTriggers);
		
// 		std::cout << "  Module trigger multiplicities:";
// 		for(const auto& trigger : moduleTriggers)
// 			std::cout << ' ' << trigger.multiplicity;
// 		std::cout << '\n';
// 		std::cout << "  Module triggers span " << triggerSpan << " ns\n";
// 		std::cout << "  Probability of global trigger: " << globalProbFullCorr << '\n';
		
		insertGlobalTriggerCorrelated(muon.GetEnergy(),globalProbFullCorr);
		
		totalReadoutTime+=totalReadoutTimeMcElroy(moduleTriggers);
	};
	tray.AddModule("I3NullSplitter","NullSplit");
	tray.AddModule(getTriggerMultiplicity,"GetTMult");
	tray.Execute();
	
	std::cout << "# Global trigger would lose " << lostEvents << " events out of " << viableEvents << std::endl;
	std::cout << "Average readout time per event: " << totalReadoutTime/viableEvents << " ns\n";
	
	std::cout << "----\n";
	std::ofstream moduleTrigMultData(moduleTrigMultFile);
	moduleTrigMultData << "#Emin\tEmax\tNEvents\tM=1     \tM=2     \tM=3     \tM=4     \tM=5" << std::endl;
	nusquids::marray<double,1> energySum({hist.extent(0)});
	for(unsigned int eIdx=0; eIdx<hist.extent(0); eIdx++){
		double elow=emin*pow(10.,double(eIdx)/espacing);
		double ehigh=emin*pow(10.,double(eIdx+1)/espacing);
		//accumulate in multiplicity
		double sum=0;
		for(int m=hist.extent(1)-1; m>=0; m--){
			sum+=hist[eIdx][m];
			hist[eIdx][m]=sum;
		}
		//print normalized
		if(sum){
			moduleTrigMultData << elow << '\t' << ehigh << '\t' << (int)sum;
			for(int m=0; m<std::min(5UL,hist.extent(1)); m++)
				moduleTrigMultData << " \t" << std::setw(8) << hist[eIdx][m]/sum;
			moduleTrigMultData << std::endl;
		}
		energySum[eIdx]=sum;
	}
	
	std::cout << "----\n";
	std::ofstream moduleTrigEffData(moduleTrigEffFile);
	moduleTrigEffData << "#Emin\tEmax\tP(module trigger)" << std::endl;
	for(unsigned int eIdx=0; eIdx<hist.extent(0); eIdx++){
		double elow=emin*pow(10.,double(eIdx)/espacing);
		double ehigh=emin*pow(10.,double(eIdx+1)/espacing);
		
		if(!energySum[eIdx])
			continue;
		
		double sum=0;
		for(int m=1; m<hist.extent(1); m++)
			sum+=moduleMultiplicityDutyCycles[m+1]*hist[eIdx][m]/energySum[eIdx];
		
		moduleTrigEffData << elow << '\t' << ehigh << '\t' << sum << std::endl;
	}
	
	std::cout << "----\n";
	std::cout << "#Emin\tEmax\tP(global trigger)" << std::endl;
	for(unsigned int eIdx=0; eIdx<globalTriggerCorrelated.extent(0); eIdx++){
		double elow=emin*pow(10.,double(eIdx)/espacing);
		double ehigh=emin*pow(10.,double(eIdx+1)/espacing);
		
		if(!energySum[eIdx])
			continue;
		
		std::cout << elow << '\t' << ehigh << '\t' << globalTriggerCorrelated[eIdx]/energySum[eIdx] << std::endl;
	}
}
