/**
 * (c) 2023
 * the p-one collaboration
 * 
 * @date Nov 29 2023
 * @author vparrish
*/

#include "src/l1-mockup.h"
#include "src/mcl0retrival.cpp"


// time comparison of triggered modules
bool compareTriggerTimes(const ModuleTrigger& a, const ModuleTrigger& b) {
    return a.time < b.time;
}

// retrieve vector of passed l0 triggers and order them in time
void timeOrderedL0Triggers(const std::vector<ModuleTrigger>& l0_triggers){
    std::vector<ModuleTrigger> orderedLevelZeroTriggers = l0_triggers;
    std::sort(orderedLevelZeroTriggers.begin(), orderedLevelZeroTriggers.end(), compareTriggerTimes);
}

/** use light cone geometry and different signal track projections to select modules of interest in event
 * currently assuming fixed vector of events (obv in triggering, this will be a live stream)
 * the analogous documentation outlining this geometric light-cone calculation 
*/
// make this a map of module keys to a vector 


//what does this do...
if(!geometry){
if(frame->Has("I3Geometry"))
    geometry=frame->Get<boost::shared_ptr<const I3Geometry>>("I3Geometry");
else
    return;
}

boost::shared_ptr<const I3Geometry> geometry;
//this is where I need to implement the light cone algorithms
//<moduletrigger> has modulekey, mult, and time 
double levelOneLightConeReconstruction(const ModuleTrigger& l0_triggers_i, const ModuleTrigger& pulse){

    //std::map<ModuleKey, std::vector<EventCluster> &eventCluster> l1_skimmedEventClusters;
    //define variables that we need for the comparison 
    double dist;
    double t_initial; 
    double tmin; 
    double tmax;
    double d_max; 
    
    //get distance between initial module and the potential neighbor  
	auto geoIt1=geometry->omgeo.find(OMKey(l0_triggers_i.GetString(),l0_triggers_i.GetOM(),1));
	double d1 = geoIt1->second.position;

    auto geoIt2=geometry->omgeo.find(OMKey(pulse.GetString(),pulse.GetOM(),1));
	double d2 = geoIt2->second.position;

    dist = sqrt(pow(d1[0]-d2[0], 2), pow(d1[1]-d2[1], 2), pow(d1[2]-d2[2], 2)); 
    
    //also need to get the time of the initial pulse
    t_initial = pulse.time();

    //const variables 
    double  c = 0.299792458; //m/ns
    double n = 1.34; //index of refraction
    double d_atten = 25; //i don't actually know what this is but it's close to 25m
    double theta_c = 40.5; //idk what it is but its like around 40 I think 

    //start the neighbor calculations 
    //is d_max attenuation length? attenuation length is used in the formulas tho   

    // lower limit
    if (dist < 2*d_atten*sin(theta_c)){
        tmin = 0;
    }else {
        double tmin_int = 1/c*sqrt(pow(dist, 2)- (2*d_atten*sin(theta_c)));
        double tmin_large = (d_atten/c)*((1/n) + sqrt((pow(dist/d_atten,2))-pow(sin(theta_c), 2)) - n);
        //determine if we're in the intermediate or large range- computer both eq34 and eq48 and use min of two
        if(tmin_int < tmin_large){
            tmin = tmin_int;
        } else {
            tmin = tmin_large;
        }
    }

    //upper limit
    //I actually don't know what this distance condition should be here??
    if(dist < 2*d_atten*sin(theta_c)){
        tmax = (n*dist)/c;
    } else {
        tmax = (d_atten/c)*(sqrt((pow(dist/d_atten,2))-pow(sin(theta_c), 2)) - (1/n) + n);
    }
    
    if(tmax!=0){
        //add t_min and t_max to t_initial to get the time window
        //return the time window for that module? 
        //all these module triggers can then be added to a map of neighbors/loosely clustered events in the main part of the script.
            //need to create event clusters in the main function as well and then just add to them
    }


    

    const std::vector<ModuleTrigger>& l1_lightConeReconstructedEvents;{

        auto modulePosition=[&](ModuleKey m){
		auto geoIt=geometry->omgeo.find(OMKey(m.GetString(),m.GetOM(),1));
        // write something to get the time of the m
		if(geoIt==geometry->omgeo.end())
			log_fatal_stream("Module " << m << " not in geometry?");
// 			throw std::runtime_error("Module not in geometry?");
		return geoIt->second.position;
	};
        
    } return ModuleKey.push_back(geoIt);

}

// Function to create event clusters within a time window and remove duplicates based on moduleID
std::vector<std::vector<ModuleTrigger>> createEventClusters(const std::vector<ModuleTrigger>& l1_lightConeReconstructedEvents, double timeWindow) {
    std::vector<std::vector<ModuleTrigger>> l1_skimmedEventClusters;

    // Iterate through sorted events to create clusters
    for (const auto& event : l1_lightConeReconstructedEvents) {
        bool addedToExistingCluster = false;

        // Check if the event can be added to an existing cluster within the time window
        for (auto& cluster : l1_skimmedEventClusters) {
            if (std::abs(event.time - cluster.back().time) <= timeWindow) {
                // Check for duplicate based on ModuleKey
                auto duplicates = std::find_if(cluster.begin(), cluster.end(),
                    [&event](const ModuleTrigger& clusterEvent) {
                        return event.module == clusterEvent.module;
                    });

                if (duplicates == cluster.end()) {
                    cluster.push_back(event);
                }

                addedToExistingCluster = true;
                break;
            }
        }


        // If the event couldn't be added to an existing cluster, create a new cluster
        if (!addedToExistingCluster) {
            l1_skimmedEventClusters.push_back({event});
        }
    }

    return l1_skimmedEventClusters;
}

//need to add lines to read in I3 files 

int main(int argc, char* argv[]){
    // timewindow for streamline module waveform reports
    double timeWindow = 1600; //ns
    // time order the l0 triggers... for now this assumes trigger fixed length :/ yikes... a later issue... 
    timeOrderedL0Triggers(const std::vector<ModuleTrigger>& l0_triggers);
    // perform the light cone algorithm 
    for(i = 0, i < = l0_triggers.length, i ++) {
        //create an event cluster for each initial trigger 
        for (int j = i + 1; j <= l0_triggers.length; j ++){
            std::vector<ModuleTrigger> pulse = l0_triggers[j];
            std::vector<ModuleTrigger> l1_reconstructedEvent = levelOneLightConeReconstruction(l0_triggers[i], pulse); 
            //add the neighbors to the event cluster 
            //also probably somewhere need to group by modules we're requesting ,,,,
        } 
    }
    
    // take outcome reconstructed event per module and request final set of module waveforms
    std::vector<ModuleTrigger> l1_requestedModules = removeDuplicateRequests(std::vector<l1_reconstructedEvent>, timeWindow);


}