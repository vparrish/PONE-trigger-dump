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

if(!geometry){
if(frame->Has("I3Geometry"))
    geometry=frame->Get<boost::shared_ptr<const I3Geometry>>("I3Geometry");
else
    return;
}

boost::shared_ptr<const I3Geometry> geometry;

double levelOneLightConeReconstruction(const std::vector<ModuleTrigger>& l0_triggers){

    std::map<ModuleKey, std::vector<EventCluster> &eventCluster> l1_skimmedEventClusters;

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

int main(int argc, char* argv[]){
    // timewindow for streamline module waveform reports
    double timeWindow = 1600; //ns
    // time order the l0 triggers... for now this assumes trigger fixed length :/ yikes... a later issue... 
    timeOrderedL0Triggers(const std::vector<ModuleTrigger>& l0_triggers);
    // perform the light cone algorithm 
    std::vector<ModuleTrigger> l1_reconstructedEvent = levelOneLightConeReconstruction();
    // take outcome reconstructed event per module and request final set of module waveforms
    std::vector<ModuleTrigger> l1_requestedModules = removeDuplicateRequests(std::vector<l1_reconstructedEvent>, timeWindow);


}

