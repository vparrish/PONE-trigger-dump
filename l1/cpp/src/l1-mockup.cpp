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
void timeOrderedL0Triggers(const std::vector<ModuleTrigger>& levelZeroTriggers){
    std::vector<ModuleTrigger> orderedLevelZeroTriggers = levelZeroTriggers;
    std::sort(orderedLevelZeroTriggers.begin(), orderedLevelZeroTriggers.end(), compareTriggerTimes);
}

// perform geometry and return vector of modules to retrieve 
double std::vector<ModuleTrigger> levelOneEventReconstruction(){
    std::vector<ModuleTrigger> levelOneReconstructedEvents;{
        
    } return levelOneReconstructedEvents;

}



// Function to create event clusters within a time window and remove duplicates
std::vector<std::vector<ModuleTrigger>> removeDuplicateRequests(const std::vector<ModuleTrigger>& levelOneReconstructedEvents, double timeWindow) {
    std::vector<std::vector<ModuleTrigger>> eventClusters;

    // Iterate through sorted events to create clusters
    for (const auto& event : sortedEvents) {
        bool addedToExistingCluster = false;

        // Check if the event can be added to an existing cluster within the time window
        for (auto& cluster : eventClusters) {
            if (std::abs(event.time - cluster.back().time) <= timeWindow) {
                // Check for duplicate based on moduleID
                auto duplicateIt = std::find_if(cluster.begin(), cluster.end(),
                    [&event](const ModuleTrigger& clusterEvent) {
                        return event.module == clusterEvent.module;
                    });

                if (duplicateIt == cluster.end()) {
                    cluster.push_back(event);
                }

                addedToExistingCluster = true;
                break;
            }
        }

        // If the event couldn't be added to an existing cluster, create a new cluster
        if (!addedToExistingCluster) {
            eventClusters.push_back({event});
        }
    }

    return eventClusters;
}

int main(int argc, char* argv[]){
    // timewindow for streamline module waveform reports
    double timeWindow = 1600; //ns

    // time order the l0 triggers... for now this assumes trigger fixed length :/ yikes... a later issue... 
    timeOrderedL0Triggers(const std::vector<ModuleTrigger>& levelZeroTriggers);
    // perform the light cone algorithm 
    std::vector<ModuleTrigger> reconstructedEvent = levelOneEventReconstruction();
    // take outcome reconstructed event per module and request final set of module waveforms
    std::vector<ModuleTrigger> requestedModules = removeDuplicateRequests(std::vector<reconstructedEvent>, timeWindow);


}

