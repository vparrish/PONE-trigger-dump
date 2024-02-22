// l1-mockup.h
// Forward declaration of ModuleKey
struct ModuleKey;


struct ModuleTrigger {
    ModuleKey module;
    unsigned int multiplicity;
    double time;
};

struct EventCluster {
    ModuleKey module;
    double time;
}

struct LCEvent {
    ModuleKey module;
    unsigned int multiplicity;
    double time;
    std::vector pos_TimeWindow;
    std::vector neg_TimeWindow;
};

#endif // L1_MOCKUP_H