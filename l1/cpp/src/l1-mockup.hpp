// l1-mockup.h
// Forward declaration of ModuleKey
struct ModuleKey;


struct ModuleTrigger {
    ModuleKey module;
    unsigned int multiplicity;
    double time;
};

struct NeighborModule {
    ModuleKey module;
    auto time;
}

#endif // L1_MOCKUP_H
