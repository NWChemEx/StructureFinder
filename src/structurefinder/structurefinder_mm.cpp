#include "structurefinder/structurefinder_mm.hpp"

namespace structurefinder {

inline void set_defaults(pluginplay::ModuleManager& mm) {
    // Default submodules between collections can be set here
}

DECLARE_PLUGIN(structurefinder) {
    // Add subcollection load calls here

    // Assign default submodules
    set_defaults(mm);
}

} // namespace structurefinder
