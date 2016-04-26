#include "particles.h"

using namespace CGL::StaticScene;

namespace CGL {

    bool Particles::simulateToTime(double t) {
        if (t <= pt->fluid_simulate_time) {
            return false;
        }
        return true;
    }

}  // namespace CGL
