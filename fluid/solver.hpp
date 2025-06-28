#include "config/config.hpp"
namespace fluid {
  class FluidModel {
    public:
      virtual void step(real_t force_field, real_t dt);
      virtual void get_velocity();
      virtual void get_pressure();
      virtual void get_vorticity();
      virtual void get_grid_info();
  };

  class FluidModel : 
}

