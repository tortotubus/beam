#include "../config/config.hpp"

namespace beam {
  namespace fem {
  class Coefficient {
    protected:
      real_t time;
    public:
      Coefficient() { time = 0.; }
      virtual void SetTime(real_t t) { time = t; }
      virtual real_t Eval() = 0;
  };

  class ConstantCoefficient : public Coefficient {
    public:
      real_t constant;
      explicit ConstantCoefficient(real_t c = 1.0) { constant = c; }
      real_t Eval() override { return (constant); }
  };
}
}