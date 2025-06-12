#pragma once
#include <string>

namespace beam {

class Beam {
public:
  Beam(double length, double EI, double load, double area);
  virtual std::string get_name();
  virtual void get_centerline();
};

// class StaticBeamViewer {
// public:
//   StaticBeamViewer(Beam beam);
// private:
//   Beam *beam;
// };

}