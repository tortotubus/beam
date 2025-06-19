#include "mesh.hpp"
#include "../general/error.hpp"

namespace beam {
namespace fem {
void Mesh::Loader(std::istream &input, int generate_edges,
                  std::string parse_tag) {
  int curved = 0, read_gf = 1;
  bool finalize_topo = true;

  if (!input) {
    BEAM_ABORT("Input stream is not open");
  }
}
}
} // namespace beam
