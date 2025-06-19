#include "../interface/basilisk_interface.h"

int main() {
  int base = 1;
  int sq = 2;
  int result = 0;

  beam_handle_t *dummy = beam_dummy_create(base);
  result = beam_dummy_compute(dummy, sq);
  beam_dummy_destroy(dummy);

  printf("Result: %d\n", result);

  return 0;
}