#include "tests/basilisk_tests.h"
#include "library/ibm/uhllman/centered.h"

int main(void) {
  mu_test_case tests[] = {
    // {"test_addition", test_addition},
    // {"test_division", test_division},
  };

  return mu_run_all(tests, (int)(sizeof(tests) / sizeof(tests[0])));
}