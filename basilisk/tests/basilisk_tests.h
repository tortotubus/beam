#ifndef MINUNIT_H
#define MINUNIT_H

typedef int (*mu_test_fn)(void);

typedef struct
{
  const char* name;
  mu_test_fn fn;
} mu_test_case;

static int mu_failures = 0;

#define MU_ASSERT(expr)                                                        \
  do {                                                                         \
    if (!(expr)) {                                                             \
      fprintf(                                                                 \
        stderr, "  FAIL: %s:%d: assert(%s)\n", __FILE__, __LINE__, #expr);     \
      return 1;                                                                \
    }                                                                          \
  } while (0)

#define MU_ASSERT_EQ_INT(a, b)                                                 \
  do {                                                                         \
    int _a = (a), _b = (b);                                                    \
    if (_a != _b) {                                                            \
      fprintf(stderr,                                                          \
              "  FAIL: %s:%d: %s == %s (%d != %d)\n",                          \
              __FILE__,                                                        \
              __LINE__,                                                        \
              #a,                                                              \
              #b,                                                              \
              _a,                                                              \
              _b);                                                             \
      return 1;                                                                \
    }                                                                          \
  } while (0)

#define MU_RUN_TEST(tc)                                                        \
  do {                                                                         \
    fprintf(stderr, "RUN  %s\n", (tc).name);                                   \
    int _rc = (tc).fn();                                                       \
    if (_rc == 0) {                                                            \
      fprintf(stderr, "PASS %s\n", (tc).name);                                 \
    } else {                                                                   \
      fprintf(stderr, "FAIL %s\n", (tc).name);                                 \
      mu_failures++;                                                           \
    }                                                                          \
  } while (0)

static inline int
mu_run_all(const mu_test_case* tests, int ntests)
{
  mu_failures = 0;
  for (int i = 0; i < ntests; i++)
    MU_RUN_TEST(tests[i]);
  fprintf(stderr, "\nRESULT: %d/%d failed\n", mu_failures, ntests);
  return mu_failures ? 1 : 0; // good for CI / Makefile
}

#endif // MINUNIT_H
