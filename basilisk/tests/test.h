// tinytest.h
#pragma once
// #include <stdio.h>
// #include <stdlib.h>
// #include <stdarg.h>
// #include <string.h>

typedef void (*tt_test_fn)(void);

typedef struct {
  int tests_run;
  int tests_failed;
  int asserts_run;
  int asserts_failed;
  const char *current_test;
} tt_state_t;

static tt_state_t tt_g;

static void tt_vfail_(const char *file, int line,
                      const char *expr,
                      const char *fmt, va_list ap)
{
  tt_g.asserts_failed++;
  fprintf(stderr, "\n[%s] FAIL at %s:%d\n", tt_g.current_test ? tt_g.current_test : "(unknown)", file, line);
  if (expr && *expr) fprintf(stderr, "  expr: %s\n", expr);
  if (fmt && *fmt) {
    fprintf(stderr, "  ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
  }
}

static void tt_fail_(const char *file, int line, const char *expr, const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  tt_vfail_(file, line, expr, fmt, ap);
  va_end(ap);
}

static void tt_begin_(const char *name)
{
  tt_g.current_test = name;
  tt_g.tests_run++;
  tt_g.asserts_failed = 0; // per-test counter
}

static void tt_end_(void)
{
  if (tt_g.asserts_failed != 0) tt_g.tests_failed++;
}

static int tt_summary_(void)
{
  fprintf(stderr, "\n==== tinytest summary ====\n");
  fprintf(stderr, "Tests run:    %d\n", tt_g.tests_run);
  fprintf(stderr, "Tests failed: %d\n", tt_g.tests_failed);
  fprintf(stderr, "Asserts run:  %d\n", tt_g.asserts_run);
  // (Asserts failed accumulates per-test; report total-ish by counting failures in output)
  fprintf(stderr, "==========================\n");
  return tt_g.tests_failed ? 1 : 0;
}

#define TT_RUN(test_fn) do { \
  tt_begin_(#test_fn);       \
  (test_fn)();               \
  tt_end_();                 \
} while (0)

#define TT_ASSERT_TRUE(expr) do { \
  tt_g.asserts_run++;             \
  if (!(expr))                    \
    tt_fail_(__FILE__, __LINE__, #expr, "expected true"); \
} while (0)

#define TT_ASSERT_FALSE(expr) do { \
  tt_g.asserts_run++;              \
  if ((expr))                      \
    tt_fail_(__FILE__, __LINE__, #expr, "expected false"); \
} while (0)

#define TT_ASSERT_EQ_LONG(a,b) do { \
  tt_g.asserts_run++;               \
  long _a = (long)(a);              \
  long _b = (long)(b);              \
  if (_a != _b)                     \
    tt_fail_(__FILE__, __LINE__, #a " == " #b, "got %ld vs %ld", _a, _b); \
} while (0)

#define TT_ASSERT_NE_LONG(a,b) do { \
  tt_g.asserts_run++;               \
  long _a = (long)(a);              \
  long _b = (long)(b);              \
  if (_a == _b)                     \
    tt_fail_(__FILE__, __LINE__, #a " != " #b, "both were %ld", _a); \
} while (0)

#define TT_ASSERT_NONNULL(p) do { \
  tt_g.asserts_run++;             \
  const void *_p = (const void*)(p); \
  if (_p == NULL)                 \
    tt_fail_(__FILE__, __LINE__, #p, "expected non-NULL"); \
} while (0)

#define TT_ASSERT_PTR_EQ(a,b) do { \
  tt_g.asserts_run++;              \
  const void *_a = (const void*)(a); \
  const void *_b = (const void*)(b); \
  if (_a != _b)                    \
    tt_fail_(__FILE__, __LINE__, #a " == " #b, "got %p vs %p", _a, _b); \
} while (0)

#define TT_ASSERT_PTR_NE(a,b) do { \
  tt_g.asserts_run++;              \
  const void *_a = (const void*)(a); \
  const void *_b = (const void*)(b); \
  if (_a == _b)                    \
    tt_fail_(__FILE__, __LINE__, #a " != " #b, "both were %p", _a); \
} while (0)
