#include "library/io/params/params-txt.h"
 

typedef struct InputFileCLIOptions {
  const char* program_name;
  const char* input_path;
  int show_help;
  int had_error;
} InputFileCLIOptions;

static inline void
input_file_cli_options_init (InputFileCLIOptions* opts) {
  if (!opts)
    return;
  opts->program_name = "program";
  opts->input_path = NULL;
  opts->show_help = 0;
  opts->had_error = 0;
}

static inline void
input_file_print_help (const char* program_name, FILE* out) {
  if (!out)
    return;
  if (!program_name || !program_name[0])
    program_name = "program";

  fprintf (out, "Usage: %s [--help|-h] [input-file]\n", program_name);
  fprintf (out, "\n");
  fprintf (out, "Arguments:\n");
  fprintf (out, "  input-file    Optional params file in [group] key = value format.\n");
  fprintf (out, "\n");
  fprintf (out, "Options:\n");
  fprintf (out, "  -h, --help    Show this help and exit.\n");
}

/**
 * @brief Parse CLI arguments into a small options struct.
 *
 * Recognized options:
 * - `-h`, `--help`
 *
 * Positional handling:
 * - The first positional argument is treated as the params input file.
 * - Additional positional arguments are currently rejected as errors.
 *
 * @return 0 on success, -1 on parse error.
 *
 * @note This assumes regular `main(argc, argv)` style and parses user args
 *       from `argv[1..argc-1]`.
 */
static inline int
input_file_parse_cli_args (int argc, char** argv, InputFileCLIOptions* opts) {
  if (!opts) {
    fprintf (stderr, "params-cli: options pointer is NULL\n");
    return -1;
  }
  input_file_cli_options_init (opts);

  if (!argv || argc < 0) {
    fprintf (stderr, "params-cli: invalid argv/argc\n");
    opts->had_error = 1;
    return -1;
  }

  if (argc > 0 && argv[0] && argv[0][0] != '\0')
    opts->program_name = argv[0];

  int start = argc > 0 ? 1 : 0;
  for (int i = start; i < argc; i++) {
    const char* arg = argv[i];
    if (!arg || !arg[0])
      continue;

    if (strcmp (arg, "--help") == 0 || strcmp (arg, "-h") == 0) {
      opts->show_help = 1;
      continue;
    }

    if (arg[0] == '-') {
      fprintf (stderr, "params-cli: unrecognized option '%s'\n", arg);
      opts->had_error = 1;
      return -1;
    }

    if (!opts->input_path) {
      opts->input_path = arg;
    } else {
      fprintf (stderr, "params-cli: unexpected extra positional argument '%s'\n", arg);
      opts->had_error = 1;
      return -1;
    }
  }

  return 0;
}

/**
 * @brief Convenience wrapper:
 *        parse CLI, handle --help, and parse params file if provided.
 *
 * @return 0 on success/no-op, 1 if help printed, -1 on error.
 */
static inline int
input_file_parse_cli (int argc, char** argv) {
  InputFileCLIOptions opts;
  if (input_file_parse_cli_args (argc, argv, &opts) != 0)
    return -1;

  if (opts.show_help) {
    input_file_print_help (opts.program_name, stdout);
    return 1;
  }

  if (opts.input_path) {
    if (input_file_parse (opts.input_path) != 0) {
      fprintf (stderr, "params-cli: failed to parse '%s'\n", opts.input_path);
      return -1;
    }
  }

  return 0;
}
