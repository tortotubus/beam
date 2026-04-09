#include "library/io/params/params.h"

#include <ctype.h> 
 
static char* _params_ltrim (char* s) {
  while (*s && isspace ((unsigned char) *s))
    s++;
  return s;
}

static void _params_rtrim (char* s) {
  size_t n = strlen (s);
  while (n > 0 && isspace ((unsigned char) s[n - 1])) {
    s[n - 1] = '\0';
    n--;
  }
}

static char* _params_trim (char* s) {
  s = _params_ltrim (s);
  _params_rtrim (s);
  return s;
}

static int _params_ieq (const char* a, const char* b) {
  while (*a && *b) {
    if (tolower ((unsigned char) *a) != tolower ((unsigned char) *b))
      return 0;
    a++;
    b++;
  }
  return *a == '\0' && *b == '\0';
}

/**
 * @brief Parse a plain-text input file using:
 *        [group]
 *        option = value
 *
 * @param path File path to parse
 *
 * @return 0 on success, -1 on failure (e.g. cannot open file, conversion errors)
 *
 * @note The most recent [group] header is used for subsequent options.
 * @note Unregistered options print warnings and are ignored.
 */
int input_file_parse (const char* path) {
  if (!path)
    return -1;

  FILE* fp = fopen (path, "r");
  if (!fp) {
    fprintf (stderr, "params: failed to open input file '%s'\n", path);
    return -1;
  }

  int had_errors = 0;
  int line_no = 0;
  char line[4096];
  char current_group[256] = {0};

  while (fgets (line, sizeof (line), fp)) {
    line_no++;
    char* s = _params_trim (line);

    if (*s == '\0' || *s == '#' || *s == ';')
      continue;

    if (*s == '[') {
      char* end = strchr (s, ']');
      if (!end) {
        fprintf (stderr, "params:%d: malformed group header (missing ']')\n", line_no);
        had_errors = 1;
        continue;
      }
      *end = '\0';
      char* gname = _params_trim (s + 1);
      if (*gname == '\0') {
        fprintf (stderr, "params:%d: empty group name\n", line_no);
        had_errors = 1;
        continue;
      }
      snprintf (current_group, sizeof (current_group), "%s", gname);
      continue;
    }

    if (current_group[0] == '\0') {
      fprintf (stderr,
               "params:%d: option defined before any [group], ignoring line\n",
               line_no);
      had_errors = 1;
      continue;
    }

    char* eq = strchr (s, '=');
    if (!eq) {
      fprintf (stderr, "params:%d: expected 'option = value'\n", line_no);
      had_errors = 1;
      continue;
    }

    *eq = '\0';
    char* option_name = _params_trim (s);
    char* value_str = _params_trim (eq + 1);

    if (*option_name == '\0') {
      fprintf (stderr, "params:%d: empty option name\n", line_no);
      had_errors = 1;
      continue;
    }

    int gi = _input_file_find_group (current_group);
    if (gi < 0) {
      fprintf (stderr,
               "params:%d: warning: unregistered group '%s' (option '%s' ignored)\n",
               line_no,
               current_group,
               option_name);
      continue;
    }

    InputGroup* g = &input_file.input_groups[gi];
    int oi = _input_group_find_option (g, option_name);
    if (oi < 0) {
      fprintf (stderr,
               "params:%d: warning: unregistered option '%s.%s' ignored\n",
               line_no,
               current_group,
               option_name);
      continue;
    }

    ParamValue* v = &g->option_values[oi];
    switch (v->type) {
    case 0: {
      char* endp = NULL;
      long parsed = strtol (value_str, &endp, 10);
      endp = _params_trim (endp ? endp : value_str);
      if (*value_str == '\0' || *endp != '\0') {
        fprintf (stderr,
                 "params:%d: invalid integer for '%s.%s': '%s'\n",
                 line_no,
                 current_group,
                 option_name,
                 value_str);
        had_errors = 1;
        continue;
      }
      v->as.i = parsed;
      v->unset = false;
      break;
    }
    case 1: {
      char* endp = NULL;
      double parsed = strtod (value_str, &endp);
      endp = _params_trim (endp ? endp : value_str);
      if (*value_str == '\0' || *endp != '\0') {
        fprintf (stderr,
                 "params:%d: invalid double for '%s.%s': '%s'\n",
                 line_no,
                 current_group,
                 option_name,
                 value_str);
        had_errors = 1;
        continue;
      }
      v->as.d = parsed;
      v->unset = false;
      break;
    }
    case 2: {
      int parsed = -1;
      if (_params_ieq (value_str, "1") || _params_ieq (value_str, "true") ||
          _params_ieq (value_str, "yes") || _params_ieq (value_str, "on")) {
        parsed = 1;
      } else if (_params_ieq (value_str, "0") || _params_ieq (value_str, "false") ||
                 _params_ieq (value_str, "no") || _params_ieq (value_str, "off")) {
        parsed = 0;
      }
      if (parsed < 0) {
        fprintf (stderr,
                 "params:%d: invalid bool for '%s.%s': '%s'\n",
                 line_no,
                 current_group,
                 option_name,
                 value_str);
        had_errors = 1;
        continue;
      }
      v->as.b = parsed;
      v->unset = false;
      break;
    }
    case 3: {
      char* copy_start = value_str;
      char* copy_end = value_str + strlen (value_str);
      if (copy_end > copy_start + 1 &&
          ((copy_start[0] == '"' && copy_end[-1] == '"') ||
           (copy_start[0] == '\'' && copy_end[-1] == '\''))) {
        copy_start++;
        copy_end--;
      }
      size_t n = (size_t) (copy_end - copy_start);
      char* parsed = (char*) malloc (n + 1);
      if (!parsed) {
        fprintf (stderr, "params:%d: allocation failed for string value\n", line_no);
        had_errors = 1;
        continue;
      }
      memcpy (parsed, copy_start, n);
      parsed[n] = '\0';

      if (!v->unset && v->as.s)
        free (v->as.s);
      v->as.s = parsed;
      v->unset = false;
      break;
    }
    default:
      fprintf (stderr,
               "params:%d: unsupported value type for '%s.%s'\n",
               line_no,
               current_group,
               option_name);
      had_errors = 1;
      break;
    }
  }

  fclose (fp);
  return had_errors ? -1 : 0;
}
