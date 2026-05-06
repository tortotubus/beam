

/* Type definitions */

typedef enum ParamValueType {
  PARAM_VALUE_INT = 0,
  PARAM_VALUE_DOUBLE,
  PARAM_VALUE_BOOL,
  PARAM_VALUE_STRING
} ParamValueType;

typedef struct ParamValue {
  ParamValueType type;
  union {
    long i;
    double d;
    int b;   /* 0 or 1 */
    char* s; /* owned by parser/container */
  } as;
  void* ptr;
  bool unset;
} ParamValue;

typedef struct InputGroup {
  char* group_name;
  int option_count;
  char** option_names;
  ParamValue* option_values;
} InputGroup;

typedef struct InputFile {
  int input_group_count;
  InputGroup* input_groups;
} InputFile;

/* Globals */

InputFile input_file = {0};

/* Macros */

/**
 * @def input_file_register_option
 */
#define input_file_register_option(group_name, option_variable, type_name)               \
  _input_file_register_option (                                                          \
    (group_name), #option_variable, (void*) &(option_variable), (type_name))

/**
 * @def input_file_register_option_named
 * @param group_name (string) The group name
 * @param option_name (string) The option name
 * @param option_variable (int, double, bool, or string) The variable holding the option
 * value
 * @param type_name (ParamValueType) The data type of the option value
 *
 * @note The option_variable is automatically cast to a void pointer
 */
#define input_file_register_option_named(                                                \
  group_name, option_name, option_variable, type_name)                                   \
  _input_file_register_option (                                                          \
    (group_name), option_name, (void*) &(option_variable), (type_name))

/* Function Declarations */

void input_file_free ();
void _input_file_register_group (const char* group_name);
void _input_file_register_option (const char* group_name,
                                  const char* option_name,
                                  void* value_ptr,
                                  ParamValueType option_value_type);
int input_file_apply_option (const char* group_name, const char* option_name);
int input_file_apply_options ();
void input_file_print_options ();
static int _input_file_find_group (const char* group_name);
static int _input_group_find_option (const InputGroup* group, const char* option_name);
static char* _input_strdup (const char* s);

/* Function Definitions */

static char* _input_strdup (const char* s) {
  if (!s)
    return NULL;
  size_t n = strlen (s);
  char* out = (char*) malloc (n + 1);
  if (!out)
    return NULL;
  memcpy (out, s, n + 1);
  return out;
}

/**
 * @brief Free the input file
 */
void input_file_free () {
  for (int gi = 0; gi < input_file.input_group_count; gi++) {
    InputGroup* g = &input_file.input_groups[gi];
    free (g->group_name);
    g->group_name = NULL;

    for (int oi = 0; oi < g->option_count; oi++) {
      free (g->option_names[oi]);
      g->option_names[oi] = NULL;

      if (!g->option_values[oi].unset && g->option_values[oi].type == 3) {
        free (g->option_values[oi].as.s);
      }
      g->option_values[oi].as.s = NULL;
      g->option_values[oi].ptr = NULL;
      g->option_values[oi].unset = true;
    }

    free (g->option_names);
    free (g->option_values);
    g->option_names = NULL;
    g->option_values = NULL;
    g->option_count = 0;
  }

  free (input_file.input_groups);
  input_file.input_groups = NULL;
  input_file.input_group_count = 0;
}

/**
 * @brief Register a new option group name
 */
void _input_file_register_group (const char* group_name) {
  if (!group_name)
    return;
  if (_input_file_find_group (group_name) >= 0)
    return; // Already registered

  int new_count = input_file.input_group_count + 1;
  InputGroup* groups = (InputGroup*) realloc (input_file.input_groups,
                                              (size_t) new_count * sizeof (InputGroup));
  if (!groups)
    return;

  input_file.input_groups = groups;
  InputGroup* g = &input_file.input_groups[input_file.input_group_count];
  g->group_name = _input_strdup (group_name);
  if (!g->group_name)
    return;
  g->option_count = 0;
  g->option_names = NULL;
  g->option_values = NULL;
  input_file.input_group_count = new_count;
}

/**
 * @brief Register a new option
 *
 * @param group_name The name of the group where the option belongs
 * @param option_name The name of the actual option
 * @param value_ptr The pointer to the variable where the option's value should be applied
 * @param option_value_type The value type definitions
 *
 * @note The group name will be automatically created if it does not exist
 */
void _input_file_register_option (const char* group_name,
                                  const char* option_name,
                                  void* value_ptr,
                                  ParamValueType option_value_type) {
  if (!group_name || !option_name)
    return;

  int gi = _input_file_find_group (group_name);
  if (gi < 0) {
    _input_file_register_group (group_name);
    gi = _input_file_find_group (group_name);
    if (gi < 0)
      return;
  }

  InputGroup* g = &input_file.input_groups[gi];
  int oi = _input_group_find_option (g, option_name);

  if (oi >= 0) {
    if (!g->option_values[oi].unset && g->option_values[oi].type == 3) {
      free (g->option_values[oi].as.s);
      g->option_values[oi].as.s = NULL;
    }
    g->option_values[oi].type = option_value_type;
    g->option_values[oi].ptr = value_ptr;
    g->option_values[oi].unset = true;
    return;
  }

  int new_count = g->option_count + 1;
  char** option_names =
    (char**) realloc (g->option_names, (size_t) new_count * sizeof (char*));
  if (!option_names)
    return;
  g->option_names = option_names;

  ParamValue* option_values =
    (ParamValue*) realloc (g->option_values, (size_t) new_count * sizeof (ParamValue));
  if (!option_values)
    return;
  g->option_values = option_values;

  g->option_names[g->option_count] = _input_strdup (option_name);
  if (!g->option_names[g->option_count])
    return;

  ParamValue* v = &g->option_values[g->option_count];
  memset (v, 0, sizeof (*v));
  v->type = option_value_type;
  v->ptr = value_ptr;
  v->unset = true;

  g->option_count = new_count;
}

/**
 * @brief Apply a single registered option
 *
 * @param group_name The name of the group where the option belongs
 * @param option_name The name of the actual option
 *
 * @returns 0 if successful, -1 on failure
 */
int input_file_apply_option (const char* group_name, const char* option_name) {
  int gi = _input_file_find_group (group_name);
  if (gi < 0)
    return -1;

  InputGroup* g = &input_file.input_groups[gi];
  int oi = _input_group_find_option (g, option_name);
  if (oi < 0)
    return -1;

  ParamValue* v = &g->option_values[oi];
  if (v->unset || !v->ptr)
    return -1;

  switch (v->type) {
  case 0:
    *((int*) v->ptr) = (int) v->as.i;
    return 0;
  case 1:
    *((double*) v->ptr) = v->as.d;
    return 0;
  case 2:
    *((bool*) v->ptr) = v->as.b != 0;
    return 0;
  case 3:
    *((char**) v->ptr) = v->as.s;
    return 0;
  default:
    return -1;
  }
}

/**
 * @brief Apply all registered options
 *
 * @returns 0 if successful, -1 on failure
 */
int input_file_apply_options () {
  for (int gi = 0; gi < input_file.input_group_count; gi++) {
    const char* group_name = input_file.input_groups[gi].group_name;
    InputGroup* g = &input_file.input_groups[gi];
    for (int oi = 0; oi < g->option_count; oi++) {
      if (input_file_apply_option (group_name, g->option_names[oi]) != 0)
        return -1;
    }
  }
  return 0;
}

static int _input_file_find_group (const char* group_name) {
  if (!group_name)
    return -1;
  for (int gi = 0; gi < input_file.input_group_count; gi++) {
    if (input_file.input_groups[gi].group_name &&
        strcmp (input_file.input_groups[gi].group_name, group_name) == 0)
      return gi;
  }
  return -1;
}

static int _input_group_find_option (const InputGroup* group, const char* option_name) {
  if (!group || !option_name)
    return -1;
  for (int oi = 0; oi < group->option_count; oi++) {
    if (group->option_names[oi] && strcmp (group->option_names[oi], option_name) == 0)
      return oi;
  }
  return -1;
}

/**
 * @brief Prints proccessed options
 */
void input_file_print_options () {
  if (pid () == 0) {
    printf ("Processed options:\n");
    for (int gi = 0; gi < input_file.input_group_count; gi++) {
      InputGroup* g = &input_file.input_groups[gi];
      printf ("[%s]\n", g->group_name ? g->group_name : "(null)");

      for (int oi = 0; oi < g->option_count; oi++) {
        const char* option_name = g->option_names[oi] ? g->option_names[oi] : "(null)";
        ParamValue* v = &g->option_values[oi];

        printf ("  %s = ", option_name);
        if (v->unset) {
          if (!v->ptr) {
            printf ("<unset>\n");
            continue;
          }

          switch (v->type) {
          case 0:
            printf ("%d <default>\n", *((int*) v->ptr));
            break;
          case 1:
            printf ("%g <default>\n", *((double*) v->ptr));
            break;
          case 2:
            printf ("%s <default>\n", *((bool*) v->ptr) ? "true" : "false");
            break;
          case 3: {
            char* s = *((char**) v->ptr);
            printf ("%s <default>\n", s ? s : "(null)");
            break;
          }
          default:
            printf ("<unknown-type>\n");
            break;
          }
          continue;
        }

        switch (v->type) {
        case 0:
          printf ("%ld\n", v->as.i);
          break;
        case 1:
          printf ("%g\n", v->as.d);
          break;
        case 2:
          printf ("%s\n", v->as.b ? "true" : "false");
          break;
        case 3:
          printf ("%s\n", v->as.s ? v->as.s : "(null)");
          break;
        default:
          printf ("<unknown-type>\n");
          break;
        }
      }
    }
  }
}
