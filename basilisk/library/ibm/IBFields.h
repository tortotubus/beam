// ============================================================================
// Type Definitions
// ============================================================================

/**
 * @struct IBscalar
 */
typedef struct {
  int i;
} IBscalar;

/**
 * @struct IBvector
 */
typedef struct {
#if dimension == 1
  IBscalar x;
#elif dimension == 2
  IBscalar x;
  IBscalar y;
#else // dimension == 3
  IBscalar x;
  IBscalar y;
  IBscalar z;
#endif
} IBvector;

/**
 * @struct IBAttributes
 */
typedef struct {
  char* name;
  IBvector v;
  bool nodump;
  bool dirty;
} IBAttributes;

// ============================================================================
// Global Declarations
// ============================================================================

IBAttributes* _ibattribute = NULL;
IBscalar* iball = NULL;
size_t _ibattribute_len = 0;
size_t nibvar = 0;

IBvector npos;
IBvector nvel;
IBvector nforce;

// ============================================================================
// Function Declarations
// ============================================================================

void init_ibsolver ();
inline double* _ibval (IBscalar s, IBNode* n);

IBscalar _init_ibscalar (const char* name);

int iblist_len (IBscalar* list);
IBscalar* iblist_append (IBscalar* list, IBscalar sc);
IBscalar* iblist_prepend (IBscalar* list, IBscalar sc);
IBscalar* iblist_add (IBscalar* list, IBscalar sc);
int iblist_lookup (IBscalar* l, IBscalar s1);
IBscalar* iblist_copy (IBscalar* l);
IBscalar* iblist_concat (IBscalar* l1, IBscalar* l2);
void iblist_print (IBscalar* l, FILE* fp);

IBvector _init_ibvector (const char* name);

int ibvectors_len (IBvector* list);
IBvector* ibvectors_append (IBvector* list, IBvector v);
IBvector* ibvectors_add (IBvector* list, IBvector vv);
IBvector* ibvectors_copy (IBvector* l);

// ============================================================================
// Macros
// ============================================================================

/**
 * @def
 */
// clang-format off
@define new_ibscalar(name) name = _init_ibscalar (#name);
// clang-format on

/**
 * @def
 */
// clang-format off
@define new_ibvector(name) name = _init_ibvector (#name);
// clang-format on

/**
 * @def
 * @note Assume node exists in scope and is a pointer of IBNode
 */
// clang-format off
#define ibval(s) (*_ibval((s), node))
// clang-format on

/**
 * @def
 * @note Assume node exists in scope and is a pointer of IBNode
 */
// clang-format off
#define ibname(s) (char*)(_ibattribute[s.i].name)
// clang-format on

/**
 * @def
 * @note Assume node exists in scope and is a pointer of IBNode
 */
// clang-format off
#define ibdirty(s) (_ibattribute[(s).i].dirty)
// clang-format on

/**
 * @def
 * @note Assume node exists in scope and is a pointer of IBNode
 */
// clang-format off
#define ibnodump(s) (_ibattribute[(s).i].nodump)
// clang-format on

/**
 * @def
 * @note
 */
macro foreach_ibscalar (IBscalar* list = iball) {
  {
    IBscalar* _i = (IBscalar*) (list);
    if (_i)
      for (IBscalar s = *_i; (&s)->i >= 0; s = *++_i) {
        // clang-format off
        {...}
        // clang-format on
      }
  }
}

/**
 * @def
 * @note
 */
macro foreach_ibvector (IBvector* list) {
  {
    IBvector* _i = (IBvector*) (list);
    if (_i)
      for (IBvector v = *_i; (&v)->x.i >= 0; v = *++_i) {
        // clang-format off
        {...}
        // clang-format on
      }
  }
}

// ============================================================================
// Function definitions
// ============================================================================

/**
 * @brief
 */
void init_ibsolver () {
  new_ibvector (nvel);
  new_ibvector (nforce);
  new_ibvector (npos);

  int n = nibvar;
  if (iball)
    free (iball);
  iball = (IBscalar*) malloc (sizeof (IBscalar) * (n + 1));
  for (int i = 0; i < n; i++)
    iball[i].i = i;
  iball[n].i = -1;
}

// ============================================================================
// IBScalar
// ============================================================================

/**
 * @brief
 */
IBscalar _init_ibscalar (const char* name) {
  IBscalar s = {.i = nibvar++};

  if ((size_t) (s.i + 1) > _ibattribute_len) {
    size_t old_len = _ibattribute_len;
    _ibattribute_len = (size_t) (s.i + 1);
    if (_ibattribute == NULL)
      _ibattribute =
        (IBAttributes*) calloc (_ibattribute_len, sizeof (IBAttributes));
    else
      _ibattribute = (IBAttributes*) realloc (
        _ibattribute, _ibattribute_len * sizeof (IBAttributes));
    assert (_ibattribute);
    for (size_t i = old_len; i < _ibattribute_len; i++)
      _ibattribute[i] = (IBAttributes) {0};
  }

  if (_ibattribute[s.i].name)
    free (_ibattribute[s.i].name);
  _ibattribute[s.i].name = strdup (name ? name : "");

  foreach_dimension() {
    _ibattribute[s.i].v.x.i = -1;
  }

  return s;
}

/**
 * @brief
 */
inline double* _ibval (IBscalar s, IBNode* n) {
  // if (set_dirty)
  //   _ibattribute[s.i].dirty = true;

  return &((double*) ((char*) n + sizeof (IBNode)))[s.i];
}

/**
 * @brief
 */
int iblist_len (IBscalar* list) {
  if (!list)
    return 0;
  int ns = 0;
  foreach_ibscalar (list) ns++;
  return ns;
}

/**
 * @brief
 */
IBscalar* iblist_append (IBscalar* list, IBscalar sc) {
  int len = iblist_len (list);
  qrealloc (list, len + 2, IBscalar);
  list[len] = sc;
  list[len + 1].i = -1;
  return list;
}

/**
 * @brief
 */
IBscalar* iblist_prepend (IBscalar* list, IBscalar sc) {
  int len = iblist_len (list);
  qrealloc (list, len + 2, IBscalar);
  for (int i = len; i >= 1; i--)
    list[i] = list[i - 1];
  list[0] = sc;
  list[len + 1].i = -1;
  return list;
}

/**
 * @brief
 */
IBscalar* iblist_add (IBscalar* list, IBscalar sc) {
  foreach_ibscalar (list) {
    if (s.i == sc.i)
      return list;
  }
  return iblist_append (list, sc);
}

/**
 * @brief
 */
int iblist_lookup (IBscalar* l, IBscalar s1) {
  if (l != NULL)
    foreach_ibscalar (l) if (s1.i == s.i) return true;
  return false;
}

/**
 * @brief
 */
IBscalar* iblist_copy (IBscalar* l) {
  IBscalar* list = NULL;
  if (l != NULL)
    foreach_ibscalar (l) list = iblist_append (list, s);
  return list;
}

/**
 * @brief
 */
IBscalar* iblist_concat (IBscalar* l1, IBscalar* l2) {
  IBscalar* l3 = iblist_copy (l1);
  foreach_ibscalar (l2) l3 = iblist_append (l3, s);
  return l3;
}

/**
 * @brief
 */
void iblist_print (IBscalar* l, FILE* fp) {
  int i = 0;
  foreach_ibscalar (l) fprintf (fp, "%s%s", i++ == 0 ? "{" : ",", ibname (s));
  fputs (i > 0 ? "}\n" : "{}\n", fp);
}

// ============================================================================
// IBvector
// ============================================================================

/**
 * @brief
 */
IBvector _init_ibvector (const char* name) {
  struct {
    char *x, *y, *z;
  } ext = {".x", ".y", ".z"};

  IBvector v = {0};

  foreach_dimension () {
    if (name) {
      char cname[strlen (name) + 3];
      strcat (strcpy (cname, name), ext.x);
      v.x = _init_ibscalar (cname);
    } else {
      v.x = _init_ibscalar (NULL);
    }
  }

  foreach_dimension () {
    _ibattribute[v.x.i].v = v;
  }

  return v;
}

/**
 * @brief
 */
int ibvectors_len (IBvector* list) {
  if (!list)
    return 0;
  int nv = 0;
  foreach_ibvector (list) nv++;
  return nv;
}

/**
 * @brief
 */
IBvector* ibvectors_append (IBvector* list, IBvector v) {
  int len = ibvectors_len (list);
  qrealloc (list, len + 2, IBvector);
  list[len] = v;
  list[len + 1] = (IBvector) {{-1}};
  return list;
}

/**
 * @brief
 */
IBvector* ibvectors_add (IBvector* list, IBvector vv) {
  foreach_ibvector (list) {
    bool id = true;
    foreach_dimension () {
      if (v.x.i != vv.x.i)
        id = false;
      if (id)
        return list;
    }
  }
  return ibvectors_append (list, vv);
}

/**
 * @brief
 */
IBvector* ibvectors_copy (IBvector* l) {
  IBvector* list = NULL;
  if (l != NULL)
    foreach_ibvector (l) list = ibvectors_append (list, v);

  return list;
}
