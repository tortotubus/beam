
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
  bool dirty;
} IBAttributes;

/* Globals */
IBAttributes* _ibattribute = NULL;
IBscalar* iball = NULL;
size_t _ibattribute_len = 0;
size_t nibvar = 0;

/* Function decalarations */
void init_ibsolver ();
IBscalar _init_ibscalar (const char* name);
IBvector _init_ibvector (const char* name);
inline double * _ibval (IBscalar s, IBNode *n, bool set_dirty);

/* Macros */

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
 *
 * @note Assume node exists in scope and is a pointer of IBNode
 */

// clang-format off
#define ibval(s) (*_ibval((s), node, ib_set_dirty))
// #define ibval(s) (((double*) ((char*) (node) + sizeof (IBNode)))[(s).i])
// clang-format on

/**
 * @def
 *
 * @note Assume node exists in scope and is a pointer of IBNode
 */
// clang-format off
#define ibname(s) (char*)(_ibattribute[s.i].name)
// clang-format on

/**
 * @def
 *
 * @note Assume node exists in scope and is a pointer of IBNode
 */
// clang-format off
// #define ibdirty(s) (bool)(_ibattribute[s.i].dirty)
#define ibdirty(s) (_ibattribute[(s).i].dirty)
// clang-format on


/**
 * @def
 *
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
 *
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

/* Function definitions */

/**
 * @brief
 */
void init_ibsolver () {
  int n = nibvar;
  iball = (IBscalar*) malloc (sizeof (IBscalar*) * (n + 1));
  for (int i = 0; i < n; i++)
    iball[i].i = i;
  iball[n].i = -1;
}

/**
 * @brief
 *
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
      _ibattribute[i].name = NULL;
  }

  if (_ibattribute[s.i].name)
    free (_ibattribute[s.i].name);
  _ibattribute[s.i].name = strdup (name ? name : "");

  return s;
}

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
    }
  }

  return v;
}


/**
 * @brief
 */
inline double * _ibval (IBscalar s, IBNode * n, bool set_dirty)
{
  if (set_dirty)
    _ibattribute[s.i].dirty = true;

  return &((double *)((char *)n + sizeof(IBNode)))[s.i];
}

/**
 * @brief
 */
int iblist_len(IBscalar * list) {
  if (!list) return 0;
  int ns = 0;
  foreach_ibscalar(list) ns++;
  return ns;
}



/**
 * @brief
 */
IBscalar * iblist_append (IBscalar * list, IBscalar sc)
{
  int len = iblist_len (list);
  qrealloc (list, len + 2, IBscalar);
  list[len] = sc;
  list[len + 1].i = -1;
  return list;
}

/**
 * @brief
 */
IBscalar * iblist_prepend(IBscalar *list, IBscalar sc) {
  int len = iblist_len (list);
  qrealloc (list, len + 2, IBscalar);
  for (int i = len; i >= 1; i--)
    list[i] = list[i-1];
  list[0] = sc;
  list[len + 1].i = -1;
  return list;
} 

/**
 * @brief
 */
IBscalar * iblist_add (IBscalar * list, IBscalar sc)
{
  foreach_ibscalar(list) {
    if (s.i == sc.i)
      return list;
  }
  return iblist_append(list,sc);
}