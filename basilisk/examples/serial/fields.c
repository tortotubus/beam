
// /**
//  * @struct IBscalar
//  */
// typedef struct
// {
//   int i;
// } IBscalar;

// /**
//  * @struct IBvector
//  */
// typedef struct
// {
// #if dimension == 1
//   IBscalar x;
// #elif dimension == 2
//   IBscalar x;
//   IBscalar y;
// #else // dimension == 3
//   IBscalar x;
//   IBscalar y;
//   IBscalar z;
// #endif
// } IBvector;

// /**
//  * @struct IBAttributes
//  */
// typedef struct
// {
//   const char* name;
// } IBAttributes;

// /* Globals */
// IBAttributes _ibattributes = NULL;
// int _ibattributes_count = 0;

// /* Function decalarations */
// IBscalar
// _init_ibscalar(IBscalar s, const char* name);
// IBvector
// _init_ibvector(IBvector v, const char* name);

// /* Macros */

// /**
//  *
//  */
// #define new_ibscalar(name)                                                     \
//   IBscalar name =                                                              \
//     _init_ibscalar((IBscalar){ .i = _ibattributes_count + 1 }, #name);

// #if dimension == 1
// #define new_ibvector(name)                                                     \
//   IBvector name = _init_ibvector(                                              \
//     (IBvector){                                                                \
//       .x = { .i = _ibattributes_count + 1 },                                   \
//     },                                                                         \
//     #name);
// #elif dimension == 2
// #define new_ibvector(name)                                                     \
//   IBvector name = _init_ibvector(                                              \
//     (IBvector){                                                                \
//       .x = { .i = _ibattributes_count + 1 },                                   \
//       .y = { .i = _ibattributes_count + 2 },                                   \
//     },                                                                         \
//     #name);
// #else
// #define new_ibvector(name)                                                     \
//   IBvector name = _init_ibvector(                                              \
//     (IBvector){                                                                \
//       .x = { .i = _ibattributes_count + 1 },                                   \
//       .y = { .i = _ibattributes_count + 2 },                                   \
//       .z = { .i = _ibattributes_count + 3 },                                   \
//     },                                                                         \
//     #name);
// #endif

// /* Function definitions */

// /**
//  * @brief
//  *
//  */
// IBscalar
// _init_ibscalar(IBscalar s, const char* name)
// {
//   (void)name;
//   return s;
// }
