#include "IBNodeList.h"

macro foreach_ibnodelist (IBNodeList* list, bool local_only = true) {
  for (size_t _node_index = 0; _node_index < (list)->size; _node_index++) {
    IBNode* node = &(list)->nodes[_node_index];
    if (local_only) {
      if (node->pid == pid ()) {
        // clang-format off
        {...}
        // clang-format on
      }
    } else {
      // clang-format off
      {...}
      // clang-format on      
    }
  }

}
