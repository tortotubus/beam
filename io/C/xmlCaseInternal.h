#pragma once
/* Internal libxml2 representation shared by document loading and rendering. */
#include <libxml/tree.h>
struct elff_xml_document_handle
{
  xmlDoc* document;
};
