/* Generic event layer. Include after selecting the Basilisk grid.
 * A driver populates BasiliskXMLCase, constructs the ELFF runtime, then calls
 * generic_ibm_case_run(). This template deliberately has no XML-version-specific
 * main(): the complete v2 reader/model-construction migration is separate.
 */ 
#include "library/io/xml/xml.h"

int main() {
    generic_ibm_case_run();
    return 0;
}