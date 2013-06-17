#ifndef BEAGLE_UTILS_H
#define BEAGLE_UTILS_H

#include "libhmsbeagle-1/libhmsbeagle/beagle.h"
#include "utilities.h"


/* This counts the number of args */
#define NARGS_SEQ(_1,_2,_3,_4,_5,_6,_7,_8,N,...) N
#define NARGS(...) NARGS_SEQ(__VA_ARGS__, 8, 7, 6, 5, 4, 3, 2, 1)

/* This will let macros expand before concating them */
#define PRIMITIVE_CAT(x, y) x ## y
#define CAT(x, y) PRIMITIVE_CAT(x, y)

/* This will call a macro on each argument passed in */
#define FOR_EACH(macro, ...) CAT(FOR_EACH_, NARGS(__VA_ARGS__))(macro, __VA_ARGS__)
#define FOR_EACH_1(m, x1) m(x1)
#define FOR_EACH_2(m, x1, x2) m(x1) m(x2)
#define FOR_EACH_3(m, x1, x2, x3) m(x1) m(x2) m(x3)
#define FOR_EACH_4(m, x1, x2, x3, x4) m(x1) m(x2) m(x3) m(x4)
#define FOR_EACH_5(m, x1, x2, x3, x4, x5) m(x1) m(x2) m(x3) m(x4) m(x5)
#define FOR_EACH_6(m, x1, x2, x3, x4, x5, x6) m(x1) m(x2) m(x3) m(x4) m(x5) m(x6)
#define FOR_EACH_7(m, x1, x2, x3, x4, x5, x6, x7) m(x1) m(x2) m(x3) m(x4) m(x5) m(x6) m(x7)
#define FOR_EACH_8(m, x1, x2, x3, x4, x5, x6, x7, x8) m(x1) m(x2) m(x3) m(x4) m(x5) m(x6) m(x7) m(x8)
#define FOR_EACH_9(m, x1, x2, x3, x4, x5, x6, x7, x8, x9) m(x1) m(x2) m(x3) m(x4) m(x5) m(x6) m(x7) m(x8) m(x9)

#define DUMP_EACH_INT(v)    fprintf(stderr,"\n\t\tDEBUG:%s:\t\t%s--->%i\n",__PRETTY_FUNCTION__,#v,(v));
#define DUMP_EACH_STRING(v) fprintf(stderr,"\n\t\tDEBUG:%s:\t\t%s--->%s\n",__PRETTY_FUNCTION__,#v,(v));
#define DUMP_I(...) FOR_EACH(DUMP_EACH_INT, __VA_ARGS__)
#define DUMP_S(...) FOR_EACH(DUMP_EACH_STRING, __VA_ARGS__)

void print_beagle_resource_list();
void print_beagle_instance_details(BeagleInstanceDetails* inst);
int  create_beagle_instance(t_tree* tree, int quiet);
int  finalize_beagle_instance(t_tree* tree);
#endif // BEAGLE_UTILS_H
