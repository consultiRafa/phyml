/*
 * author: Imran Fanaswala
 */

#ifndef BEAGLE_UTILS_H
#define BEAGLE_UTILS_H

#include "assert.h"
#include "libhmsbeagle-1/libhmsbeagle/beagle.h"
#include "utilities.h"

void print_beagle_resource_list();
void print_beagle_instance_details(BeagleInstanceDetails* inst);
int  create_beagle_instance(t_tree* tree, int quiet);
int  finalize_beagle_instance(t_tree* tree);
void update_beagle_partials(t_tree* tree, t_edge* b, t_node* d);
#endif // BEAGLE_UTILS_H
