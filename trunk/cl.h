/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#ifndef CL_H
#define CL_H

#include <unistd.h>
#include <getopt.h>
#include "utilities.h"
#include "options.h"
#include "models.h"
#include "free.h"
#include "interface.h"


void Read_Command_Line(option *input, int argc, char **argv);

#endif
