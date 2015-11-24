  * _12 November 2013._ Bug in the constrained tree option that prevented some constrained tree files to be read properly. Affects development and stable versions. For fixing this issue in the stable version, please re-install version 20120412, and copy the patch file 20131112.patch in the src/ directory. Type 'patch -i 20131112.patch' to apply the changes. You will then need to compile the code as per usual (see manual if required).

  * _16 October 2013._ Minor bug in bootstrap function prevent the analysis to proceed normally under GTR (and probably Custom model for nucleotide data). Affects all development versions. Please upgrade to 20131016 version.

  * _05 October 2013._ Serious bug in bootstrap function prevent the analysis to proceed normally. Affects all development versions. Please upgrade to 20131005 version.

  * _27 September 2013._ Bugs in bootstrap functions prevent the analysis to proceed normally when using GTR or custom models of substitution. Affects development 20130708 version. Please upgrade to 20130927 version.

  * _13 May 2013._ Bugs in bootstrap functions prevent the analysis to proceed normally. Only affect the development version 20130412. Please use version 20130513 instead.

  * _19 Fevrier 2013._ The part of the code dealing with reading data in NEXUS format and processing discrete characters was broken. A new patch (20130219.patch) provides a quick fix to this. It has not been tested extensively and might need refinement. Please do not hesitate to post feedback on the forum (https://groups.google.com/forum/?fromgroups#!forum/phyml-forum). To apply this patch, copy the file 20130219.patch in the src/ directory and enter the command 'patch -i 20130219.patch'. You will need to re-compile the program after that.

  * _18 Fevrier 2013._ There was an error in the last patch (20130211.patch). Please replace it with 20130218.patch and follow the instructions below to apply it.

  * _11 Fevrier 2013._ The optimisation of the FreeRate model parameters sometimes converged to sub-optimal solutions. Adding multi-dimension optimisation routines (instead of sequential optimisation) seems to be helping. To apply the corresponding patch, copy the file 20130211.patch in the src/ directory and type 'patch -i 20130211.patch'. Note that this patch, like all the others listed here, incorporates all the bug fixes mentioned below.

  * _08 December 2012._  When analysing multiple input trees and using one of the fast approximate edge support methods, edge lengths in the output ML tree are incorrect. Thanks again to Lars Jermiin for pointing this out. To apply the corresponding patch, copy the file 20121208.patch in the src/ directory and type 'patch -i 20121208.patch'.

  * _07 December 2012._ Changed maximum length of a branch from 100 to 10. Many thanks to Lars Jermiin for running thorough tests on this. To apply the corresponding patch, copy the file 20121207.patch in the src/ directory and type 'patch -i 20121207.patch'.

  * _08 November 2012._ Minor bug in version 20120412 in the printing of edge lengths after calculating branch length supports using aLRT, aBayes or SH. To apply the patch, copy the file 20121108.patch in the src/ directory and type 'patch -i 20121108.patch'.


  * _02 November 2012._ Minor bug in version 20120412 that prevented PhyML to print out the statistics file properly when processing multiple data sets. To apply the patch, copy the file 20121102.patch in the src/ directory and type 'patch -i 20121102.patch'.


  * _11 July 2012._ A bug in version 20120412 prevented PhyTime to read in calibration information at tip nodes properly. This bug only affects serial sample analyses. To apply the patch, copy the file 20120711.patch in the src/ directory and type 'patch -i 20120711.patch'. Note that this patch incorporate the modifications from the previous patch (20120312.patch), making this last one obsolete.


  * _15 May 2012._ A bug in version 20120412 prevented constrained tree analysis to proceed. This patch fixes it. It also merges with the previous patch for calculating likelihood scores accurately on very big tree (see just below). To apply it, copy the corresponding file in the src/ directory and type 'patch -i 20120312.patch'.


  * _27 April 2012._ With very large data sets (usually >1,000 taxa), numerical precision issues sometimes cause the likelihood to be 0. The
program then stops and prints out a message that looks like the one below:
```
. Site = 267
. invar = -1
. scale_left = 0 scale_rght = 0
. Lk = 0 LOG(Lk) = -inf < -1.79769E+308
. rr=0.136954 p=0.250000
. rr=0.476752 p=0.250000
. rr=1.000000 p=0.250000
. rr=2.386294 p=0.250000
. pinv = 0.2
. Err in file lk.c at line 739
```
If this occurs with your data set, please use the patch biglk.patch available from the Downloads page on this site. To apply this patch, just go in the src/ directory and type 'patch -i biglk.patch'. You will then need to compile the program again. Use this patch with version 20120412 of PhyML only.

  * _19 September 2011._ Versions prior to phyml-20110919.tar.gz have an important bug affecting PhyTime. **Priors on node ages are not recorded properly** (except for the root node). Analysis conducted with versions of PhyML older that this one need to be re-run. Sorry for the inconvenience.

  * _04 March 2011._ Due to a bug introduced in [revision 568](https://code.google.com/p/phyml/source/detail?r=568) (commited the 3rd of Nov. 2010), the bootstrap analysis would not complete. It has been fixed in [revision 596](https://code.google.com/p/phyml/source/detail?r=596).

  * _03 Feb 2010._ Bug found in **bootstrap calculation**. The comparison of bipartition of equal size could fail when the first character of taxon names were identical. I expect the impact of this bug to be minor.

  * _03 Feb 2010._ The calculation of **unconstrained likelihood** values ignore gap and indels. It is not a bug but a more sophisticated treatment of those characters would be more appropriate.

  * _14 Sep 2010._ The calculation of **bootstrap values with the MPI version of PhyML development version** is flawed. This bug is expected to have been introduced in the program early in 2010. It does not affect the non-development version of PhyML nor the binaries available on http://www.atgc-montpellier.fr.