PhyML is a software that estimates maximum likelihood phylogenies from alignments of nucleotide or amino acid sequences. The main strength of PhyML lies in the large number of substitution models coupled to various options to search the space  of phylogenetic tree topologies, going from very fast and efficient methods to slower but  generally more accurate approaches. PhyML was designed to process moderate to large data sets. In theory, alignments with up to 4,000 sequences 2,000,000 character-long can be processed.


---


**NEWS!**

// 28/09/2014. **The source code of PhyML has moved to GitHub and is now accessible at the following address**: https://github.com/stephaneguindon/phyml/. The sources available on the svn server provided by Google Code is no longer up to date.

// 11/09/2014. PhyML reconstructs ancestral sequences. It also implements an efficient simulated annealing algorithm to maximize the likelihood function. These features are available in version 20140911 onwards.

// 23/02/2014. PhyML downloads are now available from https://github.com/stephaneguindon/phyml-downloads/releases

// 14/04/2013. I have released a new, **unstable** version of PhyML (see http://code.google.com/p/phyml/downloads/list, PhyML Development Version). This release provides options to analyse **multi-gene data sets under a very wide variety of new substitution models** (including LG4X amongst many others). Any feedback (to s.guindon@auckland.ac.nz or on the PhyML discussion group, see link on the left) is much appreciated.

// 12/04/2012. **Improved versions of the NNI and SPR search algorithms** are now available in PhyML. In particular, NNI shows much improved performance compared to previous releases. As usual, any feedback is welcome.

// 08/03/2012. PhyML implements a **covarion model**. It is not documented yet but has been tested thoroughly by Salvador Capella. Please feel free to send me an email (s.guindon at auckland.ac.nz) for more information about setting up a PhyML analysis using this model.

// 08/03/2012. The tree **topology constraint** feature of PhyML went through a few tests (thanks to Jaime Huerta Cepas) and is now more stable. It is also documented in PhyML user manual (see download on the left).

// 09/02/2012. PhyML implements tree topology estimation under **user-defined clade constraints**. This option is not documented yet and is only available as a beta version. For people interested in testing it, please feel free to send an email to s.guindon at auckland.ac.nz for more information.