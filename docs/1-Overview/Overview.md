Overview {#overview}
=========

Fluctuating Finite Element Analysis is a new molecular modelling algorithm, designed to support systems that are larger and more complex than those modelled by atomistic molecular dynamics. Instead of modelling biological systems as a collection of connected atoms, it models them as 3D volumes comprised of tetrahedrons. Unlike previous coarse-grained models, the models FFEA generates are visco-elastic continuum solids. Unlike other applications of Finite Element Analysis, these systems are subject to thermal fluctuations.

This technique has the potential to model large, complex systems, made of many molecules, and complex processes at the frontiers of molecular biology. As it does not not require an atomistic level of detail, it can also be used to simulate biological molecules that cannot be imaged using X-ray crystallography.

Features  {#features}
========

 * Protein Interactions:
  * A 6-12 Lennard-Jones Potential.
  * A repulsive potential that is proportional 
        to the overlapping volume.
  * Specific interactions defined using precomputed potentials.
        More documentation can be found [here](\ref fmApproach).
 * Kinetic state changes can be simulated together with the continuum model to
    account for conformational changes and binding events.
 * Conversion tools for EM density data and atomistic structures into FFEA simulations.
 * A plugin for PyMOL, allowing the visualisation of FFEA systems and trajectories.
 * Initialisation and analysis tools available on the command line and under a Python API.
 * [KOBRA model](\ref rods} for slender biological objects such as coiled-coils.
 * Lees-Edwards boundary conditions.
 * An extremely high degree of reproducity - check out our integration tests!


Publications  {#publications}
============

### Methodology: ###
* Welch R. J., Harris S. A., Harlen O. G. & Read D. J. "["KOBRA: A Fluctuating Elastic Rod Model for Slender Biological Macromolecules"](https://doi.org/10.1039/D0SM00491J)" (2020), Soft Matter 16, 7544-7555.
* Solernou A., Hanson B. S., Richardson R. A., Welch R., Harris S. A., Read D. J., Harlen O. G. "["Fluctuating Finite Element Analysis (FFEA): A continuum mechanics software tool for mesoscale simulation of biomolecules"](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005897)" (2018), PLoS Comput. Biol. 14(3): e1005897.
* Oliver R., Read D. J., Harlen O. G. & Harris S. A. "["A Stochastic finite element model for the dynamics of globular macromolecules"](http://www.sciencedirect.com/science/article/pii/S0021999112007589)" (2013), J. Comp. Phys. 239, 147-165.
* Patargias G. N., Harris S. A. & Harding J. "["A demonstration of the inhomogeneity of the local dielectric response of proteins by molecular dynamics simulations"](https://www.ncbi.nlm.nih.gov/pubmed/20572740)" (2010), J. Chem. Phys. 132, 235103.

### Applications: ###
* Richardson R. A., Hanson B. S., Read D. J., Harlen O. G. & Harris S. A. "["Exploring the dynamics of flagellar dynein within the axoneme with Fluctuating Finite Element Analysis"](https://doi.org/10.1017/S0033583520000062)", Q. Rev. Biophys. 53, E9 (2020).
* Lee S. C., Collins R., Lin Y., Jamshad M., Broughton C., Harris S. A., Hanson B. S., Tognoloni C., Parslow R. A., Terry A. E., Rodger A., Smith C. J., Edler K. J., Ford R., Roper D. I. & Dafforn T. R. "["Nano-encapsulated Escherichia coli Divisome Anchor ZipA, and in Complex with FtsZ"](https://doi.org/10.1038/s41598-019-54999-x)" (2019), Sci. Rep. 9, 18712.
* Richardson R., Papachristos K., Read D. J., Harlen O. G., Harrison M. A., Paci E., Muench S. P. & Harris S. A "["Understanding the apparent stator-rotor connections in the rotary ATPase family using coarse-grained computer modelling"](https://www.ncbi.nlm.nih.gov/pubmed/25174610)" (2014), Proteins: Struct., Funct., Bioinf. 82, 3298-3311.

### Reviews: ###
* Gray A., Harlen O. G., Harris S. A., Khalid S., Leung Y. M., Lonsdale R., Mulholland A. J., Pearson A. R., Read D. J. & Richardson R. A. "["In pursuit of an accurate spatial and temporal model of biomolecules at the atomistic level: a perspective on computer simulation"](https://www.ncbi.nlm.nih.gov/pubmed/25615870)" (2015), Acta Cryst. D71, 162-172.
* Hanson B., Richardson R., Oliver R., Read D. J., Harlen O. & Harris S. "["Modelling biomacromolecular assemblies with continuum mechanics"](https://www.ncbi.nlm.nih.gov/pubmed/25849915)" (2015), Biochem. Soc. Trans. 43, 186-192.
* Oliver R. , Richardson R. A., Hanson B., Kendrick K., Read D. J., Harlen O. G. & Harris S. A. "["Modelling the Dynamic Architecture of Biomaterials Using Continuum Mechanics"](http://link.springer.com/chapter/10.1007%2F978-3-319-09976-7_8)" (2014), Protein Modelling, G. Náray-Szabó, Editor. Springer International Publishing. p. 175-197.
       

Getting Started  {#gettingstarted}
===============

FFEA is free to download and use under the GPLv3 software license. We provide binary releases, and building from source is relatively painless.

* Download the most recent [binary release](https://bitbucket.org/FFEA/ffea/downloads/) of FFEA from BitBucket, or [compile from source](cloning the repository) by cloning the repository. For the latest bleeding-edge features, switch to the [development branch](https://bitbucket.org/FFEA/ffea/src/superdev/) instead.
* If you're compiling from source, install FFEA according to the instructions found in the [installation guide](\ref install).
* Once FFEA is installed, consult the [first-time user tutorial](\ref Tutorial). For KOBRA rods, try the [rods tutorial](\ref Tutorial) instead.

       
Technology  {#technology}
============
 
   * [Boost](http://www.boost.org) (>=1.54.0) is used 
     for ease of programming 
     at the initialisation phase. Modules "system", "filesystem" and 
     "program-options" are required. A bundle of version 1.63 is shipped with FFEA.
   * [Eigen](http://eigen.tuxfamily.org) (>=3.2.1).   
     FFEA uses Eigen to calculate and solve linear approximations to the model i.e. Elastic / Dynamic Network Models. CMake will download and use Eigen 3.3.7 if not told otherwise.
   * [RngStreams](http://www.iro.umontreal.ca/~lecuyer/myftp/streams00/)<sup>[1](#RngStreams1)</sup><sup>,[2](#RngStreams2)</sup>
        is shipped with FFEA and used as Random Number Generator (RNG). RngStreams 
        allows the FFEA to **safely** generate random numbers when running 
        on a number of threads, as well as safe restarts, recovering the state 
        of the RNGs in the last saved time step.
   * [Tet_a_tet](https://github.com/erich666/jgt-code/blob/master/Volume_07/Number_2/Ganovelli2002/tet_a_tet.h)<sup>[3](#tetatetpaper)</sup>
        is shipped with FFEA and used to detect element overlapping 
        in the steric repulsion module. 
   * [Doxygen](http://www.doxygen.org) (>= 1.8) [OPTIONAL] 
        is used to generate the documentation. 
   * [PyMOL](https://www.pymol.org) (>=1.8, [though 1.8 is recommended](https://anaconda.org/mw/pymol)) [OPTIONAL] can 
        be used, using the plugin we provide,
        to visualise FFEA systems and trajectories
        as well as molecular and EM systems.
   * [mtTkinter](http://tkinter.unpythonic.net/wiki/mtTkinter) (0.4) is shipped 
        with FFEA and used in the PyMOL plugin, allowing safe threading. 
   * [GTS](http://gts.sourceforge.net) (>=0.7.6)[OPTIONAL]. The
     GNU Triangulated Surface Libraries
     allowing the manipulation and coarsening of surface profiles.
   * [NETGEN](https://sourceforge.net/projects/netgen-mesher/) 
   or [TETGEN](http://wias-berlin.de/software/tetgen/) [OPTIONAL]. 
     Programs which convert surface profile into volumetric meshes 
        to be used by FFEA.
   * [pyPcazip](https://pypi.python.org/pypi/pyPcazip)<sup>[4](#pyPCApaper)</sup>  [OPTIONAL]
     Some of the Python FFEA analysis tools interact with the pyPcazip 
     Principal Component Analysis libraries in order to generate standard
     PCA output(eigensystems, projections, animations etc)
     equivalent to those obtained from equivalent MD simulations.
     See [here](https://bitbucket.org/ramonbsc/pypcazip/wiki/Home) for their wiki.

<a name="RngStreams1">1</a>: P L'Ecuyer, "Good Parameter Sets for Combined Multiple Recursive Random Number Generators" (1999), Oper. Res., 47(1):159-164. <br> 
<a name="RngStreams2">2</a>: P L'Ecuyer et al., "An Objected-Oriented Random-Number Package with Many Long Streams and Substreams" (2002), Oper. Res., 50(6):1073-1075. <br>
<a name="tetatetpaper">3</a>:  F Ganovelli, et al., "Fast tetrahedron-tetrahedron overlap algorithm" (2002), J. Graph. Tools, 7(2):17-25. <br>
<a name="pyPCApaper">4</a>:  A Shkurti, et al., "pyPcazip: A PCA-based toolkit for compression and analysis of molecular simulation data" (2016), SoftwareX, 7:44-50.


Contribute
==========

Do you have a research question that FFEA could help to answer?

   * Try FFEA and let us know how you're using the software.
   * Send bug reports, questions and feature requests to our [issue tracker](https://bitbucket.org/FFEA/ffea/issues?status=new&status=open)
   * [Fork us](https://bitbucket.org/FFEA/ffea/fork) and help develop FFEA!

If you have questions and comments, please contact us! For biophysics, contact Sarah Harris ([S.A.Harris@leeds.ac.uk](mailto:S.A.Harris@leeds.ac.uk)). For finite elements, contact Oliver Harlen ([O.G.Harlen@leeds.ac.uk](mailto:O.G.Harlen@leeds.ac.uk)). For KOBRA and stochastic enquiries, contact Daniel Read ([D.J.Read@leeds.ac.uk](mailto:D.J.Read@leeds.ac.uk)). For software engineering, contact Joanna Leng ([J.Leng@leeds.ac.uk](mailto:J.Leng@leeds.ac.uk)).


FFEA Team  {#FFEAteam}
==========

   * [Oliver Harlen](https://www.maths.leeds.ac.uk/index.php?id=263&uid=1025)
   * [Sarah Harris](http://www.comp-bio.physics.leeds.ac.uk/)
   * [Daniel Read](http://www1.maths.leeds.ac.uk/~djread/)
   * [Joanna Leng](https://www.eps.leeds.ac.uk/computing/staff/1509/joanna-leng)
   * Albert Solernou
   * Ben Hanson
   * Robin Richardson
   * Robin Oliver
   * [Rob Welch](http://robwel.ch/)
   * Tom Ridley
   * Molly Gravett
   * [Ryan Cocking](mailto:bsrctb@leeds.ac.uk)
   * Jarvellis Rogers


Thanks {#thanks}
=======
We want to thank everybody who has helped in making FFEA possible, from 
 summer students, to experimental professors, we would not be here without you:
   * Stan Burgess
   * Stephen Muench
   * Kerrie Smith
   * Thijs van der Heijden
   * Kees Storm
   * Paul van der Schoot 
   * Toni Collis
   * Neelofer Banglawala
   * Jana Boltersdorf
   * Ondřej Vysocký
   * Guanhao Lu
   * Jonathan Boyle
   * Mike Croucher
   * Christopher Woods
   * Katrina Goldman
   * Matthew Faulkner 
   * Ashley Fenton