Installation {#install}
============

<!-- > **Warning!** FFEA is still a work in progress, and may have some rough edges. If you find one of these edges, please get in touch with us using our [issue tracker](https://bitbucket.org/sohpc-ffea/ffea/issues) so that we can round it down. -->

Install instructions for the FFEA package, consisting of the FFEA runner 
 and the FFEA tools can be found in [quickstart](\ref quickstart).
 In addition to that, [troubleshooting](\ref troubleshooting) provides 
  detailed notes on how to install dependencies, as well as some extra help. 

- @subpage quickstart
- @subpage troubleshooting

Once FFEA has been installed, users may want to provide themselves with some 
 extra packages that have proved to be useful at setting up the system
 to simulate, as well as at analysing the results: 

   * [PyMOL](https://www.pymol.org) (>=1.6) can 
        be used to visualise FFEA systems and trajectories
        as well as molecular and EM systems. Alternatives 
        to visualise and worked with molecular systems
        include [Chimera](https://www.cgl.ucsf.edu/chimera/)
        and [VMD](http://www.ks.uiuc.edu/Research/vmd/).


   * [GTS](http://gts.sourceforge.net) (>=0.7.6)[OPTIONAL]. The
     GNU Triangulated Surface Libraries
     allowing the manipulation and coarsening of surface profiles.


   * [NETGEN](https://sourceforge.net/projects/netgen-mesher/) 
   or [TETGEN](http://wias-berlin.de/software/tetgen/) [OPTIONAL]. 
     Programs which convert surface profile into volumetric meshes 
        to be used by FFEA.


   * [pyPcazip](https://pypi.python.org/pypi/pyPcazip) [OPTIONAL]
     Some of the Python FFEA analysis tools interact with these 
     Principal Component Analysis library in order to generate the standard
     PCA output (eigensystems, projections, animations etc)
     obtained from standard from equivalent MD simulations.

Some notes on how to use these tools in relation to FFEA can be found 
 in the [Tutorial](\ref Tutorial). However, mastering these tools 
 may imply consulting the documentation provided by the packages themselves.
