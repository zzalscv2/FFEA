ChangeLog {#changelog}
============

This document aims to track all the notable changes to the FFEA project.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).



2.6.0 - 2017-11-28 {#v260}
=========================


Added
-----

* Performance significantly improved for the PyMOL plugin.

* New short range potential introduced, see [General Soft Potential](\ref gPotential).

* ChangeLog introduced.


Changed
-------

* Keyword ` calc_vdw ` has been changed to ` calc_ssint `,
	the new name reflecting that they are surface-surface interactions.

* Keyword ` calc_steric ` was introduced, to make possible 
	the combination of surface-surface interactions 
   with steric interactions, without constantly extending the 
	entries for keyword ` ssint_type `.


Deprecated
-----------


* Using ` vdw ` instead of ` ssint ` in keywords is deprecated. 
	This way, ` calc_vdw ` should be ` calc_ssint `, 
	` vdw_type ` should be ` ssint_type `, ` vdw_in_fname ` should 
   be ` ssint_in_fname `, et cetera.


