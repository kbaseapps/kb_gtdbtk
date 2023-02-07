# kb_gtdbtk release notes
=========================

1.2.2
_____
* record taxon_assignment as GTDB_R07-RS207 in genome obj

1.2.1
_____
* added Krona chart to report

1.2.0
_____
* Update GTDB-Tk to v2.1.0 (released May 11, 2022)
* Update refdata to GTDB-R07_RS207 (released Apr 8, 2022)
* Update base Docker image to kbase/sdkpython:3.8.10

1.1.0
______
* Return all output from GTDB-Tk Classify workflow as downloadable ZIP archive, especially classify.tree
* Don't accept individual Genomes or Assemblies to avoid wasteful runs.  Require GenomeSet or AssemblySet.

1.0.0
______
* Update GTDB-Tk to v1.7.0
* Added option to retain old Genome taxonomy
* Added updated Genome objects to report

0.1.6
______
* Update GTDB-Tk to v1.6.0
* Update refdata to GTDB R06-RS202
* Update FASTANI to v1.33
* renamed GTDB-Tk Classify App method to run_kb_gtdbtk_classify_wf()
* removed empty values from output html to make datatables happy
* added unit tests for binnedcontigs, faster archaeal assemblySet, genome, and genomeSet
* made input types clear
* update Genome object taxon_assignment and taxonomy (if empty) with GTDB Classification

0.1.5
______
* Update GTDB-Tk to v1.3.0
* Update refdata to GTDB to R05-RS95
* Update FastANI to v1.31

0.1.4
______
* Fix typo in citations

0.1.3
______
* Added citations to the publications section
* Updated description

0.1.2
______
* Update GTDB-Tk to version 1.1.0
* Added GTDB-Tk version number to tooltip


0.1.1
_____
* Update GTDB-Tk version to 1.0.2
* remove python 2.7 environment
* changed input variable inputObjectRef to input_object_ref

0.1.0
_____
* Fixed a bug where workspace objects with pipe in the name would cause GTDB-Tk to fail.
* Fixed a bug where binned contigs with arbitrary characters in the name could cause GTDB-Tk to
  fail.
* Fixed copy of a copy / lost permissions issues.
* Drastically reduced the size of the report zip file, which should speed up load time for
  the HTML report's first load.

0.0.6
_____
* Fix failing build

0.0.5
_____
* Updated to GTDB-Tk v0.3.3

0.0.4
-----
* update app description
* added Donavan Parks as author & owner

0.0.3
-----
* fixed output display, file name scheme changed in new version of gtdbtk
* updated author list of GTDB-Tk manuscript

0.0.2
-----
* minor change

0.0.1
-----
* Initial dev release
* data version 89
* gtdbtk version 0.3.2

0.0.0
-----
* Module created by kb-sdk init
