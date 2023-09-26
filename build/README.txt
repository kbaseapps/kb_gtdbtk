Steps in processing GTDB genomes into KBase Genome objects


1. Download latest GTDB species reps data
-----------------------------------------

gtdb_tarball_downloads

untar into
gtdb_downloads/reps/assembly
gtdb_downloads/reps/proteome


2. Create IDs.list files for species reps
-----------------------------------------

    $ ls -1 gtdb_downloads/reps/proteome/archaea/ | grep '^RS_' | sed 's/^RS_//' | sed 's/_protein\.faa//' | sort | uniq > lists/reps/IDs-RS-Archaea.list
    $ ls -1 gtdb_downloads/reps/proteome/archaea/ | grep '^GB_' | sed 's/^GB_//' | sed 's/_protein\.faa//' | sort | uniq > lists/reps/IDs-GB-Archaea.list
    $ ls -1 gtdb_downloads/reps/proteome/bacteria/ | grep '^RS_' | sed 's/^RS_//' | sed 's/_protein\.faa//' | sort | uniq > lists/reps/IDs-RS-Bacteria.list
    $ ls -1 gtdb_downloads/reps/proteome/bacteria/ | grep '^GB_' | sed 's/^GB_//' | sed 's/_protein\.faa//' | sort | uniq > lists/reps/IDs-GB-Bacteria.list


3. Download NCBI GFF+FNA data
-----------------------------

note: s__Fenollaria sp905332505 (GB_GCA_905332505.2) was suppressed so removing from lists.  Fortunately not the rep for a cluster

    e.g.
    $ ./scripts/download_ncbi_gff+fna_files.pl lists/newreps-r214-chunks/IDs-Bacteria-newR214-GB.list.16001-17651 ncbi_downloads-GFF+FNA/GB/bacteria


4. Check GTDB assembly and proteome against NCBI scaffold sequences
-------------------------------------------------------------------
Run scripts/check_gtdb_proteome_gene_coords.pl

    e.g.
    $ ./scripts/check_gtdb_proteome_gene_coords.py -i lists/newreps-r214/IDs-Archaea-newR214-RS.list -a gtdb_downloads_extracted/reps/assembly/ -p gtdb_downloads_extracted/reps/proteome/ -n ncbi_downloads-GFF+FNA/RS/archaea > logs/check_gtdb_proteome_gene_coords/Archaea-newR214-RS.log


5. Add in non-overlapping features from NCBI GFFs
-------------------------------------------------
Run scripts/merge_GTDB_NCBI_genomes.pl

    e.g.
    $ ./scripts/merge_GTDB_NCBI_genomes.pl gtdb_downloads_extracted/reps/proteome/bacteria/ ncbi_downloads-GFF+FNA/RS/bacteria/ gtdb_downloads_extracted/reps/assembly/ merged_GFF/RS/bacteria/ lists/reps/IDs-RS-Bacteria.list.05001-10000 RefSeq overwrite


6. Get Isolate, MAG, and SAG lists
----------------------------------

  A. download metadata etc.

  B. parse out by genome_type (column 56) of metadata
    $ cat gtdb_tarball_downloads/ar53_metadata_r207.tsv | awk -F "\t" '{ print $1 "\t" $56 }' | fgrep 'none' | awk '{ print $1 }' > lists/all/archaea-Isolates.list
    $ cat gtdb_tarball_downloads/ar53_metadata_r207.tsv | awk -F "\t" '{ print $1 "\t" $56 }' | fgrep 'derived from single cell' | awk '{ print $1 }' > lists/all/archaea-SAGs.list
    $ cat gtdb_tarball_downloads/ar53_metadata_r207.tsv | awk -F "\t" '{ print $1 "\t" $56 }' | fgrep 'derived from metagenome' | awk '{ print $1 }' > lists/all/archaea-MAGs.list
    $ cat gtdb_tarball_downloads/ar53_metadata_r207.tsv | awk -F "\t" '{ print $1 "\t" $56 }' | fgrep 'derived from environmental sample' | awk '{ print $1 }' >> lists/all/archaea-MAGs.list

    $ cat gtdb_tarball_downloads/bac120_metadata_r207.tsv | awk -F "\t" '{ print $1 "\t" $56 }' | fgrep 'none' | awk '{ print $1 }' > lists/all/bacteria-Isolates.list
    $ cat gtdb_tarball_downloads/bac120_metadata_r207.tsv | awk -F "\t" '{ print $1 "\t" $56 }' | fgrep 'derived from single cell' | awk '{ print $1 }' > lists/all/bacteria-SAGs.list
    $ cat gtdb_tarball_downloads/bac120_metadata_r207.tsv | awk -F "\t" '{ print $1 "\t" $56 }' | fgrep 'derived from metagenome' | awk '{ print $1 }' > lists/all/bacteria-MAGs.list
    $ cat gtdb_tarball_downloads/bac120_metadata_r207.tsv | awk -F "\t" '{ print $1 "\t" $56 }' | fgrep 'derived from environmental sample' | awk '{ print $1 }' >> lists/all/bacteria-MAGs.list

  C. standardize format and sort order for comm command
    $ cat lists/all/archaea-Isolates.list | sed 's/^GB_//' | sed 's/^RS_//' | sort | uniq > foo; mv foo lists/all/archaea-Isolates.list
    $ cat lists/all/archaea-MAGs.list | sed 's/^GB_//' | sed 's/^RS_//' | sort | uniq > foo; mv foo lists/all/archaea-MAGs.list
    $ cat lists/all/archaea-SAGs.list | sed 's/^GB_//' | sed 's/^RS_//' | sort | uniq > foo; mv foo lists/all/archaea-SAGs.list
    $ cat lists/all/bacteria-Isolates.list | sed 's/^GB_//' | sed 's/^RS_//' | sort | uniq > foo; mv foo lists/all/bacteria-Isolates.list
    $ cat lists/all/bacteria-MAGs.list | sed 's/^GB_//' | sed 's/^RS_//' | sort | uniq > foo; mv foo lists/all/bacteria-MAGs.list
    $ cat lists/all/bacteria-SAGs.list | sed 's/^GB_//' | sed 's/^RS_//' | sort | uniq > foo; mv foo lists/all/bacteria-SAGs.list    

    $ sort lists/reps/IDs-RS-Archaea.list | uniq > foo; mv foo lists/reps/IDs-RS-Archaea.list
    $ sort lists/reps/IDs-GB-Archaea.list | uniq > foo; mv foo lists/reps/IDs-GB-Archaea.list
    $ sort lists/reps/IDs-RS-Bacteria.list | uniq > foo; mv foo lists/reps/IDs-RS-Bacteria.list
    $ sort lists/reps/IDs-GB-Bacteria.list | uniq > foo; mv foo lists/reps/IDs-GB-Bacteria.list

  D. Get lists of the species reps by genome type
    $ comm -12 lists/reps/IDs-RS-Archaea.list lists/all/archaea-Isolates.list > lists/reps/IDs-RS-Archaea-Isolates.list
    $ comm -12 lists/reps/IDs-RS-Archaea.list lists/all/archaea-MAGs.list > lists/reps/IDs-RS-Archaea-MAGs.list
    $ comm -12 lists/reps/IDs-RS-Archaea.list lists/all/archaea-SAGs.list > lists/reps/IDs-RS-Archaea-SAGs.list
    $ comm -12 lists/reps/IDs-GB-Archaea.list lists/all/archaea-Isolates.list > lists/reps/IDs-GB-Archaea-Isolates.list
    $ comm -12 lists/reps/IDs-GB-Archaea.list lists/all/archaea-MAGs.list > lists/reps/IDs-GB-Archaea-MAGs.list
    $ comm -12 lists/reps/IDs-GB-Archaea.list lists/all/archaea-SAGs.list > lists/reps/IDs-GB-Archaea-SAGs.list

    $ comm -12 lists/reps/IDs-RS-Bacteria.list lists/all/bacteria-Isolates.list > lists/reps/IDs-RS-Bacteria-Isolates.list
    $ comm -12 lists/reps/IDs-RS-Bacteria.list lists/all/bacteria-MAGs.list > lists/reps/IDs-RS-Bacteria-MAGs.list
    $ comm -12 lists/reps/IDs-RS-Bacteria.list lists/all/bacteria-SAGs.list > lists/reps/IDs-RS-Bacteria-SAGs.list
    $ comm -12 lists/reps/IDs-GB-Bacteria.list lists/all/bacteria-Isolates.list > lists/reps/IDs-GB-Bacteria-Isolates.list
    $ comm -12 lists/reps/IDs-GB-Bacteria.list lists/all/bacteria-MAGs.list > lists/reps/IDs-GB-Bacteria-MAGs.list
    $ comm -12 lists/reps/IDs-GB-Bacteria.list lists/all/bacteria-SAGs.list > lists/reps/IDs-GB-Bacteria-SAGs.list

  E. Move GCF_000018425.1 (Candidatus Desulforudis audaxviator MP104C) from IDs-RS-Bacteria-Isolates.list to IDs-RS-Bacteria-MAGs.list


7. Generate and Run KBase Genome object import scripts
------------------------------------------------------

  A. create workspace(s)

    for GTDB r207
    * r207 Archaea: wsid=117549, wsname=dylan:narrative_1653154088731
    * r207 Bacteria (no GB MAGs): wsid=117550, wsname=dylan:narrative_1653154121485
    * r207 Bacterial (all GB MAGs): wsid=117551, wsname=dylan:narrative_1653154144334
    * r214 Archaea: wsid=117552, wsname=dylan:narrative_1653154178433
    * r214 Bacteria (no GB MAGs): wsid=117553, wsname=dylan:narrative_1653154252637
    * r214 Bacteria (all GB MAGs): wsid=117554, wsname=dylan:narrative_1653154292963

  B. 1) $ mkdir -p direct_imports/images
	$ mkdir -p direct_imports/params
	$ mkdir -p direct_imports/scripts
	$ mkdir -p direct_imports/staging
     2) create direct_imports/sdk.cfg (don't ever commit this to github!) with:
         kbase_endpoint=https://kbase.us/services
         token=<a fresh token to PROD>
         auth_service_url=https://kbase.us/services/auth/api/legacy/KBase/Sessions/Login
         auth_service_url_allow_insecure=false

  C. for each set, run
      $ ./scripts/prepare_batch_genome_import-GFF+FNA.pl

      e.g.
      $ ./scripts/prepare_batch_genome_import-GFF+FNA.pl -gen lists/reps/IDs-GB-Archaea-MAGs.list -set GTDB_Arc-MAG-GB -obj GTDB_Arc -type MAG -gff /home/ac.dchivian/proj/genome_ref_dbs/GTDB/GTDB_r207/merged_GFF -fna /home/ac.dchivian/proj/genome_ref_dbs/GTDB/GTDB_r207/gtdb_downloads/reps/assembly -imp /home/ac.dchivian/proj/genome_ref_dbs/GTDB/GTDB_r207/direct_imports -rele GTDB_r207 -d archaea -work dylan:narrative_1653154088731

  D. for each set, run the import script
      $ ./direct_imports/scripts/GTDB_Bac-MAG-GB.sh


8. Update KBase Genome objects with additional fields
-----------------------------------------------------

  A. check the header for the metadata file to see that the columns are the same
     $  ./scripts/check_metadata_header.pl gtdb_downloads/ar53_metadata_r214.tsv ../GTDB_r207/gtdb_downloads/ar53_metadata_r207.tsv

  B. 1) $ mkdir -p direct_updates/images
	$ mkdir -p direct_updates/params
	$ mkdir -p direct_updates/scripts
	$ mkdir -p direct_updates/fields
     2) copy direct_imports/sdk.cfg to direct_updates/sdk.cfg

  C. for each set, run
      $ ./scripts/prepare_batch_genome_field_update.pl

      e.g.
      $ ./scripts/prepare_batch_genome_field_update.pl -gen lists/reps/IDs-GB-Archaea-MAGs.list -set GTDB_Arc-MAG-GB -obj GTDB_Arc -metadata gtdb_tarball_downloads/ar53_metadata_r207.tsv -updatedir /home/ac.dchivian/proj/genome_ref_dbs/GTDB/GTDB_r207/direct_updates -rele GTDB_r207 -d archaea -work dylan:narrative_1653154088731

  D. for each set, run the update script
      $ ./direct_updates/scripts/GTDB_Bac-MAG-GB.sh


9. Add functional annotations to features in Genome objects
-----------------------------------------------------------

  A. if not already done, determine which genomes are new in this release that need gene functional annotation.

    $ cat GTDB_r95.0/lists/reps/*.list | sort | uniq > GTDB_r95.0-ALL.list
    $ cat GTDB_r207/lists/reps/*.list | sort | uniq > GTDB_r207-ALL.list
    $ comm -13 GTDB_r95.0-ALL.list GTDB_r207-ALL.list > GTDB_r207-NEW.list

  B. Run eggnog mapper App on GenomeSets in GTDB workspaces
     1) mkdir -p annotations/reps/archaea/grouped
        mkdir -p annotations/reps/bacteria/grouped
     2) Run eggnog mapper App and transfer 6 output files to annotations dir

  C. 1) $ mkdir -p direct_updates_features/images
        $ mkdir -p direct_updates_features/params
        $ mkdir -p direct_updates_features/scripts
        $ mkdir -p direct_updates_features/features
     2) copy direct_imports/sdk.cfg to direct_updates_features/sdk.cfg

  D. for each set, run
      $ ./scripts/prepare_batch_genome_feature_update.pl

      e.g.
      $ ./scripts/prepare_batch_genome_feature_update.pl -updatedir /home/ac.dchivian/proj/genome_ref_dbs/GTDB/GTDB_r214.1/direct_updates_features -set GTDB_Arc-SAGs-RS -egg annotations/reps/archaea/grouped/GTDB_Arc-SAGs-RS.GenomeSet-eggnog5.emapper.annotations.gz -nov annotations/reps/archaea/grouped/GTDB_Arc-SAGs-RS.GenomeSet-novel_fams.emapper.annotations.gz

  E. for each set, run the update script
      $ ./direct_updates_features/scripts/GTDB_Arc-SAGs-RS.sh


10. Add Lineage to Genome and Assembly objects
----------------------------------------------

  add new release taxonomy and taxon_assignments info to genome objs from previous release and add new release std_lineage to genome and assembly objs

  A. 1) $ mkdir -p direct_updates_lineage/images
	$ mkdir -p direct_updates_lineage/params
	$ mkdir -p direct_updates_lineage/scripts
	$ mkdir -p direct_updates_lineage/lineage
     2) copy direct_imports/sdk.cfg to direct_updates/sdk.cfg

  B. for each set, run
      $ ./scripts/prepare_batch_genome_lineage_update.pl

      e.g.
      $ ./scripts/prepare_batch_genome_lineage_update.pl -gen lists/reps/IDs-GB-Archaea-MAGs.list -set GTDB_r214_Arc-MAG-GB -metadata gtdb_tarball_downloads/ar53_metadata_r214.tsv -updatedir /home/ac.dchivian/proj/genome_ref_dbs/GTDB/GTDB_r214.1/direct_updates_lineage -rele GTDB_r214 -work dylan:narrative_1653154178433

  C. for each set, run the update script
      $ ./direct_updates_lineage/scripts/GTDB_Bac-MAG-GB.sh


11. Make genome UPA mapping for GTDB-Tk Classify lookup
-------------------------------------------------------

  run count_ws_obj() to get obj name to UPA map

    for each workspace with GTDB genome objs, do A-C
    A. set the workspace name and run count_ws_objs() from kb_ObjectUtilities for Genomes
    B. capture the STDOUT and reduce to just the mapping lines
    C. cat r207_r214-1.map | sed 's/^.*GTDB_Arc-//' | sed 's/^.*GTDB_Bac-//' | sed 's/\.Genome -> /\t/' | sed 's/\r//' >> foo-1

    D. cat foo-? | sort > /<PATH>/kb_gtdbtk/data/Genome_UPAs-GTDB_r207_r214.tsv
    E. cd /<PATH>/kb_gtdbtk/data
    F. ln -s Genome_UPAs-GTDB_r207_r214.tsv Genome_UPAs-GTDB.tsv
  

12. Make BLAST dbs (for BLASTp App in kb_blast module)
------------------------------------------------------

  A. get mapping from Genome ID to UPA and save as Genome_UPAs-GTDB.tsv (from above)

  B. run scripts/add_UPA_to_genome_faa_header.pl for each group
     (e.g. RS-Archaea, GB-Archaea, RS-Bacteria, GB-Bacteria)

    $ ./scripts/add_UPA_to_genome_faa_header.pl -genomeidsfile lists/reps/IDs-RS-Archaea.list -idtoupafile tables/Genome_UPAs-GTDB_r207.tsv -faaindir gtdb_downloads/reps/proteome/archaea -stagingdir proteomes/reps/Archaea-RS

  C. merge upa-decorated faa files into one per group
    $ cat proteomes/reps/Archaea-RS/*faa > proteomes/reps/Archaea-RS.faa
    $ gzip proteomes/reps/Archaea-RS.faa
    $ rm -rf proteomes/reps/Archaea-RS

  D. copy faa files to public download site at NERSC
    $ scp proteomes/reps/*faa.gz dylan@cori.nersc.gov:/global/cfs/cdirs/kbase/www/GTDB/r207.0/blast_dbs/v1.0.0/


(13. Add circular contig flags to Assembly objects)
---------------------------------------------------

  A. would be NICE TO HAVE.  haven't implemented yet.  would also be nice to not have to break features that go across "edge" of circular contig, but we can't have everything, can we?





