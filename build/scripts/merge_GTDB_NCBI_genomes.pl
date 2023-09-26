#!/usr/bin/perl
##
## Copyright 2020-2022 Dylan Chivian
##
## This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
##
###############################################################################


###############################################################################
# conf
###############################################################################

# general
$| = 1;                                              # disable stdout buffering
$debug = 1;                                             # chatter while running

# paths
#if (! defined $ENV{'CHIVIAN_HOME'}) {
#    &abort ("must define \$CHIVIAN_HOME in environment");
#}
#if (! defined $ENV{'SCRATCH_DIR'}) {
#    &abort ("must define \$SCRATCH_DIR in environment");
#}
#$src_dir = $ENV{'CHIVIAN_HOME'}."/src";
#$dat_dir = $ENV{'CHIVIAN_HOME'}."/dat";
#$tmp_dir = (-d $ENV{'SCRATCH_DIR'}."/tmp") ? $ENV{'SCRATCH_DIR'}."/tmp" : "/tmp";
$run_dir = $0;
$run_dir =~ s!^(.*)/[^/]+$!$1!;
$run_dir = '.'  if ($run_dir eq $0);
$cwd_dir = `pwd`;  chomp $cwd_dir;  $cwd_dir =~ s!/+$!!;
$src_dir =~ s!/+!/!g;
$dat_dir =~ s!/+!/!g;
$tmp_dir =~ s!/+!/!g;
$run_dir =~ s!/+!/!g;
$cwd_dir =~ s!/+!/!g;

$GTDB_annot_ver = 'GTDB_r207.0';

###############################################################################
# init
###############################################################################

# argv
#my %opts = &getCommandLineOptions ();
#my $file = $opts{file};

if ($#ARGV < 5) {
    print STDERR "usage: $0 <proteome_in_dir> <gff_in_dir> <fna_in_dir> <gff_out_dir> <id_list_file> <GenBank/RefSeq> [FORCE_OVERWRITE] [CROP_OVERHANG]\n";
    exit -1;
}
$pro_dir = shift @ARGV;
$gff_dir = shift @ARGV;
$fna_dir = shift @ARGV;
$out_dir = shift @ARGV;
$id_list_file = shift @ARGV;
$db_src = shift @ARGV;
$force_overwrite = (shift @ARGV) ? 'true' : undef;
$crop_overhang = (shift @ARGV) ? 'true' : undef;

$id_list_buf = &bigFileBufArray ($id_list_file);

@out = ();

###############################################################################
# main
###############################################################################

@genomeIDs = ();
foreach $line (@$id_list_buf) {
    next if ($line =~ /^\s*\#/);
    $line =~ s/\s*\#.*$//;
    next if ($line =~ /^\s*$/);

    ($genomeID) = $line;
    push (@genomeIDs, $genomeID);
}

$cnt = 0;
for $genomeID (@genomeIDs) {
    print STDERR (++$cnt)." READING $genomeID\n";

    @out = ();
    push (@out, "##gff-version 3");
    push (@out, "#!genome-build-accession NCBI_Assembly:".$genomeID);
    push (@out, "#!annotation-source ".$GTDB_annot_ver);

    $genomeID =~ /^(GC\w)\_(\d{3})(\d{3})(\d{3})\.?\d?$/;
    $split_genomeID[0] = $1; 
    $split_genomeID[1] = $2; 
    $split_genomeID[2] = $3; 
    $split_genomeID[3] = $4; 
    
    $gtdb_prefix = ($db_src =~ /^g/i) ? 'GB' : 'RS';
    $fna_file = join ("/", $fna_dir, $split_genomeID[0], $split_genomeID[1], $split_genomeID[2], $split_genomeID[3], $genomeID."_genomic.fna.gz");
    @gff_matches = glob ($gff_dir.'/'.$genomeID."*.gff.gz");
    $gff_file = $gff_matches[0];
    $pro_file = join ("/", $pro_dir, $gtdb_prefix."_".$genomeID."_protein.faa.gz");
    $gff_out_file = join ("/", $out_dir, $gtdb_prefix."_".$genomeID.".gff");

    if (-s $gff_out_file) {
	print STDERR "\tALREADY HAVE $gff_out_file\n";
	if (! $force_overwrite) {
	    print "\tSKIPPING...\n";
	    next;
	} else {
	    print "\tOVERWRITING...\n";
	}
    }

    
    # read NCBI fna scaff IDs and get lengths
    %fna_scaff_IDs = ();
    %fna_scaff_lens = ();
    %alt_fna_scaff_IDs = ();
    %alt2_fna_scaff_IDs = ();
    $last_id = undef;
    $seq = '';
    for $fna_line (@{&bigFileBufArray($fna_file)}) {
	if ($fna_line =~ /^\s*\>\s*(\S+)/) {
	    $scaff_ID = $1;
	    $fna_scaff_IDs{$scaff_ID} = 1;
	    if (! $last_id) {
		$last_id = $scaff_ID;
	    } else {
		if ($seq) {
		    $fna_scaff_lens{$last_id} = length($seq);
		}
		$seq = '';
		$last_id = $scaff_ID;
	    }
	    
	    # need both forms of alternate that GTDB sometimes uses
	    if ($fna_line =~ /(\S*contig_\d+)/) {
		$alt_scaff_ID = $1;
		$alt_fna_scaff_IDs{$alt_scaff_ID} = $scaff_ID;
	    }
	    if ($fna_line =~ /(\S*) (contig_\d+)/) {
		$alt2_scaff_ID = $1.'_'.$2;
		$alt2_fna_scaff_IDs{$alt2_scaff_ID} = $scaff_ID;
	    }
	}
	else {
	    $fna_line =~ s/\s+//g;
	    $seq .= $fna_line;
	}
    }
    if ($seq) {
	$fna_scaff_lens{$last_id} = length($seq);
    }
    $seq = '';
    $last_id = undef;

    
    # read GTDB proteome scaff IDs
    %pro_scaff_IDs = ();
    for $faa_line (@{&bigFileBufArray($pro_file)}) {
	if ($faa_line =~ /^\s*\>\s*(\S+)/) {
	    $scaff_ID = $1;
	    $scaff_ID =~ s/_\d+$//;
	    $pro_scaff_IDs{$scaff_ID} = 1;
	}
    }
    
    # check for unmatched scaff IDs between GTDB proteome and NCBI fna
    %GTDB_to_NCBI_scaff_map = ();
    for $scaff_ID (sort keys %pro_scaff_IDs) {
	if (! $fna_scaff_IDs{$scaff_ID}) {
	    print STDERR "MISSING standard scaffold ID ".$scaff_ID." for genome ".$genomeID."\n";
	    if (! $alt_fna_scaff_IDs{$scaff_ID}) {
		print STDERR "MISSING alternate scaffold ID ".$scaff_ID." for genome ".$genomeID."\n";
		if (! $alt2_fna_scaff_IDs{$scaff_ID}) {
		    print STDERR "MISSING alternate2 scaffold ID ".$scaff_ID." for genome ".$genomeID."\n";
		    exit (-2);
		} else {
		    print STDERR "alternate2 scaffold ID ".$scaff_ID." is ".$alt2_fna_scaff_IDs{$scaff_ID}." for genome ".$genomeID."\n";
		    #push (@out, join ("\t", $genomeID, $scaff_ID, $alt2_fna_scaff_IDs{$scaff_ID}));
		    $GTDB_to_NCBI_scaff_map{$scaff_ID} = $alt2_fna_scaff_IDs{$scaff_ID};
		}
	    } else {
		print STDERR "alternate scaffold ID ".$scaff_ID." is ".$alt_fna_scaff_IDs{$scaff_ID}." for genome ".$genomeID."\n";
		#push (@out, join ("\t", $genomeID, $scaff_ID, $alt_fna_scaff_IDs{$scaff_ID}));
		$GTDB_to_NCBI_scaff_map{$scaff_ID} = $alt_fna_scaff_IDs{$scaff_ID};
	    }
	    #exit (-2);
	}
	#else {
	#    $GTDB_to_NCBI_scaff_map{$scaff_ID} = $scaff_ID;
	#    #print "FOUND scaffold ID ".$scaff_ID." for genome ".$genomeID."\n";
	#    #exit (0);
	#}
    }


    # read CDS from GTDB proteome
    %CDS_from_GTDB_pro = ();
    for $faa_line (@{&bigFileBufArray($pro_file)}) {
	if ($faa_line =~ /^\s*\>\s*(\S+)/) {
	    $scaff_ID = $1;
	    $scaff_ID =~ s/_\d+$//;

	    if ($GTDB_to_NCBI_scaff_map{$scaff_ID}) {
		$scaff_ID = $GTDB_to_NCBI_scaff_map{$scaff_ID};
	    }

	    ($gtdb_CDS_ID, $beg, $end, $strand_num, $info_str) = split (/\s*\#\s*/, $faa_line);
	    $gtdb_CDS_ID =~ s/^\>//;
	    $strand = ($strand_num == -1) ? '-' : '+';
	    @info = split (/\;/, $info_str);
	    $info[0] = 'ID='.$gtdb_CDS_ID;
	    $info_str = join (';', @info);
	    
	    $gff_formatted_line = join ("\t",
					$scaff_ID,
					$GTDB_annot_ver,
					'CDS',
					$beg,
					$end,
					'.',
					$strand,
					'0',
					$info_str);

	    $pos = "$beg,$end,$strand";
	    if (! $CDS_from_GTDB_pro{$scaff_ID}) {
		$CDS_from_GTDB_pro{$scaff_ID} = +{};
	    }
	    if (! $CDS_from_GTDB_pro{$scaff_ID}->{$beg}) {
		$CDS_from_GTDB_pro{$scaff_ID}->{$beg} = +{};
	    }
	    $CDS_from_GTDB_pro{$scaff_ID}->{$beg}->{$pos} = $gff_formatted_line;
	}
    }

    # read CDS from NCBI gff (and capture any species line
    $species_gff_line = undef;
    $ncbi_taxID = undef;
    if (! -s $gff_file) {
	print STDERR "\tNO GFF file from NCBI for $genomeID\n";
    } else {
	%CDS_from_NCBI_gff = ();
	for $gff_line (@{&bigFileBufArray($gff_file)}) {
	    if ($gff_line =~ /^##species .+id=(\d+)$/) {
		$ncbi_taxID = $1;
		$species_gff_line = $gff_line;
	    }
	    next if ($gff_line =~ /^\s*\#/);
	    $gff_line =~ s/\s*\#.*$//;
	    next if ($gff_line =~ /^\s*$/);

	    ($scaff_ID, $annot_src, $feature_type, $beg, $end, $d1, $strand, $d2, $info) = split (/\t/, $gff_line);

	    next  if ($feature_type eq 'region');

	    if ($feature_type eq 'CDS') {
		$pos = "$beg,$end,$strand";
		if (! $CDS_from_NCBI_gff{$scaff_ID}) {
		    $CDS_from_NCBI_gff{$scaff_ID} = +{};
		}
		if (! $CDS_from_NCBI_gff{$scaff_ID}->{$beg}) {
		    $CDS_from_NCBI_gff{$scaff_ID}->{$beg} = +{};
		}
		$CDS_from_NCBI_gff{$scaff_ID}->{$beg}->{$pos} = $gff_line;
	    }
	}
	if ($species_gff_line) {
	    push (@out, $species_gff_line);
	}

	# reconcile NCBI CDS with GTDB CDS  (NOTE: Many many do not perfectly match coords)
	%NCBI_to_GTDB_perfect = ();
	%GTDB_to_NCBI_perfect = ();
	%NCBI_to_GTDB = ();
	%GTDB_to_NCBI = ();
	$min_overlap_frac = 0.5;
	for $scaff_ID (sort keys %CDS_from_NCBI_gff) {
	    if (! $NCBI_to_GTDB{$scaff_ID}) {
		$NCBI_to_GTDB_perfect{$scaff_ID} = +{};
		$NCBI_to_GTDB{$scaff_ID} = +{};
	    }
	    if (! $GTDB_to_NCBI{$scaff_ID}) {
		$GTDB_to_NCBI_perfect{$scaff_ID} = +{};
		$GTDB_to_NCBI{$scaff_ID} = +{};
	    }
	    for $NCBI_beg (sort {$a <=> $b} keys %{$CDS_from_NCBI_gff{$scaff_ID}}) {
		for $NCBI_pos (sort keys %{$CDS_from_NCBI_gff{$scaff_ID}->{$NCBI_beg}}) {
		    if ($CDS_from_GTDB_pro{$scaff_ID}->{$NCBI_beg} &&
			$CDS_from_GTDB_pro{$scaff_ID}->{$NCBI_beg}->{$NCBI_pos}
			) {
			$GTDB_pos = $NCBI_pos;
			$NCBI_to_GTDB_perfect{$scaff_ID}->{$NCBI_pos} = $GTDB_pos;
			$GTDB_to_NCBI_perfect{$scaff_ID}->{$GTDB_pos} = $NCBI_pos;
			$NCBI_to_GTDB{$scaff_ID}->{$NCBI_pos} = $GTDB_pos;
			$GTDB_to_NCBI{$scaff_ID}->{$GTDB_pos} = $NCBI_pos;
		    } else {
			($NCBI_beg2, $NCBI_end, $NCBI_strand) = split (/,/, $NCBI_pos); 
			$best_GTDB_pos = undef;
			$best_overlap_score = 10000000000000;
			for $GTDB_beg (sort {$a <=> $b} keys %{$CDS_from_GTDB_pro{$scaff_ID}}) {
			    for $GTDB_pos (sort keys %{$CDS_from_GTDB_pro{$scaff_ID}->{$GTDB_beg}}) {
				next if ($GTDB_to_NCBI_perfect{$scaff_ID}->{$GTDB_pos});
				($GTDB_beg2, $GTDB_end, $GTDB_strand) = split (/,/, $GTDB_pos); 
				next if ($GTDB_strand ne $NCBI_strand);
				next if ($GTDB_beg > $NCBI_end or $NCBI_beg > $GTDB_end);

				$GTDB_len = $GTDB_end - $GTDB_beg + 1;
				$NCBI_len = $NCBI_end - $NCBI_beg + 1;

				$GTDB_overlap = $GTDB_len;
				if ($NCBI_beg > $GTDB_beg) {
				    $GTDB_overlap -= $NCBI_beg - $GTDB_beg;
				}
				if ($NCBI_end < $GTDB_end) {
				    $GTDB_overlap -= $GTDB_end - $NCBI_end;
				}

				$NCBI_overlap = $NCBI_len;
				if ($GTDB_beg > $NCBI_beg) {
				    $NCBI_overlap -= $GTDB_beg - $NCBI_beg;
				}
				if ($GTDB_end < $NCBI_end) {
				    $NCBI_overlap -= $NCBI_end - $GTDB_end;
				}

				next if ($NCBI_overlap/$NCBI_len < $min_overlap_frac);
				next if ($GTDB_overlap/$GTDB_len < $min_overlap_frac);
				
				if ((abs($NCBI_len-$GTDB_len) + abs($NCBI_overlap-$GTDB_overlap)) < $best_overlap_score) {
				    $best_overlap_score = abs($NCBI_len-$GTDB_len) + abs($NCBI_overlap-$GTDB_overlap);
				    $best_GTDB_pos = $GTDB_pos;
				}
			    }
			}
			$NCBI_to_GTDB{$scaff_ID}->{$NCBI_pos} = ($best_GTDB_pos) ? $best_GTDB_pos : undef;
		    }
		}
	    }
	}

	# reconcile GTDB CDS with NCBI CDS  (NOTE: may not map to same gene)
	for $scaff_ID (sort keys %CDS_from_GTDB_pro) {
	    if (! $GTDB_to_NCBI{$scaff_ID}) {
		$GTDB_to_NCBI{$scaff_ID} = +{};
	    }
	    for $GTDB_beg (sort {$a<=>$b} keys %{$CDS_from_GTDB_pro{$scaff_ID}}) {
		for $GTDB_pos (sort keys %{$CDS_from_GTDB_pro{$scaff_ID}->{$GTDB_beg}}) {
		    next if ($GTDB_to_NCBI_perfect{$scaff_ID}->{$GTDB_pos});

		    ($GTDB_beg2, $GTDB_end, $GTDB_strand) = split (/,/, $GTDB_pos); 
		    $best_NCBI_pos = undef;
		    $best_overlap_score = 10000000000000;
		    for $NCBI_beg (sort {$a<=>$b} keys %{$CDS_from_NCBI_gff{$scaff_ID}}) {
			for $NCBI_pos (sort keys %{$CDS_from_NCBI_gff{$scaff_ID}->{$NCBI_beg}}) {
			    next if ($NCBI_to_GTDB_perfect{$scaff_ID}->{$NCBI_pos});
			    ($NCBI_beg2, $NCBI_end, $NCBI_strand) = split (/,/, $NCBI_pos); 
			    next if ($GTDB_strand ne $NCBI_strand);
			    next if ($GTDB_beg > $NCBI_end or $NCBI_beg > $GTDB_end);

			    $GTDB_len = $GTDB_end - $GTDB_beg + 1;
			    $NCBI_len = $NCBI_end - $NCBI_beg + 1;

			    $GTDB_overlap = $GTDB_len;
			    if ($NCBI_beg > $GTDB_beg) {
				$GTDB_overlap -= $NCBI_beg - $GTDB_beg;
			    }
			    if ($NCBI_end < $GTDB_end) {
				$GTDB_overlap -= $GTDB_end - $NCBI_end;
			    }

			    $NCBI_overlap = $NCBI_len;
			    if ($GTDB_beg > $NCBI_beg) {
				$NCBI_overlap -= $GTDB_beg - $NCBI_beg;
			    }
			    if ($GTDB_end < $NCBI_end) {
				$NCBI_overlap -= $NCBI_end - $GTDB_end;
			    }
				
			    next if ($NCBI_overlap/$NCBI_len < $min_overlap_frac);
			    next if ($GTDB_overlap/$GTDB_len < $min_overlap_frac);

			    if ((abs($NCBI_len-$GTDB_len) + abs($NCBI_overlap-$GTDB_overlap)) < $best_overlap_score) {
				$best_overlap_score = abs($NCBI_len-$GTDB_len) + abs($NCBI_overlap-$GTDB_overlap);
				$best_NCBI_pos = $NCBI_pos;
			    }
			}
		    }
		    $GTDB_to_NCBI{$scaff_ID}->{$GTDB_pos} = ($best_NCBI_pos) ? $best_NCBI_pos : undef;
		}
	    }
	}

	# check that they are reciprocal best matches, otherwise set to null
	for $scaff_ID (sort keys %GTDB_to_NCBI) {
	    for $GTDB_pos (sort keys %{$GTDB_to_NCBI{$scaff_ID}}) {
		next if ($GTDB_to_NCBI_perfect{$scaff_ID}->{$GTDB_pos});
		next if (! $GTDB_to_NCBI{$scaff_ID}->{$GTDB_pos});
		$NCBI_pos = $GTDB_to_NCBI{$scaff_ID}->{$GTDB_pos};
		next if (! $NCBI_to_GTDB{$scaff_ID}->{$NCBI_pos});
		if ($GTDB_to_NCBI{$scaff_ID}->{$GTDB_pos} ne $NCBI_pos or
		    $NCBI_to_GTDB{$scaff_ID}->{$NCBI_pos} ne $GTDB_pos
		    ) {
		    $GTDB_to_NCBI{$scaff_ID}->{$GTDB_pos} = undef;
		    $NCBI_to_GTDB{$scaff_ID}->{$NCBI_pos} = undef;
		}
	    }
	}
	for $scaff_ID (sort keys %NCBI_to_GTDB) {
	    for $NCBI_pos (sort keys %{$NCBI_to_GTDB{$scaff_ID}}) {
		next if ($NCBI_to_GTDB_perfect{$scaff_ID}->{$NCBI_pos});
		next if (! $NCBI_to_GTDB{$scaff_ID}->{$NCBI_pos});
		$GTDB_pos = $NCBI_to_GTDB{$scaff_ID}->{$NCBI_pos};
		next if (! $GTDB_to_NCBI{$scaff_ID}->{$GTDB_pos});
		if ($GTDB_to_NCBI{$scaff_ID}->{$GTDB_pos} ne $NCBI_pos or 
		    $NCBI_to_GTDB{$scaff_ID}->{$NCBI_pos} ne $GTDB_pos
		    ) {
		    $GTDB_to_NCBI{$scaff_ID}->{$GTDB_pos} = undef;
		    $NCBI_to_GTDB{$scaff_ID}->{$NCBI_pos} = undef;
		}
	    }
	}
	# DEBUG
	#for $scaff_ID (sort keys %NCBI_to_GTDB) {
	#    for $NCBI_pos (sort keys %{$NCBI_to_GTDB{$scaff_ID}}) {
	#	if (! $NCBI_to_GTDB{$scaff_ID}->{$NCBI_pos}) {
	#	    print STDERR "MISSING $genomeID $scaff_ID no match for NCBI_pos $NCBI_pos\n";
	#	}
	#	$GTDB_pos = $NCBI_to_GTDB{$scaff_ID}->{$NCBI_pos};
	#	if ($GTDB_pos ne $NCBI_pos) {
	#	    print STDERR "MAPPED $genomeID $scaff_ID NCBI_pos: $NCBI_pos GTDB_pos: $GTDB_pos\n";
	#	}
	#    }
	#}
	#for $scaff_ID (sort keys %GTDB_to_NCBI) {
	#    for $GTDB_pos (sort keys %{$GTDB_to_NCBI{$scaff_ID}}) {
	#	if (! $GTDB_to_NCBI{$scaff_ID}->{$GTDB_pos}) {
	#	    print STDERR "MISSING $genomeID $scaff_ID no match for GTDB_pos $GTDB_pos\n";
	#	}
	#    }
	#}

	# Now that we have the equivalenced CDSes, we can build the merged GFF
    }


    # no NCBI gff.  just build GFF from GTDB CDSes (value stored in CDS_from_GTDB_pro is gff_line
    %out_gff_recs = ();
    if (! -s $gff_file) {
	for $scaff_ID (sort keys %CDS_from_GTDB_pro) {
	    if (! $out_gff_recs{$scaff_ID}) {
		$out_gff_recs{$scaff_ID} = +{};
	    }
	    for $beg (sort {$a <=> $b} keys %{$CDS_from_GTDB_pro{$scaff_ID}}) {
		if (! $out_gff_recs{$scaff_ID}->{$beg}) {
		    $out_gff_recs{$scaff_ID}->{$beg} = +{};
		}
		for $pos (sort keys %{$CDS_from_GTDB_pro{$scaff_ID}->{$beg}}) {
		    $out_gff_recs{$scaff_ID}->{$beg}->{$pos} = $CDS_from_GTDB_pro{$scaff_ID}->{$beg}->{$pos};
		}
	    }
	}
    }
    # merge GTDB and NCBI features, favoring GTDB CDSes but adding CDS info from NCBI.
    else {
	%out_gff_recs = ();
	%parent_gene_info = ();
	
	# start with GTDB CDS records (we'll add info also later from those with NCBI matches)
	for $scaff_ID (sort keys %CDS_from_GTDB_pro) {
	    if (! $out_gff_recs{$scaff_ID}) {
		$out_gff_recs{$scaff_ID} = +{};
	    }
	    for $beg (sort {$a <=> $b} keys %{$CDS_from_GTDB_pro{$scaff_ID}}) {
		if (! $out_gff_recs{$scaff_ID}->{$beg}) {
		    $out_gff_recs{$scaff_ID}->{$beg} = +{};
		}
		for $pos (sort keys %{$CDS_from_GTDB_pro{$scaff_ID}->{$beg}}) {
		    $out_gff_recs{$scaff_ID}->{$beg}->{$pos} = $CDS_from_GTDB_pro{$scaff_ID}->{$beg}->{$pos};
		}
	    }
	}
	
	# add in features from NCBI GFF
	#  Note: skipping scaffold region, parent genes (except GeneID dbxref), and exons
	#  Note: replace NCBI CDS with GTDB CDS, unless no mapping, then keep NCBI CDS
	#  Note: add info from NCBI CDS to equivalent GTDB CDS
	#
	for $gff_line (@{&bigFileBufArray($gff_file)}) {
	    next if ($gff_line =~ /^\s*\#/);
	    $gff_line =~ s/\s*\#.*$//;
	    next if ($gff_line =~ /^\s*$/);

	    ($scaff_ID, $annot_src, $feature_type, $beg, $end, $d1, $strand, $d2, $info_str) = split (/\t/, $gff_line);
	    $pos = join (",", $beg, $end, $strand);

	    # skip region, exon, and signal_peptide_region_of_CDS
	    next  if ($feature_type eq 'region');
	    next  if ($feature_type eq 'exon');
	    next  if ($feature_type eq 'signal_peptide_region_of_CDS');

	    # instantiate
	    if (! $out_gff_recs{$scaff_ID}) {
		$out_gff_recs{$scaff_ID} = +{};
	    }

	    # just add line for all features other than gene and CDS
	    if ($feature_type ne 'gene' and $feature_type ne 'CDS') {
		if (! $out_gff_recs{$scaff_ID}->{$beg}) {
		    $out_gff_recs{$scaff_ID}->{$beg} = +{};
		}
		$out_gff_recs{$scaff_ID}->{$beg}->{$pos} = $gff_line;
		next;
	    }
	    
	    # capture info from gene to add to CDS (usually precedes CDS, so should be populated)
	    if ($feature_type eq 'gene') {
		if (! $parent_gene_info{$scaff_ID}) {
		    $parent_gene_info{$scaff_ID} = +{};
		}
		$parent_gene_info{$scaff_ID}->{$pos} = $info_str;
		next;
	    }

	    # adjust CDS rec to any matching GTDB CDS and add info from parent gene
	    if ($feature_type eq 'CDS') {
		$this_gff_line = $gff_line;
		
		# info from NCBI CDS
		@ncbi_cds_info = split (/\;/, $info_str);

		# info from NCBI parent
		@ncbi_parent_info = ();
		if ($parent_gene_info{$scaff_ID}->{$pos}) {
		    @ncbi_parent_info = split (/\;/, $parent_gene_info{$scaff_ID}->{$pos});
		}

		# info from GTDB CDS and replace NCBI gff_line with GTDB gff_line
		@gtdb_cds_info = ();
		$gtdb_id = undef;
		if ($NCBI_to_GTDB{$scaff_ID}->{$pos}) {
		    $GTDB_pos = $NCBI_to_GTDB{$scaff_ID}->{$pos};
		    ($GTDB_beg, $GTDB_end, $strand) = split (/,/, $GTDB_pos);

		    # key line to reset gff CDS rec to GTDB's version
		    $this_gff_line = $CDS_from_GTDB_pro{$scaff_ID}->{$GTDB_beg}->{$GTDB_pos};
		    ($scaff_ID, $annot_src, $feature_type, $beg, $end, $d1, $strand, $d2, $info_str) = split (/\t/, $this_gff_line);
		    @gtdb_cds_info = split (/\;/, $info_str);
		    $gtdb_id = $gtdb_cds_info[0];
		}

		# adjust info
		@adjusted_info = ();
		%db_xref = ();
		%names = ();
		$locus_tag = undef;
		$ncbi_cds_id = undef;
		$ncbi_parent_id = undef;
		if (@gtdb_cds_info) {
		    push (@adjusted_info, @gtdb_cds_info);
		}
		if (@ncbi_cds_info) {
		    for $info_pair (@ncbi_cds_info) {
			($info_key, $info_val) = split (/=/, $info_pair);
			if ($info_key eq 'ID' or $info_key eq 'Parent') {
			    if ($info_key eq 'ID') {
				$ncbi_cds_id = $info_val;
				$ncbi_cds_id =~ s/^cds-//i;
			    }
			    next;
			} elsif ($info_key =~ /^dbxref$/i) {
			    foreach $db_xref_pair (split(/,/, $info_val)) {
				($db_xref_name, $db_xref_id) = split (/:/, $db_xref_pair);
				$db_xref{$db_xref_name} = $db_xref_id;
			    } 
			} elsif ($info_key =~ /^name$/i) {
			    $names{$info_val} = 1;
			} elsif ($info_key =~ /^locus_tag$/i) {
			    $locus_tag = $info_val;
			} else {
			    push (@adjusted_info, $info_pair);
			}
		    }
		}
		if (@ncbi_parent_info) {
		    for $info_pair (@ncbi_parent_info) {
			($info_key, $info_val) = split (/=/, $info_pair);
			if ($info_key eq 'ID' or $info_key eq 'gbkey' or $info_key eq 'gene_biotype') {
			    if ($info_key eq 'ID') {
				$ncbi_parent_id = $info_val;
				$ncbi_parent_id =~ s/^gene-//i;
			    }
			    next;
			} elsif ($info_key =~ /^dbxref$/i) {
			    foreach $db_xref_pair (split(/,/, $info_val)) {
				($db_xref_name, $db_xref_id) = split (/:/, $db_xref_pair);
				$db_xref{$db_xref_name} = $db_xref_id;
			    } 
			} elsif ($info_key =~ /^name$/i) {
			    $names{$info_val} = 1;
			} elsif ($info_key =~ /^locus_tag$/i) {
			    if (! $locus_tag) {
				$locus_tag = $info_val;
			    }
			} else {
			    push (@adjusted_info, $info_pair);
			}
		    }
		}

		# set name and add to @adjusted_info
		$name = undef;
		if ($locus_tag) {
		    $name = $locus_tag;
		} elsif (%names) {
		    for $this_name (keys %names) {
			$name = $this_name;
		    }
		} elsif ($ncbi_cds_id) {
		    $name = $ncbi_cds_id;
		} elsif ($ncbi_parent_id) {
		    $name = $ncbi_parent_id;
		} elsif ($gtdb_id) {
		    $name = $gtdb_id;
		} else {
		    print STDERR "unable to find name in genome $genomeID gff_line: '$gff_line'\n";
		    exit (-2);
		}
		
		# set db_xref
		$db_xref_val = undef;
		@db_xrefs = ();
		for $db_xref_name (sort keys %db_xref) {
		    push (@db_xrefs, join (':', $db_xref_name, $db_xref{$db_xref_name}));
		}
		if (@db_xrefs) {
		    $db_xref_val = join (',', @db_xrefs);
		}

		# add ID (if not from GTDB rec), db_xrefs, name, and locus_tag to 
		if (! $gtdb_id) {
		    $this_id = $name;
		    unshift (@adjusted_info, 'ID='.$this_id);
		}
		if ($db_xref_val) {
		    push (@adjusted_info, 'Dbxref='.$db_xref_val);
		}
		if ($name) {
		    push (@adjusted_info, 'Name='.$name);
		}
		if ($locus_tag) {
		    push (@adjusted_info, 'locus_tag='.$locus_tag);
		}		

		# recreate and store gff_line with adjusted info
		$adjusted_info_str = join (';', @adjusted_info);
		
		($scaff_ID, $annot_src, $feature_type, $beg, $end, $d1, $strand, $d2, $info_str) = split (/\t/, $this_gff_line);
		$pos = join (',',$beg,$end,$strand);
		$adjusted_info_gff_line = join ("\t",
						$scaff_ID,
						$annot_src,
						'CDS',
						$beg,
						$end,
						'.',
						$strand,
						'0',
						$adjusted_info_str);

		if (! $out_gff_recs{$scaff_ID}->{$beg}) {
		    $out_gff_recs{$scaff_ID}->{$beg} = +{};
		}
		$out_gff_recs{$scaff_ID}->{$beg}->{$pos} = $adjusted_info_gff_line;
	    }
	}
    }


    # for circular genomes, fix genes that span the end of the scaffold
    %new_out_gff_recs = ();
    for $scaff_ID (sort keys %out_gff_recs) {
	$new_out_gff_recs{$scaff_ID} = +{};
	for $beg (sort {$a <=> $b} keys %{$out_gff_recs{$scaff_ID}}) {
	    $new_out_gff_recs{$scaff_ID}->{$beg} = +{};

	    for $pos (sort keys %{$out_gff_recs{$scaff_ID}->{$beg}}) {

		($scaff_ID, 
		 $annot_src,
		 $feature_type,
		 $this_beg,
		 $end,
		 $d1,
		 $strand,
		 $d2,
		 $info_str
		 ) = split (/\t/, $out_gff_recs{$scaff_ID}->{$beg}->{$pos});

		if (! $fna_scaff_lens{$scaff_ID}) {
		    &abort ("$scaff_ID has no length");
		}
		if ($end <= $fna_scaff_lens{$scaff_ID}) {
		    $new_out_gff_recs{$scaff_ID}->{$beg}->{$pos} = $out_gff_recs{$scaff_ID}->{$beg}->{$pos};
		}
		elsif ($beg > $fna_scaff_lens{$scaff_ID}) {
		    &alert ("feature begins after length of scaffold (".$fna_scaff_lens{$scaff_ID}.") in genome $genomeID: gff_rec:\n'".$out_gff_recs{$scaff_ID}->{$beg}->{$pos});
		    &alert ("SKIPPING Feature");
		    next;
		}
		else {
		    &alert ("feature exceeds length of scaffold (".$fna_scaff_lens{$scaff_ID}.") in genome $genomeID: gff_rec:\n'".$out_gff_recs{$scaff_ID}->{$beg}->{$pos});
		    
		    if ($crop_overhang) {
			&alert ("CROPPING");
		    } else {
			&alert("BREAKING into two parts");
		    }
		    $part_1_beg = $beg;
		    $part_1_end = $fna_scaff_lens{$scaff_ID};
		    $part_2_beg = 1;
		    $part_2_end = $end - $fna_scaff_lens{$scaff_ID};
		    
		    # DEBUG
		    #print "$beg, $end, $strand: part1: $part_1_beg..$part_1_end , part2: $part_2_beg..$part_2_end\n";
		    
		    @info = split (/\;/, $info_str);
		    ($info_key, $info_val) = split (/=/, $info[0]);
		    if ($info_key ne 'ID') {
			&abort ("bad info string.  can't find ID in gff_line:\n'".$out_gff_recs{$scaff_ID}->{$beg}->{$pos}."'");
		    }
		    $part_1_ID = $info_val.'_p1';
		    $part_2_ID = $info_val.'_p2';
		    @info_1 = @info;
		    @info_2 = @info;
		    $info_1[0] = join("=",'ID',$part_1_ID);
		    $info_2[0] = join("=",'ID',$part_2_ID);
		    
		    if (! $new_out_gff_recs{$scaff_ID}->{$part_1_beg}) {
			$new_out_gff_recs{$scaff_ID}->{$part_1_beg} = +{};
		    }
		    if (! $new_out_gff_recs{$scaff_ID}->{$part_2_beg}) {
			$new_out_gff_recs{$scaff_ID}->{$part_2_beg} = +{};
		    }
		    
		    $pos_1 = join(",",$part_1_beg,$part_1_end,$strand);
		    $pos_2 = join(",",$part_2_beg,$part_2_end,$strand);
		    
		    # add back part 1
		    $new_out_gff_recs{$scaff_ID}->{$part_1_beg}->{$pos_1} =
			join ("\t", 
			      $scaff_ID, 
			      $annot_src, 
			      $feature_type,
			      $part_1_beg,
			      $part_1_end,
			      $d1,
			      $strand,
			      $d2,
			      join ("\;", @info_1));
		    
		    &alert ("PART 1: \n'".$new_out_gff_recs{$scaff_ID}->{$part_1_beg}->{$pos_1}."'");

		    if (! $crop_overhang) {
			# add back part 2
			$new_out_gff_recs{$scaff_ID}->{$part_2_beg}->{$pos_2} =
			    join ("\t", 
				  $scaff_ID, 
				  $annot_src, 
				  $feature_type,
				  $part_2_beg,
				  $part_2_end,
				  $d1,
				  $strand,
				  $d2,
				  join ("\;", @info_2));
			
			&alert ("PART 2: \n'".$new_out_gff_recs{$scaff_ID}->{$part_2_beg}->{$pos_2}."'")
		    }
		}
	    }
	}
    }
    %out_gff_recs = %new_out_gff_recs;


    # Finally! build gff output buffer and add back in parent gene recs
    for $scaff_ID (sort keys %out_gff_recs) {
	for $beg (sort {$a <=> $b} keys %{$out_gff_recs{$scaff_ID}}) {
	    for $pos (sort keys %{$out_gff_recs{$scaff_ID}->{$beg}}) {

		($scaff_ID, 
		 $annot_src,
		 $feature_type,
		 $this_beg,
		 $end,
		 $d1,
		 $strand,
		 $d2,
		 $info_str
		 ) = split (/\t/, $out_gff_recs{$scaff_ID}->{$beg}->{$pos});

		# check for any ends less than begins
		if ($end < $beg) {
		    &abort ("bad gff ref (end < beg) for genome $genomeID.\ngff_rec: '".$out_gff_recs{$scaff_ID}->{$beg}->{$pos}."'");
		}
		
		# prepare info to adjust
		@info = split (/\;/, $info_str);
		@this_info = ();
		@parent_info = ();
		$rec_id = undef;
		$parent_id = undef;
		$gb_key = undef;
		$transcript_gbkey = undef;
		$gene_biotype = $feature_type;

		# feature recs that doesn't need adjustment nor parent gene
		if ($feature_type ne 'CDS' and
		    $feature_type ne 'rRNA' and
		    $feature_type ne 'tRNA' and
		    $feature_type ne 'tmRNA' and
		    $feature_type ne 'ncRNA' and
		    $feature_type ne 'miRNA' and
		    $feature_type ne 'snRNA' and
		    $feature_type ne 'snoRNA' and
		    $feature_type ne 'SRP_RNA' and
		    $feature_type ne 'RNase_P_RNA' and
		    $feature_type ne 'antisense_RNA' and
		    $feature_type ne 'ribozyme' and
		    $feature_type ne 'hammerhead_ribozyme' and
		    $feature_type ne 'pseudogenic_rRNA' and
		    $feature_type ne 'pseudogenic_tRNA' and
		    $feature_type ne 'autocatalytically_spliced_intron' and
		    $feature_type ne 'transcript'
		    ) {

		    # check if it has a parent gene in info to add rec for
		    $parent_id = undef;
		    for ($i=0; $i <= $#info; ++$i) {
			($info_key, $info_val) = split (/=/, $info[$i]);
			if ($info_key eq 'Parent') {
			    $parent_id = $info_val;
			}
		    }
		    if ($parent_id) {
			@parent_info = ();
			push (@parent_info, join('=', 'ID', $parent_id));
			push (@parent_info, join('=', 'gbkey', 'Gene'));
			for ($i=0; $i <= $#info; ++$i) {
			    ($info_key, $info_val) = split (/=/, $info[$i]);
                            if ($info_key eq 'locus_tag' or
                                $info_key eq 'old_locus_tag' or
                                $info_key eq 'Dbxref' or
                                $info_key eq 'Name' or
                                $info_key eq 'gene' or
                                $info_key eq 'pseudo'
                                ) {
                                push (@parent_info, $info[$i]);
			    } elsif ($info_key eq 'gbkey') {
				push (@parent_info, join('=', 'gene_biotype', $info_val));
			    }
			}			    
			push (@out, join ("\t", 
					  $scaff_ID, 
					  $annot_src, 
					  'gene',
					  $beg,
					  $end,
					  '.',
					  $strand,
					  '.',
					  join ("\;", @parent_info)));
			@parent_info = ();
		    }
		    # add original unchanged rec
		    push (@out, $out_gff_recs{$scaff_ID}->{$beg}->{$pos});
		}

		# else it's a more common feature type we want to handle carefully
		else {
		    for ($i=0; $i <= $#info; ++$i) {
			($info_key, $info_val) = split (/=/, $info[$i]);
			if ($i==0 and $info_key ne 'ID') {
			    &abort ("missing ID in genome $genomeID gff_rec: '".$out_gff_recs{$scaff_ID}->{$beg}->{$pos});
			}
			if ($info_key eq 'ID') {
			    if ($feature_type eq 'CDS') {
				$info_val =~ s/^cds-//;
				$rec_id = 'cds-'.$info_val;
			    } 
			    elsif ($feature_type eq 'rRNA' or
				   $feature_type eq 'tRNA' or
				   $feature_type eq 'tmRNA' or
				   $feature_type eq 'ncRNA' or		
				   $feature_type eq 'miRNA' or
				   $feature_type eq 'snRNA' or
				   $feature_type eq 'snoRNA' or
				   $feature_type eq 'SRP_RNA' or
				   $feature_type eq 'RNase_P_RNA' or
				   $feature_type eq 'antisense_RNA' or
				   $feature_type eq 'ribozyme' or
				   $feature_type eq 'hammerhead_ribozyme' or
				   $feature_type eq 'pseudogenic_rRNA' or
				   $feature_type eq 'pseudogenic_tRNA' or
				   $feature_type eq 'autocatalytically_spliced_intron' or
				   $feature_type eq 'transcript'
				) {
				$info_val =~ s/^rna-//;
				$rec_id = 'rna-'.$info_val;
			    }
			    push (@this_info, join ('=', 'ID', $rec_id));
			    if ($feature_type ne 'transcript') {
				$parent_id = 'gene-'.$info_val;
				push (@parent_info, join ('=', 'ID', $parent_id));
			    }
			    next;
			} elsif ($info_key eq 'Parent') {
			    next;
			} elsif ($info_key eq 'gbkey') {
			    if ($feature_type eq 'transcript') {
				$transcript_gbkey = $info_val;
			    }
			    next;
			} elsif ($info_key eq 'gene_biotype') {
			    next;
			} else {
			    push (@this_info, $info[$i]);
			    if ($info_key eq 'locus_tag' or
				$info_key eq 'old_locus_tag' or
				$info_key eq 'Dbxref' or
				$info_key eq 'Name' or 
				$info_key eq 'gene' or
				$info_key eq 'pseudo'
				) {
				push (@parent_info, $info[$i]);
			    }
			}
		    }
		    
		    # add gbkey to both infos and gene_biotype to parent info
		    if ($feature_type eq 'CDS') {
			$gb_key = 'CDS';
			$gene_biotype = 'protein_coding';
		    } elsif ($feature_type eq 'rRNA' or
			     $feature_type eq 'pseudogenic_rRNA') {
			$gb_key = 'rRNA';
			if ($feature_type eq 'pseudogenic_rRNA') {
			    $gene_biotype = 'rRNA_pseudogene';
			}
		    } elsif ($feature_type eq 'tRNA' or 
			     $feature_type eq 'pseudogenic_tRNA') {
			$gb_key = 'tRNA';
			if ($feature_type eq 'pseudogenic_tRNA') {
			    $gene_biotype = 'tRNA_pseudogene';
			}
		    } elsif ($feature_type eq 'tmRNA' or
			     $feature_type eq 'ncRNA' or
			     $feature_type eq 'miRNA' or
			     $feature_type eq 'snRNA' or
			     $feature_type eq 'snoRNA' or
			     $feature_type eq 'SRP_RNA' or
			     $feature_type eq 'RNase_P_RNA' or
			     $feature_type eq 'antisense_RNA' or
			     $feature_type eq 'ribozyme' or
			     $feature_type eq 'hammerhead_ribozyme' or 
			     $feature_type eq 'autocatalytically_spliced_intron'
			) {
			$gb_key = 'ncRNA';
		    } elsif ($feature_type eq 'transcript') {
			$gb_key = $transcript_gbkey;
		    }
		    push (@this_info, join ('=','gbkey',$gb_key));

		    # add parent gene and connection
		    if ($feature_type ne 'transcript') {
			push (@parent_info, join ('=','gbkey','Gene'));
			push (@parent_info, join ('=','gene_biotype',$gene_biotype));
		    
			# add linkage from feature rec to parent rec
			push (@this_info, join ('=', 'Parent', $parent_id));
		    
			# add parent gene rec
			push (@out, join ("\t", 
					  $scaff_ID, 
					  $annot_src, 
					  'gene',
					  $beg,
					  $end,
					  '.',
					  $strand,
					  '.',
					  join ("\;", @parent_info)));
		    }
		    
		    # add adjusted feature rec
		    push (@out, join ("\t", 
				      $scaff_ID, 
				      $annot_src, 
				      $feature_type,
				      $beg,
				      $end,
				      $d1,
				      $strand,
				      $d2,
				      join ("\;", @this_info)));

		}
	    }
	}
    }
    push (@out, "###");

    
    # write it
    $out_file = $gff_out_file;
    if ($out_file) {
	open (OUT, '>'.$out_file);
	select (OUT);
    }
    print join ("\n", @out)."\n"  if (@out);
    if ($out_file) {
	close (OUT);
	select (STDOUT);
    }
}


# done
exit 0;

###############################################################################
# subs
###############################################################################

sub by_cnt_then_lineage {
    if ($cnt_by_tax{$b} == $cnt_by_tax{$a}) {
	return $a cmp $b;
    }
    return $cnt_by_tax{$b} <=> $cnt_by_tax{$a};
}

# getCommandLineOptions()
#
#  rets: \%opts  pointer to hash of kv pairs of command line options
#
sub getCommandLineOptions {
    use Getopt::Long;
    my $usage = qq{usage: $0 -file <file>\n};

    # Get args
    my %opts = ();
    &GetOptions (\%opts, "file=s");

    # Check for legal invocation
    if (! defined $opts{file}
        ) {
        print STDERR "$usage\n";
        exit -1;
    }
    &checkExist ('f', $opts{file});

    return %opts;
}

###############################################################################
# util
###############################################################################

# chompEnds ()
#
sub chompEnds {
    my $str = shift;
    $str =~ s/^\s+|\s+$//g;
    return $str;
}


# chip (chop for front of strings)
#
sub chip {
    my @flo = ();
    for ($i=0; $i <= $#_; ++$i) {
        $flo[$i] = substr ($_[$i], 0, 1);
        $_[$i] = substr ($_[$i], 1);                   # don't think this works
    }
    return $flo[0]  if ($#_ == 0);
    return @flo;
}


# chimp (chomp for front of strings)
#
sub chimp {
    my @flo = ();
    for ($i=0; $i <= $#_; ++$i) {
        $_[$i] =~ s/^(\s*)//;                          # don't think this works
        $flo[$i] = $1;
    }
    return $flo[0]  if ($#_ == 0);
    return @flo;
}


# cleanStr ()
#
sub cleanStr {
    my $str = shift;
    $str =~ s/[\x00-\x08\x0B-\x1F\x80-\xFF]//g;
    return $str;
}


# log base 10
#
sub log10 {
    my $n = shift;
    return log($n)/log(10);
}                  


# log base 2
#
sub log2 {
    my $n = shift;
    return log($n)/log(2);
}                  


# base36 ()
#
sub base36 {
    my $d = shift;
    return $d  if ($d =~ /^\d$/);
    return chr ($d - 10 + ord ('A'));
}


# charToHex ()
#
sub charToHex {
    my $ascii = ord($_[0]);
    my %hexMap = (  0 => '0',
                    1 => '1',
                    2 => '2',
                    3 => '3',
                    4 => '4',
                    5 => '5',
                    6 => '6',
                    7 => '7',
                    8 => '8',
                    9 => '9',
                   10 => 'a',
                   11 => 'b',
                   12 => 'c',
                   13 => 'd',
                   14 => 'e',
                   15 => 'f'
                    );

    return $hexMap{(($ascii & 0xf0) >> 4)} . $hexMap{($ascii & 0x0f)};
}
# end charToHex ()


#  hexToChar ()
#
sub hexToChar {
    my $ascii = hex($_[0]);
    return chr $ascii;
}
# end hexToChar ()


# listMember ()
#
sub listMember {
    my ($item, @list) = @_;
    my $element;
    foreach $element (@list) {
        return $item  if ($item eq $element);
    }
    return undef;
}


# iterElimSortIndexList ()
#
sub iterElimSortIndexList {
    my ($val1_list, $val2_list, $fraction, $direction) = @_;
    my $index_list        = +[];
    my $local_index_list  = +[];
    my $local_val_list    = +[];
    my $local_sorted_list = +[];
    my ($index, $i, $j);

    my $sorted_val1_list = &insertSortIndexList ($val1_list, $direction);
    for ($i=0; $i <= $#{$sorted_val1_list}; ++$i) {
	$index_list->[$i] = $sorted_val1_list->[$i];
    }

    my $done = undef;
    my $toggle = 2;
    $cut = int ($#{$index_list} * $fraction);
    $last_cut = $#{$index_list};
    while ($cut > 0) {
	# sort the right half ("discards")
	$local_index_list = +[];
	$local_val_list   = +[];
	for ($j=0; $cut+$j+1 <= $last_cut; ++$j) {
	    $index                  = $index_list->[$cut+$j+1];
	    $local_index_list->[$j] = $index;
	    $local_val_list->[$j]   = ($toggle == 1) ? $val1_list->[$index]
		                                     : $val2_list->[$index];
	}
	$local_sorted_index_list = &insertSortIndexList ($local_val_list, $direction);
	for ($j=0; $cut+$j+1 <= $last_cut; ++$j) {
	    $local_index = $local_sorted_index_list->[$j];
	    $index_list->[$cut+$j+1] = $local_index_list->[$local_index];
	}

	# sort the left half ("keeps")
	$local_index_list = +[];
	$local_val_list   = +[];
	for ($j=0; $j <= $cut; ++$j) {
	    $index                  = $index_list->[$j];
	    $local_index_list->[$j] = $index;
	    $local_val_list->[$j]   = ($toggle == 1) ? $val1_list->[$index]
		                                     : $val2_list->[$index];
	}
	$local_sorted_index_list = &insertSortIndexList ($local_val_list, $direction);
	for ($j=0; $j <= $cut; ++$j) {
	    $local_index = $local_sorted_index_list->[$j];
	    $index_list->[$j] = $local_index_list->[$local_index];
	}
	
	# update cut and toggle
	$toggle = ($toggle == 1) ? 2 : 1;
	$last_cut = $cut;
	$cut = int ($last_cut * $fraction);
    }

    return $index_list;
}

# insertSortIndexList ()
#
sub insertSortIndexList {
    my ($val_list, $direction) = @_;
    my $index_list = +[];
    my ($index, $val, $i, $i2, $assigned);

    $index_list->[0] = 0;
    for ($index=1; $index <= $#{$val_list}; ++$index) {
        $assigned = undef;
        $val = $val_list->[$index];
        for ($i=0; $i <= $#{$index_list}; ++$i) {
            if ($direction eq 'decreasing') {
                if ($val > $val_list->[$index_list->[$i]]) {
                    for ($i2=$#{$index_list}; $i2 >= $i; --$i2) {
                        $index_list->[$i2+1] = $index_list->[$i2];
                    }
                    $index_list->[$i] = $index;
                    $assigned = 'true';
                    last;
                }
            }
            else {
                if ($val < $val_list->[$index_list->[$i]]) {
                    for ($i2=$#{$index_list}; $i2 >= $i; --$i2) {
                        $index_list->[$i2+1] = $index_list->[$i2];
                    }
                    $index_list->[$i] = $index;
                    $assigned = 'true';
                    last;
                }
            }
        }
        $index_list->[$#{$index_list}+1] = $index  if (! $assigned);
    }
    return $index_list;
}

# readFiles
#
sub readFiles {
    my ($dir, $fullpath_flag) = @_;
    my $inode;
    my @inodes = ();
    my @files = ();
    
    opendir (DIR, $dir);
    @inodes = sort readdir (DIR);
    closedir (DIR);
    foreach $inode (@inodes) {
	next if (! -f "$dir/$inode");
	next if ($inode =~ /^\./);
	push (@files, ($fullpath_flag) ? "$dir/$inode" : "$inode");
    }
    return @files;
}

# createDir
#
sub createDir {
    my $dir = shift;
    if (! -d $dir && (system (qq{mkdir -p $dir}) != 0)) {
	print STDERR "$0: unable to mkdir -p $dir\n";
	exit -2;
    }
    return $dir;
}

# copyFile
#
sub copyFile {
    my ($src, $dst) = @_;
    if (-f $src) {
	if (system (qq{cp $src $dst}) != 0) {
	    print STDERR "$0: unable to cp $src $dst\n";
	    exit -2;
	}
    } else {
	print STDERR "$0: file not found: '$src'\n";
    }
    return $dst;
}

# zip
#
sub zip {
    my $file = shift;
    if ($file =~ /^\.Z/ || $file =~ /\.gz/) {
	&abort ("already a zipped file $file");
    }
    if (-s $file) {
	if (system (qq{gzip -9f $file}) != 0) {
	    &abort ("unable to gzip -9f $file");
	}
    } elsif (-f $file) {
	&abort ("file empty: '$file'");
    } else {
	&abort ("file not found: '$file'");
    }
    $file .= ".gz";
    return $file;
}

# unzip
#
sub unzip {
    my $file = shift;
    if ($file !~ /^\.Z/ && $file !~ /\.gz/) {
	&abort ("not a zipped file $file");
    }
    if (-f $file) {
	if (system (qq{gzip -d $file}) != 0) {
	    &abort ("unable to gzip -d $file");
	}
    } else {
	&abort ("file not found: '$file'");
    }
    $file =~ s/\.Z$|\.gz$//;
    if (! -s $file) {
	&abort ("file empty: '$file'");
    }
    return $file;
}

# remove
#
sub remove {
    my $inode = shift;
    if (-e $inode) {
	if (system (qq{rm -rf $inode}) != 0) {
	    print STDERR "$0: unable to rm -rf $inode\n";
	    exit -2;
	}
    } else {
	print STDERR "$0: inode not found: '$inode'\n";
    }
    return $inode;
}
     
# runCmd
#
sub runCmd {
    my ($cmd, $nodie, $silent) = @_;
    my $ret;
    my $date = `date +'%Y-%m-%d_%T'`;  $date = &chompEnds ($date);
    print "[$date][RUN][$0] $cmd\n" if ($debug && ! $silent);
    $ret = system ($cmd);
    #$ret = ($?>>8)-256;
    $ret = ($?>>8);
    if ($ret != 0) {
	$ret -= 256  if ($ret >= 128);
	$date = `date +'%Y-%m-%d_%T'`;  $date = &chompEnds ($date);
	print STDERR ("[$date][FAILURE:$ret][$0] $cmd\n");
	if ($nodie) {
	    return $ret;
	} else {
	    exit $ret;
	}
    }
    return 0;
}

# logMsg()
#
sub logMsg {
    my ($msg, $logfile) = @_;
    my $date = `date +'%Y-%m-%d_%T'`;  $date = &chompEnds ($date);

    if ($logfile) {
        open (LOGFILE, ">".$logfile);
        select (LOGFILE);
    }
    else {
	select (STDOUT);
    }
    print "[$date][LOG][$0] $msg\n";
    if ($logfile) {
        close (LOGFILE);
    }
    select (STDOUT);

    return 'true';
}

# checkExist()
#
sub checkExist {
    my ($type, $path) = @_;
    if ($type eq 'd') {
	if (! -d $path) { 
            &alert ("dirnotfound: $path");
            exit -3;
	}
    }
    elsif ($type eq 'f') {
	if (! -f $path) {
            &alert ("filenotfound: $path");
            exit -3;
	}
	elsif (! -s $path) {
            &alert ("emptyfile: $path");
            exit -3;
	}
    }
}

# alert()
#
sub alert {
    my $msg = shift;
    my $date = `date +'%Y-%m-%d_%T'`;  $date = &chompEnds ($date);
    print STDERR "[$date][ALERT][$0] $msg\n";
    return;
}

# abort()
#
sub abort {
    my $msg = shift;
    my $date = `date +'%Y-%m-%d_%T'`;  $date = &chompEnds ($date);
    print STDERR "[$date][ABORT][$0] $msg\n";
    exit -2;
}

# writeBufToFile()
#
sub writeBufToFile {
    my ($file, $bufptr) = @_;
    if (! open (FILE, '>'.$file)) {
	&abort ("unable to open file $file for writing");
    }
    print FILE join ("\n", @{$bufptr}), "\n";
    close (FILE);
    return;
}

# fileBufString()
#
sub fileBufString {
    my $file = shift;
    my $oldsep = $/;
    undef $/;
    if ($file =~ /\.gz$|\.Z$/) {
	if (! open (FILE, "gzip -dc $file |")) {
	    &abort ("unable to open file $file for gzip -dc");
	}
    }
    elsif (! open (FILE, $file)) {
	&abort ("unable to open file $file for reading");
    }
    my $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    return $buf;
}

# fileBufArray()
#
sub fileBufArray {
    my $file = shift;
    my $oldsep = $/;
    undef $/;
    if ($file =~ /\.gz$|\.Z$/) {
	if (! open (FILE, "gzip -dc $file |")) {
	    &abort ("unable to open file $file for gzip -dc");
	}
    }
    elsif (! open (FILE, $file)) {
	&abort ("unable to open file $file for reading");
    }
    my $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    @buf = split (/$oldsep/, $buf);
    pop (@buf)  if ($buf[$#buf] eq '');
    return @buf;
}

# bigFileBufArray()
#
sub bigFileBufArray {
    my $file = shift;
    my $buf = +[];
    if ($file =~ /\.gz$|\.Z$/) {
        if (! open (FILE, "gzip -dc $file |")) {
            &abort ("unable to open file $file for gzip -dc");
        }
    }
    elsif (! open (FILE, $file)) {
        &abort ("unable to open file $file for reading");
    }
    while (<FILE>) {
        chomp;
        push (@$buf, $_);
    }
    close (FILE);
    return $buf;
}     

###############################################################################
# end
1;                                                     # in case it's a package
###############################################################################
