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
# imports
###############################################################################

use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

###############################################################################
# conf
###############################################################################

# general
$| = 1;                                              # disable stdout buffering
$debug = 1;                                             # chatter while running

###############################################################################
# base init
###############################################################################

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

###############################################################################
# init
###############################################################################

# argv
my %opts = &getCommandLineOptions ();
my $genomeIDs_file  = $opts{genomeidsfile};
my $set_name        = $opts{setname};
my $object_prefix   = $opts{objectprefix};
my $metadata_file   = $opts{metadatafile};
my $update_dir      = $opts{updatedir};
my $domain          = $opts{domain};
my $release_version = $opts{releaseversion};
my $workspace_name  = $opts{workspacename};


# sdk cfg file
$sdk_cfg_file = $update_dir.'/sdk.cfg';
&checkExist ('f', $sdk_cfg_file);

# image dir
$image_dir = $update_dir."/images/$set_name";
if (-d $image_dir) {
    &abort ("already have $image_dir");
}
&createDir ($image_dir);
&copyFile ($sdk_cfg_file, $image_dir);

# params file
$params_dir = $update_dir."/params";
&createDir ($params_dir);
$params_file = $params_dir.'/'.$set_name.".json";
if (-s $params_file) {
    &abort ("already have $params_file");
}

# script file  
$scripts_dir = $update_dir."/scripts";
&createDir ($scripts_dir);
$script_file = $scripts_dir.'/'.$set_name.".sh";
if (-s $script_file) {
    &abort ("already have $script_file");
}    

# fields dir
$fields_dir = $update_dir."/fields/$set_name";
if (-d $fields_dir) {
    &abort ("already have $fields_dir");
}
&createDir ($fields_dir);


# clean up paths
$image_dir   =~ s/\/\//\//g;
$image_dir   =~ s/\/$//;
$params_file =~ s/\/\//\//g;
$script_file =~ s/\/\//\//g;
$fields_dir  =~ s/\/\//\//g;
$fields_dir  =~ s/\/$//;

# clean up args
$object_prefix =~ s/\-$//;

###############################################################################
# main
###############################################################################


# get target genome ids
print ("READING target genome IDs...\n");
%target_genomeIDs = ();
foreach $genomeID (&fileBufArray($genomeIDs_file)) {
    $genomeID =~ s/^GTDB_//;
    $genomeID =~ s/^Arc_//;
    $genomeID =~ s/^Bac_//;
    $genomeID =~ s/^GB_//;
    $genomeID =~ s/^RS_//;
    $target_genomeIDs{$genomeID} = 1;
}

# create field bufs that don't read metadata
print ("STORING direct vals...\n");
@obj_newname_buf = ();
@domain_buf = ();
@source_buf = ();
@release_buf = ();
if ($domain =~ /^a/) {
    $domain_val = 'Archaea';
} else {
    $domain_val = 'Bacteria';
}
foreach $genomeID (sort keys %target_genomeIDs) {

    #print ("GENOMEID: '$genomeID'\n");
    
    # obj newname
    $obj_newname_val = $object_prefix.'-'.$genomeID.'.'.'Genome';
    push (@obj_newname_buf, join("\t", $genomeID, $obj_newname_val));

    # domain
    push (@domain_buf, join("\t", $genomeID, $domain_val));

    # source db
    if ($genomeID =~ /^GCA/) {
	$source_val = 'GenBank';
    } elsif ($genomeID =~ /^GCF/) {
	$source_val = 'RefSeq';
    } else {
	&abort ("unknown genome source for genome ID $genomeID");
    }
    push (@source_buf, join("\t", $genomeID, $source_val));

    # release
    push (@release_buf, join("\t", $genomeID, $release_version));
}


# read metadata and fill field bufs
print ("READING metadata...\n");
@species_name_buf       = ();
@genome_type_buf        = ();
@ncbi_tax_id_buf        = ();
@tax_hierarchy_buf      = ();
@genome_qual_scores_buf = ();
$GENOME_ID_I            = 0;
$CHECKM_COMPLETENESS_I  = 2;
$CHECKM_CONTAMINATION_I = 3;
$GTDB_TAX_I             = 16;
$ASSEMBLY_STATUS_I      = 45; 
$GENOME_TYPE_I          = 55;
$NCBI_SPECIES_NAME_I    = 62;
$NCBI_TAXID_I           = 77;
foreach $metadata_line (@{&bigFileBufArray($metadata_file)}) {
    @metadata = split ("\t", $metadata_line);
    $genomeID = $metadata[$GENOME_ID_I];
    $genomeID =~ s/^GB_//;
    $genomeID =~ s/^RS_//;

    next  if (! $target_genomeIDs{$genomeID});
    
    # soecies name
    $species_name_val = $metadata[$NCBI_SPECIES_NAME_I];
    push (@species_name_buf, join ("\t", $genomeID, $species_name_val));

    # genome type
    $genome_type_string = $metadata[$GENOME_TYPE_I];
    if ($genome_type_string eq 'derived from metagenome' ||
	$genome_type_string eq 'derived from environmental sample') {
	$genome_type_val = 'mag';
    } elsif ($genome_type_string eq 'derived from single cell') {
	$genome_type_val = 'sag';
    } elsif ($genome_type_string eq 'none') {
	if ($metadata[$ASSEMBLY_STATUS_I] eq 'Complete Genome') {
	    $genome_type_val = 'finished isolate';
	} else {
	    $genome_type_val = 'draft isolate';
	}
    } else {
	&abort ("unknown genome type '$genome_type_string' for genome ID $genomeID");
    }
    # fix misclassifications
    if ($genomeID eq 'GCF_000018425.1') {  # Candidatus Desulforudis audaxviator MP104C
	$genome_type_val = 'mag';
    }
    push (@genome_type_buf, join("\t", $genomeID, $genome_type_val));

    # NCBI tax id
    $ncbi_tax_id_val = $metadata[$NCBI_TAXID_I];
    push (@ncbi_tax_id_buf, join("\t", $genomeID, $ncbi_tax_id_val));
    
    # GTDB taxonomy hierarchy
    $tax_hierarchy_val = $metadata[$GTDB_TAX_I];
    push (@tax_hierarchy_buf, join("\t", $genomeID, $tax_hierarchy_val));
    
    # CheckM scores
    $checkm_completeness  = $metadata[$CHECKM_COMPLETENESS_I];
    $checkm_contamination = $metadata[$CHECKM_CONTAMINATION_I];
    $completeness_info = join(',', 
			      join('=', 'method', 'CheckM'),
			      join('=', 'method_version', '1.1.6'),
			      join('=', 'score', $checkm_completeness),
			      join('=', 'score_interpretation', 'percent_completeness'),
			      join('=', 'timestamp', '2022-04-08_00:00:00'));
    $contamination_info = join(',', 
			       join('=', 'method', 'CheckM'),
			       join('=', 'method_version', '1.1.6'),
			       join('=', 'score', $checkm_contamination),
			       join('=', 'score_interpretation', 'percent_contamination'),
			       join('=', 'timestamp', '2022-04-08_00:00:00'));
    $genome_qual_scores_val = join(';', $completeness_info, $contamination_info);
    push (@genome_qual_scores_buf, join("\t", $genomeID, $genome_qual_scores_val));
}


# save field files
print ("SAVING field val files...\n");
$targets_list_file = join ('/', $fields_dir, 'targets'.'.list');
$obj_newname_file = join ('/', $fields_dir, 'obj_newname'.'.map');
$domain_file = join ('/', $fields_dir, 'domain'.'.map');
$source_file = join ('/', $fields_dir, 'source'.'.map');
$release_file = join ('/', $fields_dir, 'release'.'.map');
$species_name_file = join ('/', $fields_dir, 'species_name'.'.map');
$genome_type_file = join ('/', $fields_dir, 'genome_type'.'.map');
$ncbi_tax_id_file = join ('/', $fields_dir, 'ncbi_tax_id'.'.map');
$tax_hierarchy_file = join ('/', $fields_dir, 'tax_hierarchy'.'.map');
$genome_qual_scores_file = join ('/', $fields_dir, 'genome_qual_scores'.'.map');
&copyFile ($genomeIDs_file, $targets_list_file);
&writeBufToFile ($obj_newname_file, \@obj_newname_buf);
&writeBufToFile ($domain_file, \@domain_buf);
&writeBufToFile ($source_file, \@source_buf);
&writeBufToFile ($release_file, \@release_buf);
&writeBufToFile ($species_name_file, \@species_name_buf);
&writeBufToFile ($genome_type_file, \@genome_type_buf);
&writeBufToFile ($ncbi_tax_id_file, \@ncbi_tax_id_buf);
&writeBufToFile ($tax_hierarchy_file, \@tax_hierarchy_buf);
&writeBufToFile ($genome_qual_scores_file, \@genome_qual_scores_buf);


# create params file
print ("SAVING params file...\n");
@params_buf = ();
push (@params_buf, '{');
push (@params_buf, '  "workspace_name": "'.$workspace_name.'",');
push (@params_buf, '  "target_list_file": "'.'/fields/targets.list'.'",');
push (@params_buf, '  "object_newname_file": "'.'/fields/obj_newname.map'.'",');
push (@params_buf, '  "species_name_file": "'.'/fields/species_name.map'.'",');
push (@params_buf, '  "source_file": "'.'/fields/source.map'.'",');
push (@params_buf, '  "domain_file": "'.'/fields/domain.map'.'",');
push (@params_buf, '  "genome_type_file": "'.'/fields/genome_type.map'.'",');
push (@params_buf, '  "release_file": "'.'/fields/release.map'.'",');
push (@params_buf, '  "taxonomy_hierarchy_file": "'.'/fields/tax_hierarchy.map'.'",');
push (@params_buf, '  "taxonomy_ncbi_id_file": "'.'/fields/ncbi_tax_id.map'.'",');
push (@params_buf, '  "genome_qual_scores_file": "'.'/fields/genome_qual_scores.map'.'",');
push (@params_buf, '  "keep_spoofed_mRNAs": "'.'0'.'"');
push (@params_buf, '}');    
&writeBufToFile ($params_file, \@params_buf);


# create script file
print ("SAVING script file...\n");
@script_buf = ();
push (@script_buf, "#!/bin/bash\n\n");
push (@script_buf, qq{kb-sdk run -t beta --input $params_file --mount-points $fields_dir:/fields --sdk-home $image_dir kb_ObjectUtilities.KButil_update_genome_fields_from_files
});
&writeBufToFile ($script_file, \@script_buf);
chmod 0755, $script_file;


# done
print ("DONE\n");
exit 0;

###############################################################################
# subs
###############################################################################

# getCommandLineOptions()
#
#  rets: \%opts  pointer to hash of kv pairs of command line options
#
sub getCommandLineOptions {
    use Getopt::Long;
    my $usage = qq{usage: $0 
\t\t-genomeidsfile   <genomeIDs_file>
\t\t-setname         <set_name>        (e.g. 'GTDB_Arc-Iso_RS')
\t\t-objectprefix    <object_prefix>   (e.g. 'GTDB_Arc-')
\t\t-metadatafile    <metadata_file>   (e.g. 'ar53_metadata_r207.tsv')
\t\t-releaseversion  <release_version> (e.g. 'GTDB_r207')
\t\t-domain          <domain>          (e.g. 'archaea'/'bacteria')
\t\t-updatedir       <update_dir>
\t\t-workspace_name  <workspace_name>  (e.g. 'dylan:narrative_1653154088731')
};

    # Get args
    my %opts = ();
    &GetOptions (\%opts, 
		 "genomeidsfile=s",
		 "setname=s",
		 "objectprefix=s",
		 "metadatafile=s",
		 "updatedir=s",
                 "releaseversion=s",
                 "domain=s",		 
		 "workspacename=s");

    # Check for legal invocation
    if (! defined $opts{genomeidsfile} ||
	! defined $opts{setname} ||
	! defined $opts{objectprefix} ||
	! defined $opts{metadatafile} ||
	! defined $opts{updatedir} ||
        ! defined $opts{releaseversion} ||
        ! defined $opts{domain} ||	
	! defined $opts{workspacename}
	) {
        print STDERR "$usage\n";
        exit -1;
    }
    &checkExist ('f', $opts{genomeidsfile});
    &checkExist ('f', $opts{metadatafile});
    &checkExist ('d', $opts{updatedir});

    # ensure absolute paths
    if ($opts{updatedir} !~ /^\//) {
	&abort ("paths must be absolute.  -updatedir is not");
    }
    
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
