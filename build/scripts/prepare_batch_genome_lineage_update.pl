#!/usr/bin/perl
##
## Copyright 2020-2023 Dylan Chivian
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
my $metadata_file   = $opts{metadatafile};
my $update_dir      = $opts{updatedir};
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

# lineage dir
$lineage_dir = $update_dir."/lineage/$set_name";
if (-d $lineage_dir) {
    &abort ("already have $lineage_dir");
}
&createDir ($lineage_dir);


# clean up paths
$image_dir   =~ s/\/\//\//g;
$image_dir   =~ s/\/$//;
$params_file =~ s/\/\//\//g;
$script_file =~ s/\/\//\//g;
$lineage_dir  =~ s/\/\//\//g;
$lineage_dir  =~ s/\/$//;

###############################################################################
# main
###############################################################################


# get target genome ids
print ("READING target genome IDs...\n");
%candidate_target_genomeIDs = ();
foreach $genomeID (&fileBufArray($genomeIDs_file)) {
    $genomeID =~ s/^GTDB_//;
    $genomeID =~ s/^Arc_//;
    $genomeID =~ s/^Bac_//;
    $genomeID =~ s/^GB_//;
    $genomeID =~ s/^RS_//;
    $candidate_target_genomeIDs{$genomeID} = 1;

    #print ("GENOME_ID: '$genomeID'\n");
}

# read metadata to determine which target IDs are present (may not be for deprecated genomes)
print ("READING metadata for valid targets...\n");
@genomeIDs_buf = ();
%target_genomeIDs = ();
$GENOME_ID_I            = 0;
foreach $metadata_line (@{&bigFileBufArray($metadata_file)}) {
    @metadata = split ("\t", $metadata_line);
    $genomeID = $metadata[$GENOME_ID_I];
    $genomeID =~ s/^GB_//;
    $genomeID =~ s/^RS_//;

    #print ("GENOME_ID: '$genomeID'\n");

    next  if (! $candidate_target_genomeIDs{$genomeID});
    $target_genomeIDs{$genomeID} = 1;
    push (@genomeIDs_buf, $genomeID);
}
if ($#genomeIDs_buf == -1) {
    print ("No IDs in new metadata");
    exit (0);
}

# create field bufs that don't read metadata
print ("STORING direct vals...\n");
@release_buf = ();

foreach $genomeID (sort keys %target_genomeIDs) {

    #print ("GENOMEID: '$genomeID'\n");
    
    # release
    push (@release_buf, join("\t", $genomeID, $release_version));
}


# re-read metadata and fill field bufs
print ("RE-READING metadata...\n");
@tax_hierarchy_buf      = ();
$GENOME_ID_I            = 0;
$GTDB_TAX_I             = 16;
foreach $metadata_line (@{&bigFileBufArray($metadata_file)}) {
    @metadata = split ("\t", $metadata_line);
    $genomeID = $metadata[$GENOME_ID_I];
    $genomeID =~ s/^GB_//;
    $genomeID =~ s/^RS_//;

    next  if (! $target_genomeIDs{$genomeID});
    
    # GTDB taxonomy hierarchy
    $tax_hierarchy_val = $metadata[$GTDB_TAX_I];
    push (@tax_hierarchy_buf, join("\t", $genomeID, $tax_hierarchy_val));
}


# save field files
print ("SAVING field val files...\n");
$targets_list_file = join ('/', $lineage_dir, 'targets'.'.list');
$release_file = join ('/', $lineage_dir, 'release'.'.map');
$tax_hierarchy_file = join ('/', $lineage_dir, 'tax_hierarchy'.'.map');
#&copyFile ($genomeIDs_file, $targets_list_file);
&writeBufToFile ($targets_list_file, \@genomeIDs_buf);
&writeBufToFile ($release_file, \@release_buf);
&writeBufToFile ($tax_hierarchy_file, \@tax_hierarchy_buf);


# create params file
print ("SAVING params file...\n");
@params_buf = ();
push (@params_buf, '{');
push (@params_buf, '  "workspace_name": "'.$workspace_name.'",');
push (@params_buf, '  "target_list_file": "'.'/lineage/targets.list'.'",');
push (@params_buf, '  "release_file": "'.'/lineage/release.map'.'",');
if ($opts{'deleteoldtaxonassignments'}) {
    push (@params_buf, '  "taxonomy_hierarchy_file": "'.'/lineage/tax_hierarchy.map'.'",');
    push (@params_buf, '  "delete_old_taxon_assignments": 1');
} else {
    push (@params_buf, '  "taxonomy_hierarchy_file": "'.'/lineage/tax_hierarchy.map'.'"');
}
push (@params_buf, '}');    
&writeBufToFile ($params_file, \@params_buf);


# create script file
print ("SAVING script file...\n");
@script_buf = ();
push (@script_buf, "#!/bin/bash\n\n");
push (@script_buf, qq{kb-sdk run -t beta --input $params_file --mount-points $lineage_dir:/lineage --sdk-home $image_dir kb_ObjectUtilities.KButil_update_genome_lineage_from_files
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
\t\t -genomeidsfile   <genomeIDs_file>
\t\t -setname         <set_name>        (e.g. 'GTDB_Arc-Iso-RS')
\t\t -metadatafile    <metadata_file>   (e.g. 'ar53_metadata_r207.tsv')
\t\t -releaseversion  <release_version> (e.g. 'GTDB_r207')
\t\t -updatedir       <update_dir>
\t\t -workspace_name  <workspace_name>  (e.g. 'dylan:narrative_1653154088731')
\t\t[-deleteoldtaxonassignments  <TRUE/FALSE>]  (def: FALSE)
};

    # Get args
    my %opts = ();
    &GetOptions (\%opts, 
		 "genomeidsfile=s",
		 "setname=s",
		 "metadatafile=s",
		 "updatedir=s",
                 "releaseversion=s",
		 "deleteoldtaxonassignments=s",
		 "workspacename=s");

    # Check for legal invocation
    if (! defined $opts{genomeidsfile} ||
	! defined $opts{setname} ||
	! defined $opts{metadatafile} ||
	! defined $opts{updatedir} ||
        ! defined $opts{releaseversion} ||
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

    # defaults
    if (defined $opts{'deleteoldtaxonassignments'} &&
	$opts{'deleteoldtaxonassignments'} =~ /^f/i) {
	$opts{'deleteoldtaxonassignments'} = undef
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
