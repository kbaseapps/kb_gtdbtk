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
my $genome_type     = $opts{type};
my $gff_in_dir      = $opts{gffindir};
my $fna_in_dir      = $opts{fnaindir};
my $import_dir      = $opts{importdir};
my $domain          = $opts{domain};
my $release_version = $opts{releaseversion};
my $workspace_name  = $opts{workspacename};


# sdk cfg file
$sdk_cfg_file = $import_dir.'/sdk.cfg';
&checkExist ('f', $sdk_cfg_file);

# image dir
$image_dir = $import_dir."/images/$set_name";
if (-d $image_dir) {
    &abort ("already have $image_dir");
}
&createDir ($image_dir);
&copyFile ($sdk_cfg_file, $image_dir);

# params file
$params_dir = $import_dir."/params";
&createDir ($params_dir);
$params_file = $params_dir.'/'.$set_name.".json";
if (-s $params_file) {
    &abort ("already have $params_file");
}

# script file  
$scripts_dir = $import_dir."/scripts";
&createDir ($scripts_dir);
$script_file = $scripts_dir.'/'.$set_name.".sh";
if (-s $script_file) {
    &abort ("already have $script_file");
}    

# staging dir
$staging_dir = $import_dir."/staging".'/'.$set_name;
if (-d $staging_dir) {
    &abort ("already have $staging_dir");
}
&createDir ($staging_dir);

# clean up paths
$image_dir   =~ s/\/\//\//g;
$image_dir   =~ s/\/$//;
$params_file =~ s/\/\//\//g;
$script_file =~ s/\/\//\//g;
$staging_dir =~ s/\/\//\//g;
$staging_dir =~ s/\/$//;

###############################################################################
# main
###############################################################################

# put data into staging dir
$genome_cnt = 0;
foreach $genomeID (&fileBufArray($genomeIDs_file)) {
    $genome_cnt += 1;
    $db_prefix = 'RS';
    $dbid_prefix = 'GCF';
    if ($genomeID =~ /GCA/) {
	$db_prefix = 'GB';
	$dbid_prefix = 'GCA';
    }
    $genomeID =~ /(\d{3})(\d{3})(\d{3})/;
    $id_3 = $1;
    $id_6 = $2;
    $id_9 = $3;
    $domain_prefix = 'archaea';
    if ($domain =~ /^B/i) {
	$domain_prefix = 'bacteria';
    }

    if ($object_prefix !~ /\-$/ && $object_prefix !~ /\_$/) {
	$object_prefix .= '-';
    }
    $object_basename = $object_prefix.$genomeID;
    
    # src gff file
    $src_gff_file = join('/', $gff_in_dir, $db_prefix, $domain_prefix, $db_prefix.'_'.$genomeID.'.gff');

    # src fna file
    $src_fna_file = join('/', $fna_in_dir, $dbid_prefix, $id_3, $id_6, $id_9, $genomeID.'_genomic.fna.gz');

    # dst gff file
    $dst_gff_file = join ('/', $staging_dir, $object_basename.'.gff');

    # dst fna file
    $dst_fna_file = join ('/', $staging_dir, $object_basename.'.fna');
    
    print ("STAGING GFF for $genome_cnt $genomeID\n");
    &copyFile ($src_gff_file, $dst_gff_file);
    print ("STAGING FNA for $genome_cnt $genomeID\n");
    &copyFile ($src_fna_file, $dst_fna_file.'.gz');
    gunzip $dst_fna_file.'.gz' => $dst_fna_file
	or die "gunzip failed for file $dst_fna_file.gz: $GunzipError\n";
    unlink $dst_fna_file.'.gz';
}


# create params file
@params_buf = ();
$control_vocab_type = 'draft isolate';
if ($genome_type =~ /^M/i) {
    $control_vocab_type = 'mag';
} elsif ($genome_type =~ /^S/i) {
    $control_vocab_type = 'sag';
}
push (@params_buf, '{');
push (@params_buf, '  "workspace_name": "'.$workspace_name.'",');
push (@params_buf, '  "staging_subdir": "'.'/'.'",');
push (@params_buf, '  "genome_set_name": "'.$set_name.'.GenomeSet'.'",');
push (@params_buf, '  "genome_type": "'.$control_vocab_type.'",');
push (@params_buf, '  "source": "'.'GTDB'.'",');
push (@params_buf, '  "taxon_wsname": "'.''.'",');
push (@params_buf, '  "taxon_id": "'.''.'",');
push (@params_buf, '  "release": "'.$release_version.'",');
push (@params_buf, '  "genetic_code": '.'11'.',');
push (@params_buf, '  "generate_missing_genes": 1');
push (@params_buf, '}');    
&writeBufToFile ($params_file, \@params_buf);


# create script file
@script_buf = ();
push (@script_buf, "#!/bin/bash\n\n");
push (@script_buf, qq{kb-sdk run --input $params_file --mount-points $staging_dir:/data/bulk/dylan --sdk-home $image_dir kb_uploadmethods.batch_import_genomes_from_staging
});
&writeBufToFile ($script_file, \@script_buf);
chmod 0755, $script_file;


# output
if ($out_file) {
    open (OUT, '>'.$out_file);
    select (OUT);
}
print join ("\n", @out)."\n"  if (@out);
if ($out_file) {
    close (OUT);
    select (STDOUT);
}

# done
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
\t\t-genomeidsfile  <genomeIDs_file>
\t\t-setname        <set_name>        (e.g. 'GTDB_Arc-Iso_RS')
\t\t-objectprefix   <object_prefix>   (e.g. 'GTDB_Arc-')
\t\t-type           <genome_type>     (e.g. 'Isolate'/'MAG'/'SAG')
\t\t-gffindir       <gff_in_dir>
\t\t-fnaindir       <fna_in_dir>
\t\t-importdir      <import_dir>
\t\t-releaseversion <release_version> (e.g. 'GTDB_r207')
\t\t-domain         <domain>          (e.g. 'archaea'/'bacteria')
\t\t-workspace_name <workspace_name>  (e.g. 'dylan:narrative_1653154088731')
};

    # Get args
    my %opts = ();
    &GetOptions (\%opts, 
		 "genomeidsfile=s",
		 "setname=s",
		 "objectprefix=s",
		 "type=s",
		 "gffindir=s",
		 "fnaindir=s",
		 "importdir=s",
		 "releaseversion=s",
		 "domain=s",
		 "workspacename=s");

    # Check for legal invocation
    if (! defined $opts{genomeidsfile} ||
	! defined $opts{setname} ||
	! defined $opts{objectprefix} ||
	! defined $opts{type} ||
	! defined $opts{gffindir} ||
	! defined $opts{fnaindir} ||
	! defined $opts{importdir} ||
	! defined $opts{releaseversion} ||
	! defined $opts{domain} ||
	! defined $opts{workspacename}
	) {
        print STDERR "$usage\n";
        exit -1;
    }
    &checkExist ('f', $opts{genomeidsfile});
    &checkExist ('d', $opts{gffindir});
    &checkExist ('d', $opts{fnaindir});
    &checkExist ('d', $opts{importdir});

    # ensure absolute paths
    if ($opts{gffindir} !~ /^\//) {
	&abort ("paths must be absolute.  -gffindir is not");
    }
    if ($opts{fnaindir} !~ /^\//) {
	&abort ("paths must be absolute.  -fnaindir is not");
    }
    if ($opts{importdir} !~ /^\//) {
	&abort ("paths must be absolute.  -importdir is not");
    }
    
    if ($opts{type} !~ /^I/i &&
	$opts{type} !~ /^M/i &&
	$opts{type} !~ /^S/i) {
	&abort ("-type must be one of 'Isolate'/'MAG'/'SAG'");
    }
    if ($opts{domain} !~ /^A/i &&
	$opts{domain} !~ /^B/i) {
	&abort ("-domain must be one of 'archaea'/'bacteria'");
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
