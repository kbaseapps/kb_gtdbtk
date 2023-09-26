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
my $update_dir           = $opts{updatedir};
my $set_name             = $opts{setname};
my $eggnog_annot_file    = $opts{eggnogannotfile};
my $novelfams_annot_file = $opts{novelfamsannotfile};

# sdk cfg file
$sdk_cfg_file = $update_dir.'/sdk.cfg';
&checkExist ('f', $sdk_cfg_file);

# image dir
$image_dir = $update_dir."/images/$set_name";
if (-d $image_dir) {
    &alert ("already have $image_dir");
}
&createDir ($image_dir);
&copyFile ($sdk_cfg_file, $image_dir);

# params file
$params_dir = $update_dir."/params";
&createDir ($params_dir);
$params_file = $params_dir.'/'.$set_name.".json";
if (-s $params_file) {
    &alert ("already have $params_file");
}

# script file  
$scripts_dir = $update_dir."/scripts";
&createDir ($scripts_dir);
$script_file = $scripts_dir.'/'.$set_name.".sh";
if (-s $script_file) {
    &alert ("already have $script_file");
}    

# features dir
$features_dir = $update_dir."/features";
&createDir ($features_dir);
$features_update_file = $set_name.'.map';
$features_update_path = join ('/', $features_dir, $features_update_file);
if (! -f $features_update_path) {
    &alert ("already have $features_dir");
}


# clean up paths
$image_dir    =~ s/\/\//\//g;
$image_dir    =~ s/\/$//;
$params_file  =~ s/\/\//\//g;
$script_file  =~ s/\/\//\//g;
$features_dir =~ s/\/\//\//g;
$features_dir =~ s/\/$//;


###############################################################################
# config
###############################################################################

$genome_ref_feature_id_delim = '.f:';
$genome_id_feature_id_delim = '-';

###############################################################################
# main
###############################################################################

# read eggnog annotations and fill buf
$QUERY_ID_I = 0;
$SEED_ORTHOLOG_I = 1;
$EVALUE_I = 2;
$BITSCORE_I = 3;
$EGGNOG_OG_I = 4;
$MAX_ANNOT_LEVEL_I = 5;
$COG_CAT_I = 6;
$DESC_I = 7;
$PREF_NAME_I = 8;
$GO_I = 9;
$EC_I = 10;
$KEGG_KO_I = 11;
$KEGG_PATHWAY_I = 12;
$KEGG_MODULE_I = 13;
$KEGG_RXN_I = 14;
$KEGG_RCLASS_I = 15;
$BRITE_I = 16;
$KEGG_TC_I = 17;
$CAZY_I = 18;
$BIGG_RXN_I = 19;
$PFAMS_I = 20;

print ("READING eggnog annotations...\n");
$f_up = +{};
foreach $eggnog_annot_line (@{&bigFileBufArray($eggnog_annot_file)}) {
    next if ($eggnog_annot_line =~ /^\#/);
    @eggnog_line = split ("\t", $eggnog_annot_line);

    $query_id = $eggnog_line[$QUERY_ID_I];
    if ($query_id =~ /$genome_ref_feature_id_delim/) {
	($genome_ref, $fid) = split (/$genome_ref_feature_id_delim/, $query_id);
	$genome_id = $genome_ref;
    } else {
	$query_id =~ /^([^$genome_id_feature_id_delim]+)$genome_id_feature_id_delim(.*)$/;
	$genome_id = $1;
	$fid = $2;
    }

    # instantiate genome
#    $f_up->{$genome_id} = +{}  if (! defined $f_up->{$genome_id});

    # instantiate feature
#    $f_up->{$genome_id}->{$fid} = +{}  if (! defined $f_up->{$genome_id}->{$fid});

    # instantiate aliases
#    $f_up->{$genome_id}->{$fid}->{'aliases'} = +{}  if (! defined $f_up->{$genome_id}->{$fid}->{'aliases'});

    # instantiate functions
#    $f_up->{$genome_id}->{$fid}->{'functions'} = +{}  if (! defined $f_up->{$genome_id}->{$fid}->{'functions'});
    
    # PREF NAME
    $F_I = $PREF_NAME_I;
    if ($eggnog_line[$F_I] && $eggnog_line[$F_I] ne '-' && $eggnog_line[$F_I] !~ /^\s*$/) {
	#print join ("\t", $genome_id, $fid, $eggnog_line[$F_I])."\n";  # DEBUG
#	$f_up->{$genome_id}->{$fid}->{'aliases'}->{'gene'} = +[]  if (! defined $f_up->{$genome_id}->{$fid}->{'aliases'}->{'gene'});
	foreach $val (split (/,/, $eggnog_line[$F_I])) {
	    push (@{$f_up->{$genome_id}->{$fid}->{'aliases'}->{'gene'}}, $val);
	}
    }

    # EGGNOG_OG
    $F_I = $EGGNOG_OG_I;
    if ($eggnog_line[$F_I] && $eggnog_line[$F_I] ne '-' && $eggnog_line[$F_I] !~ /^\s*$/) {
	#print join ("\t", $genome_id, $fid, $eggnog_line[$F_I])."\n";  # DEBUG
#	$f_up->{$genome_id}->{$fid}->{'functions'}->{'egg'} = +[]  if (! defined $f_up->{$genome_id}->{$fid}->{'functions'}->{'egg'});

	foreach $val (split (/,/, $eggnog_line[$F_I])) {
	    ($fam, $tax) = split (/\@/, $val);
	    push (@{$f_up->{$genome_id}->{$fid}->{'functions'}->{'egg'}}, $fam);
	}
    }
    
    # DESC
    $F_I = $DESC_I;
    if ($eggnog_line[$F_I] && $eggnog_line[$F_I] ne '-' && $eggnog_line[$F_I] !~ /^\s*$/) {
#	$f_up->{$genome_id}->{$fid}->{'functions'}->{'egg_desc'} = +[]  if (! defined $f_up->{$genome_id}->{$fid}->{'functions'}->{'egg_desc'});
	push (@{$f_up->{$genome_id}->{$fid}->{'functions'}->{'egg_desc'}}, $eggnog_line[$F_I]);
    }

    # EC
    $F_I = $EC_I;
    if ($eggnog_line[$F_I] && $eggnog_line[$F_I] ne '-' && $eggnog_line[$F_I] !~ /^\s*$/) {
#	$f_up->{$genome_id}->{$fid}->{'functions'}->{'EC'} = +[]  if (! defined $f_up->{$genome_id}->{$fid}->{'functions'}->{'EC'});

	foreach $ec (split (/,/, $eggnog_line[$F_I])) {
	    push (@{$f_up->{$genome_id}->{$fid}->{'functions'}->{'EC'}}, $ec);
	}
    }

    # KO
    $F_I = $KEGG_KO_I;
    if ($eggnog_line[$F_I] && $eggnog_line[$F_I] ne '-' && $eggnog_line[$F_I] !~ /^\s*$/) {
#	$f_up->{$genome_id}->{$fid}->{'functions'}->{'KO'} = +[]  if (! defined $f_up->{$genome_id}->{$fid}->{'functions'}->{'KO'});

	foreach $ko (split (/,/, $eggnog_line[$F_I])) {
	    $ko =~ s/^ko\://;
	    push (@{$f_up->{$genome_id}->{$fid}->{'functions'}->{'KO'}}, $ko);
	}
    }

    # CAZY
    $F_I = $CAZY_I;
    if ($eggnog_line[$F_I] && $eggnog_line[$F_I] ne '-' && $eggnog_line[$F_I] !~ /^\s*$/) {
#	$f_up->{$genome_id}->{$fid}->{'functions'}->{'CAZy'} = +[]  if (! defined $f_up->{$genome_id}->{$fid}->{'functions'}->{'CAZy'});

	foreach $cazy (split (/,/, $eggnog_line[$F_I])) {
	    push (@{$f_up->{$genome_id}->{$fid}->{'functions'}->{'CAZy'}}, $cazy);
	}
    }

    # PFAM
    $F_I = $PFAMS_I;
    if ($eggnog_line[$F_I] && $eggnog_line[$F_I] ne '-' && $eggnog_line[$F_I] !~ /^\s*$/) {
	#print join ("\t", $genome_id, $fid, $eggnog_line[$F_I])."\n";  # DEBUG
#	$f_up->{$genome_id}->{$fid}->{'functions'}->{'Pfam'} = +[]  if (! defined $f_up->{$genome_id}->{$fid}->{'functions'}->{'Pfam'});

	foreach $pfam (split (/,/, $eggnog_line[$F_I])) {
	    push (@{$f_up->{$genome_id}->{$fid}->{'functions'}->{'Pfam'}}, $pfam);
	    #print join ("\t", $genome_id, $fid, $pfam)."\n";  # DEBUG
	}
    }

    # Inference
    $f_up->{$genome_id}->{$fid}->{'inference'} = +{};
    if ($eggnog_line[$SEED_ORTHOLOG_I] && $eggnog_line[$SEED_ORTHOLOG_I] ne '-' && $eggnog_line[$SEED_ORTHOLOG_I] !~ /^\s*$/) {
	$seed_ortholog = $eggnog_line[$SEED_ORTHOLOG_I];
	$e_value = $eggnog_line[$EVALUE_I];
	$bit_score = $eggnog_line[$BITSCORE_I];

	$f_up->{$genome_id}->{$fid}->{'inference'}->{'eggnog5'} = +{};
	$f_up->{$genome_id}->{$fid}->{'inference'}->{'eggnog5'}->{'category'} = 'eggnog5-base';
	$f_up->{$genome_id}->{$fid}->{'inference'}->{'eggnog5'}->{'type'} = 'sequence homology';
	$f_up->{$genome_id}->{$fid}->{'inference'}->{'eggnog5'}->{'evidence'} = join(",",
										     'seed_ortholog='.$seed_ortholog,
										     'e_value='.$e_value,
										     'bit_score='.$bit_score);
    }
}

# add novel fams
#
$QUERY_ID_I = 0;
$TARGET_I = 1;
$EVALUE_I = 2;
$BITSCORE_I = 3;
$NOVEL_FAMS_I = 4;

if ($novelfams_annot_file) {
    print ("READING novel_fam annotations...\n");
    foreach $eggnog_annot_line (@{&bigFileBufArray($novelfams_annot_file)}) {
	next if ($eggnog_annot_line =~ /^\#/);
	@eggnog_line = split ("\t", $eggnog_annot_line);

	$query_id = $eggnog_line[$QUERY_ID_I];
	if ($query_id =~ /$genome_ref_feature_id_delim/) {
	    ($genome_ref, $fid) = split (/$genome_ref_feature_id_delim/, $query_id);
	    $genome_id = $genome_ref;
	} else {
	    $query_id =~ /^([^$genome_id_feature_id_delim]+)$genome_id_feature_id_delim(.*)$/;
	    $genome_id = $1;
	    $fid = $2;
	}

	# NOVEL FAMS
	$F_I = $NOVEL_FAMS_I;
	if ($eggnog_line[$F_I] && $eggnog_line[$F_I] ne '-' && $eggnog_line[$F_I] !~ /^\s*$/) {
#	    $f_up->{$genome_id}->{$fid}->{'functions'}->{'egg'} = +[]  if (! defined $f_up->{$genome_id}->{$fid}->{'functions'}->{'egg'});

	    foreach $val (split (/,/, $eggnog_line[$F_I])) {
		($fam, $tax) = split (/\@/, $val);
		push (@{$f_up->{$genome_id}->{$fid}->{'functions'}->{'egg'}}, $fam);
	    }
	}

	# Inference
	if (! defined $f_up->{$genome_id}->{$fid}->{'inference'}) {
	    $f_up->{$genome_id}->{$fid}->{'inference'} = +{};
	}
	if ($eggnog_line[$TARGET_I] && $eggnog_line[$TARGET_I] ne '-' && $eggnog_line[$TARGET_I] !~ /^\s*$/) {
	    $novel_fams_target = $eggnog_line[$TARGET_I];
	    $e_value = $eggnog_line[$EVALUE_I];
	    $bit_score = $eggnog_line[$BITSCORE_I];
	
	    $f_up->{$genome_id}->{$fid}->{'inference'}->{'novel_fams'} = +{};
	    $f_up->{$genome_id}->{$fid}->{'inference'}->{'novel_fams'}->{'category'} = 'eggnog5-novel_fams';
	    $f_up->{$genome_id}->{$fid}->{'inference'}->{'novel_fams'}->{'type'} = 'sequence homology';
	    $f_up->{$genome_id}->{$fid}->{'inference'}->{'novel_fams'}->{'evidence'} = join(",",
											    'target='.$novel_fams_target,
											    'e_value='.$e_value,
											    'bit_score='.$bit_score);
	}
    }
}


# create output buf
@feature_buf = ();
foreach $genome_id (sort keys %{$f_up}) {
    foreach $fid (sort keys %{$f_up->{$genome_id}}) {

	# aliases
	@aliases_buf = ();
	foreach $cat (qw(gene)) {
	    foreach $val (@{$f_up->{$genome_id}->{$fid}->{'aliases'}->{$cat}}) {
		push (@aliases_buf, '["'.$cat.'","'.$val.'"]');
	    }
	}
	$aliases = '"aliases":['.join(",",@aliases_buf).']';
	
	# functions
	@functions_buf = ();
	foreach $cat (qw(egg EC KO CAZy Pfam egg_desc)) {
	    foreach $val (@{$f_up->{$genome_id}->{$fid}->{'functions'}->{$cat}}) {
		push (@functions_buf, '"'.$cat.':'.$val.'"');
	    }
	}
	$functions = '"functions":['.join(",",@functions_buf).']';

	# inference
	@inference_buf = ();
	if (defined $f_up->{$genome_id}->{$fid}->{'inference'}->{'eggnog5'}) {
	    $inference_str = '{'.
		  '"category":"'.$f_up->{$genome_id}->{$fid}->{'inference'}->{'eggnog5'}->{'category'}.'",'.
		  '"type":"'.$f_up->{$genome_id}->{$fid}->{'inference'}->{'eggnog5'}->{'type'}.'",'.
		  
		  '"evidence":"'.$f_up->{$genome_id}->{$fid}->{'inference'}->{'eggnog5'}->{'evidence'}.'"'.
		  '}';
	    push (@inference_buf, $inference_str);
	}
	if (defined $f_up->{$genome_id}->{$fid}->{'inference'}->{'novel_fams'}) {
	    $inference_str = '{'.
		  '"category":"'.$f_up->{$genome_id}->{$fid}->{'inference'}->{'novel_fams'}->{'category'}.'",'.
		  '"type":"'.$f_up->{$genome_id}->{$fid}->{'inference'}->{'novel_fams'}->{'type'}.'",'.
		  
		  '"evidence":"'.$f_up->{$genome_id}->{$fid}->{'inference'}->{'novel_fams'}->{'evidence'}.'"'.
		  '}';
	    push (@inference_buf, $inference_str);
	}
	$inferences = '"inference_data":['.join(",",@inference_buf).']';
	
	push (@feature_buf, join("\t", $genome_id, $fid, $aliases, $functions, $inferences));
    }
}

# save feature buf
print ("SAVING feature val files...\n");
&writeBufToFile ($features_update_path, \@feature_buf);


# create params file
print ("SAVING params file...\n");
@params_buf = ();
push (@params_buf, '{');
push (@params_buf, '  "feature_update_file": "'.'/features/'.$features_update_file.'"');
push (@params_buf, '}');    
&writeBufToFile ($params_file, \@params_buf);


# create script file
print ("SAVING script file...\n");
@script_buf = ();
push (@script_buf, "#!/bin/bash\n\n");
push (@script_buf, qq{kb-sdk run -t beta --input $params_file --mount-points $features_dir:/features --sdk-home $image_dir kb_ObjectUtilities.KButil_update_genome_features_from_file});
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
\t-updatedir           <update_dir>              (e.g. 'direct_updates_features')
\t-setname             <set_name>                (e.g. 'GTDB_r207_Arc-chunk_0001-0050')
\t-eggnogannotfile     <eggnog_annot_file>       (e.g. 'GTDB_r207_Arc-chunk_0001-0050-diamond.emapper.annotations')
\t[-novelfamsannotfile  <novelfams_annot_file>]  (e.g. 'GTDB_r207_Arc-chunk_0001-0050-novel_fams.emapper.annotations')
};

    # Get args
    my %opts = ();
    &GetOptions (\%opts, 
		 "updatedir=s",
		 "setname=s",
		 "eggnogannotfile=s",
		 "novelfamsannotfile=s");

    # Check for legal invocation
    if (! defined $opts{updatedir} ||
	! defined $opts{setname} ||
	! defined $opts{eggnogannotfile}
	) {
        print STDERR "$usage\n";
        exit -1;
    }
    &checkExist ('d', $opts{updatedir});
    &checkExist ('f', $opts{eggnogannotfile});
    &checkExist ('f', $opts{novelfamsannotfile})  if (defined $opts{novelfamsannotfile});

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
