#!/usr/bin/perl

open (OLD_METADATA_FILE, $ARGV[0]);
$old_header_line = <OLD_METADATA_FILE>;
close (OLD_METADATA_FILE);

open (NEW_METADATA_FILE, $ARGV[1]);
$new_header_line = <NEW_METADATA_FILE>;
close (NEW_METADATA_FILE);

@old_header = split ("\t", $old_header_line);
@new_header = split ("\t", $new_header_line);

for (my $i = 0; $i <= $#old_header; $i++) {
    if ($old_header[$i] ne $new_header[$i]) {
	print ("I: $i, old: $old_header[$i], new: $new_header[$i]\n");
    }
}
