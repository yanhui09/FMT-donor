#!/usr/bin/perl

use strict;
use warnings;

my $read_name="start";
my $ref_name="start";

#open ("IN", "<", "$ARGV[0]") or die "Cannot open file $ARGV[0]";
while (defined (my $line = <STDIN>)) {
	chomp $line;
	if ($line =~ m/^@/) { ##print header section of the bam file
		print "$line\n";
	}
	
	else {					##fields of bam file
		my @tmp=split(/\s+/,$line);
		if (($read_name eq "start") && ($ref_name eq "start")) {
			print "$line\n";
			$read_name=$tmp[0];
			$ref_name=$tmp[2];
		}
		if ($ref_name ne "*") {
		if (($read_name ne $tmp[0]) && ($ref_name ne $tmp[2])) {
				print "$line\n";
		}
		}
		else {
			print "$line\n";
		}
		
		$read_name=$tmp[0];
		$ref_name=$tmp[2];
	}
}