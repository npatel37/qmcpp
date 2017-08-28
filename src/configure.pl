#!/usr/bin/perl
=pod
Copyright (c) 2015-2016, UT-Battelle, LLC
All rights reserved

[Determinantal QMC++, Version 0.]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

=cut
use warnings;
use strict;

use lib "../../PsimagLite/scripts";
use Make;
my ($arg) = @ARGV;

if (defined($arg) and -r "$arg" and $arg ne "Config.make") {
	my $cmd = "cp Config.make Config.make.bak";
	system($cmd);
	$cmd = "cp $arg Config.make";
	system($cmd);
}

my %provenanceDriver = (name => 'Provenance', aux => 1);
my $dotos = "qmc.o Provenance.o";
my %qptDriver = (name => 'qmc', dotos => $dotos);

my @drivers = (\%provenanceDriver,\%qptDriver);

createMakefile();

sub createMakefile
{
	Make::backupMakefile();
	if (!(-r "Config.make")) {
		my $cmd = "cp Config.make.sample Config.make";
		system($cmd);
		print STDERR "$0: Executed $cmd\n";
	}

	my $fh;
	open($fh,">Makefile") or die "Cannot open Makefile for writing: $!\n";

	my %additionals;
	$additionals{"code"} = "QMC++";
	Make::newMake($fh,\@drivers,\%additionals);
	print STDERR "File Makefile has been written\n";
}

