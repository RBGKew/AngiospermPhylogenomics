#!/usr/bin/perl -w
# Given a tree topology, write out a constraint alignment

use strict;
use Bio::TreeIO;
use Bio::Tree::Tree;
use Bio::Tree::Node;
use Getopt::Long;

sub Splits;

{
    die "Run as a filter: TreeToConstraints.pl < newick_file > constraint_alignment\n"
	unless @ARGV==0;
    my $io = new Bio::TreeIO('-format' => 'newick', '-fh' => \*STDIN);
    my $tree = $io->next_tree;
    die "Cannot parse tree" if !defined $tree;

    my %leafid = map {$_->id => 1} grep {$_->is_Leaf} $tree->get_nodes();
    my %constraints = (); # leaf id => iConstraint => value
    my $nConstraints = 0;
    my $splits = Splits($tree, 'list');
    foreach my $split (@$splits) {
	my @member = @$split;
	my %member = map {$_ => 1} @member;
	my @nonmember = grep {!exists $member{$_}} (keys %leafid);
	foreach (@member) {
	    $constraints{$_}{$nConstraints} = 1;
	}
	foreach (@nonmember) {
	    $constraints{$_}{$nConstraints} = 0;
	}
	$nConstraints++;
    }

    die "No constraints selected\n" if $nConstraints == 0;

    my @names = sort (keys %constraints);
    my %seqs = ();
    while(my ($name,$leafconstraints) = each %constraints) {
	my $seq = "";
	for (my $i = 0; $i < $nConstraints; $i++) {
	    my $value = $leafconstraints->{$i};
	    $value = "-" if !defined $value;
	    $seq .= $value;
	}
	$seqs{$name} = $seq;
    }
    print STDERR "Identified $nConstraints constraints for " . scalar(@names) . " sequences\n";
    foreach my $name (@names) {
	print ">$name\n";
	my @whole = $seqs{$name} =~ /.{1,60}/g;
	print join("\n",@whole) . "\n";
    }
}

sub Splits {
    my ($tree,$mode) = @_;
    $mode = "hash" if !defined $mode;
    my $savenode = $mode eq "hashnode" ? 1 : 0;
    my %splits = ();
    my @splits = ();
    my @taxa = map {$_->id} grep {$_->is_Leaf} $tree->get_nodes;

    foreach my $node ($tree->get_nodes) {
	next if !defined $node->ancestor || $node->is_Leaf;
	my %split = map {$_->id => 1} grep {$_->is_Leaf} $node->get_all_Descendents;
	my @with0 = sort grep {!exists $split{$_}} @taxa;
	my @with1 = sort grep {exists $split{$_}} @taxa;
	
	if (@with0>1 && @with1>1) {
	    my $key0 = join(",",@with0);
	    my $key1 = join(",",@with1);
	    if (!exists $splits{$key0} && !exists $splits{$key1}) {
		push @splits, \@with0;
	    }
	    $splits{$key0} = $savenode ? $node : 1;
	    $splits{$key1} = $savenode ? $node : 1;
	}
    }
    return $mode eq "hash" || $mode eq "hashnode" ? \%splits : \@splits;
}
