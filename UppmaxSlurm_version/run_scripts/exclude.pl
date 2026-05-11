#!/usr/bin/perl -w

my %skip = ();
open (IN, "exclusionlist");
while (<IN>) {
  chomp;
  ++$skip{$_};
}

open (OUTR, ">rejectedlines");
while (<>) {
  my @bit = split;
  if ($skip{$bit[9]} ||
      $bit[10] =~ /^Simp|^Sat|^Low|RNA|Alu/ ||
      $bit[9] =~ /^L1P|^MST|^THE|^L1MA[1-8]A?$/ ) {
    print OUTR;
  } else {
    if ($bit[9] =~ /^L1M[ABCD]/ && $bit[1] < 18 && $bit[6] - $bit[5] > 1000) {
      print OUTR;
    } else {
      print;
    }
  }
}
