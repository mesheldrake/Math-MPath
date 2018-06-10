#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 1;

use Math::MPath;

## LL ##
my $mp1 = Math::MPath->newlite('M10,10L20,1',0.00001);
my $mp2 = Math::MPath->newlite('M9,1  L21,11',0.00001);
my @intersections = $mp1->getIntersections($mp2);
ok(    abs($intersections[0]->[0] - 14.71) < 0.01
    && abs($intersections[0]->[1] -  5.76) < 0.01
);



1;
