#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 2;

use Math::MPath;

my $mp = Math::MPath->new('M0,0L3,4',0.00001);
ok($mp->{pathspec});

my $mp_bez_cubic = Math::MPath::BezierCubicSegment->new([0.0,0.0],[200.0,80.0],[-100.0,80.0],[100.0,0.0],0.00001,1); # with islite flag
ok($mp_bez_cubic->{p1});

