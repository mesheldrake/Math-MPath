#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 3;

use Math::MPath;

my $mp_bez_seg = Math::MPath::BezierCubicSegment->new([0,0],[200,80],[-100,80],[100,0],0.00001,1); # with islite flag

my @ys = $mp_bez_seg->f(48.0);
my @ys_old = $mp_bez_seg->f_old(48.0);

#diag("ys:\n".join("\n",map {$ys[$_] . ' =?= ' . $ys_old[$_]} (0..$#ys))."\n");
is(scalar(@ys), 3);
is_deeply([sort {$a<=>$b} @ys], [sort {$a<=>$b} @ys_old]);

my @ts = $mp_bez_seg->t(48.0);
my @ts_old = $mp_bez_seg->solveXforTheta(48.0);

#diag("ts(48.0):\n".join("\n",map {$ts[$_] . ' =?= ' . $ts_old[$_]} (0..$#ts))."\n");
is_deeply([sort {$a<=>$b} @ts], [sort {$a<=>$b} @ts_old]);

1;