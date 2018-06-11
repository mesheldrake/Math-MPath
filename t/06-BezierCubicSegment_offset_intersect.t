#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 6;

use Math::MPath;




### SAME LOOPED BEZIER, TWO DIFFERNT OFFSETS, SHOULD INTERSECT TWO TIMES

my $mp_bez_seg = Math::MPath::BezierCubicSegment->new([0,0],[200,80],[-100,80],[100,0],0.00001,1); # with islite flag

my @intersections = Math::MPath::Intersections::intersect_CoCo($mp_bez_seg,$mp_bez_seg, 0.3, 0.2);

ok( scalar(@intersections) == 2 , 'two intersections');

ok( @intersections
    && abs($intersections[0]->[0] - 49.9052) < 0.0001
    && abs($intersections[0]->[1] - 24.2942) < 0.0001
    , '1st intersection coords check'
);
ok( @intersections
    && abs($intersections[1]->[0] - 50.0947) < 0.0001
    && abs($intersections[1]->[1] - 24.2942) < 0.0001
    , '2nd intersection coords check'
);

#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));



### INTERSECT CASE WHERE YOU HAVE TWO INTERSECTIONS WITHIN ONE XSPAN

my $mp_bez_seg_1 = Math::MPath::BezierCubicSegment->new([4.0,0],[11.0,14.0],[23.0,24.0],[35.0,26.0],0.00001,1); # with islite flag
my $mp_bez_seg_2 = Math::MPath::BezierCubicSegment->new([0,6.0],[13.0,9.0] ,[24.0,19.0],[30.0,29.0],0.00001,1); # with islite flag
# intersections somewhere around (10.35,10.0) and (25.72,22.82) (eyballing it in inkscape) (could be wrong, inkscape coord display is sometimes confusing)

my @intersection_pair = Math::MPath::Intersections::intersect_CoCo($mp_bez_seg_1,$mp_bez_seg_2, 0.3, 0.2);

ok(scalar(@intersection_pair) == 2);

ok( @intersection_pair > 0
    && abs($intersection_pair[0]->[0] -  9.9973) < 0.0001
    && abs($intersection_pair[0]->[1] - 10.0306) < 0.0001
);

ok( @intersection_pair > 1
    && abs($intersection_pair[1]->[0] - 25.6850) < 0.0001
    && abs($intersection_pair[1]->[1] - 23.2022) < 0.0001
);

#diag("intersection pair count:".scalar(@intersection_pair)."\n".join("\n",map {'['.join(',',@$_).']'} @intersection_pair));



1;
