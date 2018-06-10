#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 3;

use Math::MPath;


### SINGLE SELF INTERSECT CASE

my $mp_bez_seg = Math::MPath::BezierCubicSegment->new([0,0],[200,80],[-100,80],[100,0],0.00001,1); # with islite flag

my @intersections = Math::MPath::Intersections::bez_bez_intersect($mp_bez_seg,$mp_bez_seg);

ok(    $intersections[0]->[0] - 50 < 0.0001
    && $intersections[0]->[1] - 24 < 0.0001
);

#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));


### INTERSECT CASE WHERE YOU HAVE TWO INTERSECTIONS WITHIN ONE XSPAN

my $mp_bez_seg_1 = Math::MPath::BezierCubicSegment->new([4.0,0],[11.0,14.0],[23.0,24.0],[35.0,26.0],0.00001,1); # with islite flag
my $mp_bez_seg_2 = Math::MPath::BezierCubicSegment->new([0,6.0],[13.0,9.0] ,[24.0,19.0],[30.0,29.0],0.00001,1); # with islite flag
# intersections somewhere around (10.35,10.0) and (25.72,22.82) (eyballing it in inkscape) (could be wrong, inkscape coord display is sometimes confusing)

my @intersection_pair = Math::MPath::Intersections::bez_bez_intersect($mp_bez_seg_1,$mp_bez_seg_2);

ok(    $intersection_pair[0]->[0] - 10.3517 < 0.0001
    && $intersection_pair[0]->[1] -  9.9952 < 0.0001
);
ok(    $intersection_pair[1]->[0] - 25.6583 < 0.0001
    && $intersection_pair[1]->[1] - 22.8527 < 0.0001
);

#diag("intersection pair count:".scalar(@intersection_pair)."\n".join("\n",map {'['.join(',',@$_).']'} @intersection_pair));

1;
