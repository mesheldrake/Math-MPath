#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 8;

use Math::MPath;

my $mp1;
my $mp2;
my @intersections;

## LL ##
$mp1 = Math::MPath->newlite('M10,10L20,1',0.00001);
$mp2 = Math::MPath->newlite('M9,1  L21,11',0.00001);
@intersections = $mp1->getIntersections($mp2);
#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(    abs($intersections[0]->[0] - 14.71) < 0.01
    && abs($intersections[0]->[1] -  5.76) < 0.01
);

## AL ## general elliptical arc and line
$mp1 = Math::MPath->newlite('M10,10A30,31,0,0,0,20,1',0.00001);
$mp2 = Math::MPath->newlite('M9,1  L21,11',0.00001);
@intersections = $mp1->getIntersections($mp2);
#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(    abs($intersections[0]->[0] - 15.2972) < 0.0001
    && abs($intersections[0]->[1] -  6.2477) < 0.0001
);

## A1A1 ## special case of two circular elliptical arcs
$mp1 = Math::MPath->newlite('M10,10A30,30,0,0,0,20,1',0.00001);
$mp2 = Math::MPath->newlite('M9,1  A30,30,0,0,0,21,11',0.00001);
@intersections = $mp1->getIntersections($mp2);
#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(    abs($intersections[0]->[0] - 14.5053) < 0.0001
    && abs($intersections[0]->[1] -  6.9334) < 0.0001
);

## CL ##
$mp1 = Math::MPath->newlite('M10,10C13,10,16,8,20,1',0.00001);
$mp2 = Math::MPath->newlite('M9,1  L21,11',0.00001);
@intersections = $mp1->getIntersections($mp2);
#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(    abs($intersections[0]->[0] - 15.9727) < 0.0001
    && abs($intersections[0]->[1] -  6.8106) < 0.0001
);

## CoLo ##
my $mp_bez_seg  = Math::MPath::BezierCubicSegment->new([4.0,0],[11.0,14.0],[23.0,24.0],[35.0,26.0],0.00001,1); # with islite flag
my $mp_line_seg = Math::MPath::LineSegment->new([0,3.0],[30.0,26.0],0.00001,1); # with islite flag
@intersections = Math::MPath::Intersections::intersect_CoLo($mp_bez_seg,$mp_line_seg, 0.3, 0.2);
#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(    abs($intersections[0]->[0] - 11.5659) < 0.0001
    && abs($intersections[0]->[1] - 11.9326) < 0.0001
);
ok(    abs($intersections[1]->[0] - 27.2036) < 0.0001
    && abs($intersections[1]->[1] - 23.9215) < 0.0001
);
## CoL ##
$mp_bez_seg  = Math::MPath::BezierCubicSegment->new([4.0,0],[11.0,14.0],[23.0,24.0],[35.0,26.0],0.00001,1); # with islite flag
$mp_line_seg = Math::MPath::LineSegment->new([0,3.0],[30.0,26.0],0.00001,1); # with islite flag
@intersections = Math::MPath::Intersections::intersect_CoL($mp_bez_seg,$mp_line_seg, 0.3);
#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(    abs($intersections[0]->[0] - 11.4011) < 0.0001
    && abs($intersections[0]->[1] - 11.7408) < 0.0001
);
ok(    abs($intersections[1]->[0] - 27.4052) < 0.0001
    && abs($intersections[1]->[1] - 24.0107) < 0.0001
);

1;
