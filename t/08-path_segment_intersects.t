#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 59;

use Math::MPath;

my $mp1;
my $mp2;
my $mp_line_seg;
my $mp_arc_seg;
my $mp_bez_seg;
my @intersections;

## LL ##
$mp1 = Math::MPath->newlite('M10,10L20,1',0.00001);
$mp2 = Math::MPath->newlite('M9,1  L21,11',0.00001);
@intersections = $mp1->getIntersections($mp2);
#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(    abs($intersections[0]->[0] - 14.71) < 0.01
    && abs($intersections[0]->[1] -  5.76) < 0.01
    , 'LL line segment intersection'
);

## LoLo ## offset lines - wrapper for LL
@intersections = Math::MPath::Intersections::intersect_LoLo(
    Math::MPath::LineSegment->new([10,10],[20, 1]),
    Math::MPath::LineSegment->new([ 9, 1],[21,11]),
    0.3,0.2,0
);
#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(    abs($intersections[0]->[0] - 14.6288) < 0.0001
    && abs($intersections[0]->[1] -  5.4303) < 0.0001
    , 'LoLo offset line segment intersection'
);

## AL ## general elliptical arc and line
$mp1 = Math::MPath->newlite('M10,10A30,31,0,0,0,20,1',0.00001);
$mp2 = Math::MPath->newlite('M9,1  L21,11',0.00001);
@intersections = $mp1->getIntersections($mp2);
#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(    abs($intersections[0]->[0] - 15.2972) < 0.0001
    && abs($intersections[0]->[1] -  6.2477) < 0.0001
    , 'AL circular arc and line one intersection'
);

## A1oLo ## offset circular arc and offset line (wrapper for AL)
$mp_arc_seg = Math::MPath::EllipticalArc->new([10,10],[30,30],0,0,0,[20,1],0.00001,1);
$mp_line_seg = Math::MPath::LineSegment->new([9,1],[21,11],0.00001,1);
@intersections = Math::MPath::Intersections::intersect_A1oLo($mp_arc_seg,$mp_line_seg,0.3,0.2,0,0);
#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(    abs($intersections[0]->[0] - 15.4485) < 0.0001
    && abs($intersections[0]->[1] -  6.1134) < 0.0001
    , 'A1oLo offset circular arc and offset line one intersection'
);

## A1A1 ## special case of two circular elliptical arcs
$mp1 = Math::MPath->newlite('M10,10A30,30,0,0,0,20,1',0.00001);
$mp2 = Math::MPath->newlite('M9,1  A30,30,0,0,0,21,11',0.00001);
@intersections = $mp1->getIntersections($mp2);
#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(    abs($intersections[0]->[0] - 14.5053) < 0.0001
    && abs($intersections[0]->[1] -  6.9334) < 0.0001
    , 'A1A1 circular arcs one intersection'
);

## A1oA1o ## special case of two offset circular elliptical arcs - wrapper for A1A1
@intersections = Math::MPath::Intersections::intersect_A1oA1o(
    Math::MPath::EllipticalArc->new([10,10],[30,30],0,0,0,[20, 1],0.00001,1),
    Math::MPath::EllipticalArc->new([ 9, 1],[30,30],0,0,0,[21,11],0.00001,1),
    0.3,0.2,0
);
#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(    abs($intersections[0]->[0] - 14.5053) < 0.0001
    && abs($intersections[0]->[1] -  6.9334) < 0.0001
    , 'A1oA1o offset circular arcs one intersection'
);

## AoAo ## two offset elliptical arcs, one intersection
@intersections = Math::MPath::Intersections::intersect_AoAo(
    Math::MPath::EllipticalArc->new([10,10],[30,31],0,0,0,[20, 1],0.00001,1),
    Math::MPath::EllipticalArc->new([ 9, 1],[30,31],0,0,0,[21,11],0.00001,1),
    0.3,0.2
);
#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(    abs($intersections[0]->[0] - 14.5850) < 0.0001
    && abs($intersections[0]->[1] -  7.2472) < 0.0001
    , 'AoAo one intersection'
);

## AoAo ## two offset elliptical arcs, two intersections within one sub section
@intersections = Math::MPath::Intersections::intersect_AoAo(
    Math::MPath::EllipticalArc->new([4.0,0],[30,31],0,0,0,[35.0,26.0],0.00001,1),
    Math::MPath::EllipticalArc->new([0,6.0],[30,31],0,0,1,[30.0,29.0],0.00001,1),
    0.3,0.2
);
#diag("intersection pair count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(scalar(@intersections) == 2, 'AoAo two intersections');
ok(    abs($intersections[0]->[0] -  5.4275) < 0.0001
    && abs($intersections[0]->[1] -  6.5224) < 0.0001
    , 'AoAo intersect pair, first'
);
ok(    abs($intersections[1]->[0] - 28.8309) < 0.0001
    && abs($intersections[1]->[1] - 25.9410) < 0.0001
    , 'AoAo intersect pair, second'
);

## AoAo ## as above, but first reversed
@intersections = Math::MPath::Intersections::intersect_AoAo(
    Math::MPath::EllipticalArc->new([35.0,26.0],[30,31],0,0,1,[4.0,0],0.00001,1),
    Math::MPath::EllipticalArc->new([0,6.0],[30,31],0,0,1,[30.0,29.0],0.00001,1),
    0.3,0.2
);
#diag("intersection pair count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(scalar(@intersections) == 2, 'AoAo, first reversed, two intersections');
ok(    abs($intersections[0]->[0] -  6.1159) < 0.0001
    && abs($intersections[0]->[1] -  6.6375) < 0.0001
    , 'AoAo, first reversed, intersect pair, first'
);
ok(    abs($intersections[1]->[0] - 28.5805) < 0.0001
    && abs($intersections[1]->[1] - 25.2893) < 0.0001
    , 'AoAo, first reversed, intersect pair, second'
);

## AoAo ## as above, but second reversed
@intersections = Math::MPath::Intersections::intersect_AoAo(
    Math::MPath::EllipticalArc->new([4.0,0],[30,31],0,0,0,[35.0,26.0],0.00001,1),
    Math::MPath::EllipticalArc->new([30.0,29.0],[30,31],0,0,0,[0,6.0],0.00001,1),
    0.3,0.2
);
#diag("intersection pair count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(scalar(@intersections) == 2, 'AoAo, second reversed, two intersections');
ok(    abs($intersections[0]->[0] -  5.2672) < 0.0001
    && abs($intersections[0]->[1] -  6.0936) < 0.0001
    , 'AoAo, second reversed, intersect pair, first'
);
ok(    abs($intersections[1]->[0] - 29.2832) < 0.0001
    && abs($intersections[1]->[1] - 26.0120) < 0.0001
    , 'AoAo, second reversed, intersect pair, second'
);
## CL ##
$mp1 = Math::MPath->newlite('M10,10C13,10,16,8,20,1',0.00001);
$mp2 = Math::MPath->newlite('M9,1  L21,11',0.00001);
@intersections = $mp1->getIntersections($mp2);
#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(    abs($intersections[0]->[0] - 15.9727) < 0.0001
    && abs($intersections[0]->[1] -  6.8106) < 0.0001
    , 'CL intersect'
);

## CoLo ##
$mp_bez_seg  = Math::MPath::BezierCubicSegment->new([4.0,0],[11.0,14.0],[23.0,24.0],[35.0,26.0],0.00001,1); # with islite flag
$mp_line_seg = Math::MPath::LineSegment->new(       [0,3.0],                        [30.0,26.0],0.00001,1); # with islite flag
@intersections = Math::MPath::Intersections::intersect_CoLo($mp_bez_seg,$mp_line_seg, 0.3, 0.2);
#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(    abs($intersections[0]->[0] - 11.5659) < 0.0001
    && abs($intersections[0]->[1] - 11.9326) < 0.0001
    , 'CoLo intersect pair, first'
);
ok(    abs($intersections[1]->[0] - 27.2036) < 0.0001
    && abs($intersections[1]->[1] - 23.9215) < 0.0001
    , 'CoLo intersect pair, second'
);

## CoLo ##
$mp_bez_seg  = Math::MPath::BezierCubicSegment->new([35.0,26.0],[23.0,24.0],[11.0,14.0],[4.0,0],0.00001,1); # with islite flag
$mp_line_seg = Math::MPath::LineSegment->new(       [3.0,4.0],                        [33.0,27.0],0.00001,1); # with islite flag
@intersections = Math::MPath::Intersections::intersect_CoLo($mp_bez_seg,$mp_line_seg, 0.3, 0.2);
#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(    abs($intersections[0]->[0] - 10.7367) < 0.0001
    && abs($intersections[0]->[1] -  9.9969) < 0.0001
    , 'CoLo, Co reversed, intersect pair, first'
);
ok(    abs($intersections[1]->[0] - 29.0594) < 0.0001
    && abs($intersections[1]->[1] - 24.0443) < 0.0001
    , 'CoLo, Co reversed, intersect pair, second'
);
## CoLo ##
$mp_bez_seg  = Math::MPath::BezierCubicSegment->new([4.0,0],[11.0,14.0],[23.0,24.0],[35.0,26.0],0.00001,1); # with islite flag
$mp_line_seg = Math::MPath::LineSegment->new(       [30.0,26.0],                    [ 0,   3.0],0.00001,1); # with islite flag
@intersections = Math::MPath::Intersections::intersect_CoLo($mp_bez_seg,$mp_line_seg, 0.3, 0.2);
#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(    abs($intersections[0]->[0] - 11.2407) < 0.0001
    && abs($intersections[0]->[1] - 11.5524) < 0.0001
    , 'CoLo, Lo reversed, intersect pair, first'
);
ok(    abs($intersections[1]->[0] - 27.6025) < 0.0001
    && abs($intersections[1]->[1] - 24.0965) < 0.0001
    , 'CoLo, Lo reversed, intersect pair, second'
);

## CoL ##
$mp_bez_seg  = Math::MPath::BezierCubicSegment->new([4.0,0],[11.0,14.0],[23.0,24.0],[35.0,26.0],0.00001,1); # with islite flag
$mp_line_seg = Math::MPath::LineSegment->new([0,3.0],[30.0,26.0],0.00001,1); # with islite flag
@intersections = Math::MPath::Intersections::intersect_CoL($mp_bez_seg,$mp_line_seg, 0.3);
#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(    abs($intersections[0]->[0] - 11.4011) < 0.0001
    && abs($intersections[0]->[1] - 11.7408) < 0.0001
    , 'CoL intersect pair, first'
);
ok(    abs($intersections[1]->[0] - 27.4052) < 0.0001
    && abs($intersections[1]->[1] - 24.0107) < 0.0001
    , 'CoL intersect pair, second'
);

## CA ## cubic Bezier and EllipticalArc intersection where there are two intersections between decomposed sub sections of each
@intersections = Math::MPath::Intersections::intersect_CA(
    Math::MPath::BezierCubicSegment->new([4.0,0],[11.0,14.0],[23.0,24.0],[35.0,26.0],0.00001,1),
    Math::MPath::EllipticalArc->new([0,6.0],[30,31],0,0,1,[30.0,29.0],0.00001,1)
);
#diag("intersection pair count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(scalar(@intersections) == 2, 'CA two intersections');
ok(    abs($intersections[0]->[0] -  8.0254) < 0.0001
    && abs($intersections[0]->[1] -  6.8414) < 0.0001
    , 'CA intersect pair, first'
);
ok(    abs($intersections[1]->[0] - 28.2791) < 0.0001
    && abs($intersections[1]->[1] - 24.0560) < 0.0001
    , 'CA intersect pair, second'
);

## CoAo ## cubic Bezier and elliptical arc intersection where there are two intersections between decomposed sub sections of each OFFSET curve
@intersections = Math::MPath::Intersections::intersect_CoAo(
    Math::MPath::BezierCubicSegment->new([4.0,0],[11.0,14.0],[23.0,24.0],[35.0,26.0],0.00001,1),
    Math::MPath::EllipticalArc->new([0,6.0],[30,31],0,0,1,[30.0,29.0],0.00001,1),
    0.3, 0.2
);
#diag("intersection pair count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(scalar(@intersections) == 2, 'CoAo two intersections');
ok( @intersections > 0
    && abs($intersections[0]->[0] -  7.7580) < 0.0001
    && abs($intersections[0]->[1] -  6.9820) < 0.0001
    , 'CoAo intersect pair, first'
);

ok( @intersections > 1
    && abs($intersections[1]->[0] - 28.1865) < 0.0001
    && abs($intersections[1]->[1] - 24.3423) < 0.0001
    , 'CoAo intersect pair, second'
);

## CoAo ## as above, but just Bezier control points reversed
@intersections = Math::MPath::Intersections::intersect_CoAo(
    Math::MPath::BezierCubicSegment->new([35.0,26.0],[23.0,24.0],[11.0,14.0],[4.0,0],0.00001,1),
    Math::MPath::EllipticalArc->new([0,6.0],[30,31],0,0,1,[30.0,29.0],0.00001,1),
    0.3, 0.2
);
#diag("intersection pair count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(scalar(@intersections) == 2, 'CoAo, Co reversed, two intersections');
ok( @intersections > 0
    && abs($intersections[0]->[0] -  8.5357) < 0.0001
    && abs($intersections[0]->[1] -  7.0521) < 0.0001
    , 'CoAo, Co reversed, intersect pair, first'
);

ok( @intersections > 1
    && abs($intersections[1]->[0] - 28.1136) < 0.0001
    && abs($intersections[1]->[1] - 23.6622) < 0.0001
    , 'CoAo, Co reversed, intersect pair, second'
);

## CoAo ## as above, but just arc reversed - reverse end points AND sweep flag
@intersections = Math::MPath::Intersections::intersect_CoAo(
    Math::MPath::BezierCubicSegment->new([4.0,0],[11.0,14.0],[23.0,24.0],[35.0,26.0],0.00001,1),
    Math::MPath::EllipticalArc->new([30.0,29.0],[30,31],0,0,0,[0,6.0],0.00001,1),
    0.3, 0.2
);
#diag("intersection pair count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(scalar(@intersections) == 2, 'CoAo, Ao reversed, two intersections');
ok( @intersections > 0
    && abs($intersections[0]->[0] -  7.4282) < 0.0001
    && abs($intersections[0]->[1] -  6.4947) < 0.0001
    , 'CoAo, Ao reversed, intersect pair, first'
);

ok( @intersections > 1
    && abs($intersections[1]->[0] - 28.7116) < 0.0001
    && abs($intersections[1]->[1] - 24.5528) < 0.0001
    , 'CoAo, Ao reversed, intersect pair, second'
);

## CC ## cubic Bezier loop self-intersect case
$mp_bez_seg = Math::MPath::BezierCubicSegment->new([0,0],[200,80],[-100,80],[100,0],0.00001,1); # with islite flag
@intersections = Math::MPath::Intersections::intersect_CC($mp_bez_seg,$mp_bez_seg);
#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(    abs($intersections[0]->[0] - 50) < 0.0001
    && abs($intersections[0]->[1] - 24) < 0.0001
    , 'CC self intersect'
);

## CC ## cubic Bezier intersection where there are two intersections between decomposed sub sections of beziers
@intersections = Math::MPath::Intersections::intersect_CC(
    Math::MPath::BezierCubicSegment->new([4.0,0],[11.0,14.0],[23.0,24.0],[35.0,26.0],0.00001,1),
    Math::MPath::BezierCubicSegment->new([0,6.0],[13.0, 9.0],[24.0,19.0],[30.0,29.0],0.00001,1)
);
#diag("intersection pair count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(scalar(@intersections) == 2, 'CC pair count == 2');
ok(    abs($intersections[0]->[0] - 10.3517) < 0.0001
    && abs($intersections[0]->[1] -  9.9952) < 0.0001
    , 'CC pair first'
);
ok(    abs($intersections[1]->[0] - 25.6583) < 0.0001
    && abs($intersections[1]->[1] - 22.8527) < 0.0001
    , 'CC pair second'
);

## CC ## as above, but with control points in reverse order for second Bezier
#        so that increasing t gives decreasing x for that one.
@intersections = Math::MPath::Intersections::intersect_CC(
    Math::MPath::BezierCubicSegment->new([4.0,0],[11.0,14.0],[23.0,24.0],[35.0,26.0],0.00001,1),
    Math::MPath::BezierCubicSegment->new([30.0,29.0],[24.0,19.0],[13.0, 9.0],[0,6.0],0.00001,1)
);
#diag("intersection pair count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(scalar(@intersections) == 2, 'CC pair, 2nd reversed, count == 2');
ok(    abs($intersections[0]->[0] - 10.3517) < 0.0001
    && abs($intersections[0]->[1] -  9.9952) < 0.0001
    , 'CC pair, 2nd reversed, first'
);
ok(    abs($intersections[1]->[0] - 25.6583) < 0.0001
    && abs($intersections[1]->[1] - 22.8527) < 0.0001
    , 'CC pair, 2nd reversed, second'
);

## CoCo ## same looped bezier but two different offsets; should intersect two times
$mp_bez_seg = Math::MPath::BezierCubicSegment->new([0,0],[200,80],[-100,80],[100,0],0.00001,1); # with islite flag
@intersections = Math::MPath::Intersections::intersect_CoCo($mp_bez_seg,$mp_bez_seg, 0.3, 0.2);
#diag("intersection count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok( scalar(@intersections) == 2 , 'CoCo self intersect, two intersections');
ok( @intersections
    && abs($intersections[0]->[0] - 49.9052) < 0.0001
    && abs($intersections[0]->[1] - 24.2942) < 0.0001
    , 'CoCo self intersect, 1st intersection coords check'
);
ok( @intersections
    && abs($intersections[1]->[0] - 50.0947) < 0.0001
    && abs($intersections[1]->[1] - 24.2942) < 0.0001
    , 'CoCo self intersect, 2nd intersection coords check'
);

## CoCo ## cubic Bezier intersection where there are two intersections between decomposed sub sections of OFFSET beziers
@intersections = Math::MPath::Intersections::intersect_CoCo(
    Math::MPath::BezierCubicSegment->new([4.0,0],[11.0,14.0],[23.0,24.0],[35.0,26.0],0.00001,1),
    Math::MPath::BezierCubicSegment->new([0,6.0],[13.0,9.0] ,[24.0,19.0],[30.0,29.0],0.00001,1),
    0.3, 0.2
);
#diag("intersection pair count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(scalar(@intersections) == 2, 'CoCo pair count == 2');
ok( @intersections > 0
    && abs($intersections[0]->[0] -  9.9973) < 0.0001
    && abs($intersections[0]->[1] - 10.0306) < 0.0001
    , 'CoCo pair, first'
);

ok( @intersections > 1
    && abs($intersections[1]->[0] - 25.6850) < 0.0001
    && abs($intersections[1]->[1] - 23.2022) < 0.0001
    , 'CoCo pair, second'
);

## CoCo ## same as above, but with control points reversed on 2nd Bezier, so
#          increasing t gives decreasing x
@intersections = Math::MPath::Intersections::intersect_CoCo(
    Math::MPath::BezierCubicSegment->new([4.0,0],[11.0,14.0],[23.0,24.0],[35.0,26.0],0.00001,1),
    Math::MPath::BezierCubicSegment->new([30.0,29.0],[24.0,19.0],[13.0,9.0],[0,6.0],0.00001,1),
    0.3, 0.2
);
#diag("intersection pair count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(scalar(@intersections) == 2, 'CoCo pair, 2nd reversed, count == 2');
ok( @intersections > 0
    && abs($intersections[0]->[0] -  9.4082) < 0.0001
    && abs($intersections[0]->[1] -  9.2694) < 0.0001
    , 'CoCo pair, 2nd reversed, first'
);

ok( @intersections > 1
    && abs($intersections[1]->[0] - 26.5203) < 0.0001
    && abs($intersections[1]->[1] - 23.6083) < 0.0001
    , 'CoCo pair, 2nd reversed, second'
);

## CoCo ## this time control points reversed on 1st Bezier
@intersections = Math::MPath::Intersections::intersect_CoCo(
    Math::MPath::BezierCubicSegment->new([35.0,26.0],[23.0,24.0],[11.0,14.0],[4.0,0],0.00001,1),
    Math::MPath::BezierCubicSegment->new([0,6.0],[13.0,9.0],[24.0,19.0],[30.0,29.0],0.00001,1),
    0.3, 0.2
);
#diag("intersection pair count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(scalar(@intersections) == 2, 'CoCo pair, 1st reversed, count == 2');
ok( @intersections > 0
    && abs($intersections[0]->[0] - 11.4126) < 0.0001
    && abs($intersections[0]->[1] - 10.8243) < 0.0001
    , 'CoCo pair, 1st reversed, first'
);

ok( @intersections > 1
    && abs($intersections[1]->[0] - 24.6784) < 0.0001
    && abs($intersections[1]->[1] - 21.9952) < 0.0001
    , 'CoCo pair, 1st reversed, second'
);

## CoCo ## this time control points reversed on both Beziers
@intersections = Math::MPath::Intersections::intersect_CoCo(
    Math::MPath::BezierCubicSegment->new([35.0,26.0],[23.0,24.0],[11.0,14.0],[4.0,0],0.00001,1),
    Math::MPath::BezierCubicSegment->new([30.0,29.0],[24.0,19.0],[13.0,9.0],[0,6.0],0.00001,1),
    0.3, 0.2
);
#diag("intersection pair count:".scalar(@intersections)."\n".join("\n",map {'['.join(',',@$_).']'} @intersections));
ok(scalar(@intersections) == 2, 'CoCo pair, both reversed, count == 2');
ok( @intersections > 0
    && abs($intersections[0]->[0] - 10.7108) < 0.0001
    && abs($intersections[0]->[1] -  9.9646) < 0.0001
    , 'CoCo pair, both reversed, first'
);

ok( @intersections > 1
    && abs($intersections[1]->[0] - 25.6269) < 0.0001
    && abs($intersections[1]->[1] - 22.4997) < 0.0001
    , 'CoCo pair, both reversed, second'
);

1;
