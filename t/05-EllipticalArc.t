#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 1;

use Math::MPath;

# delta theta is negative for this one - CW in y-up frame
my $mp_arc_seg_a = Math::MPath::EllipticalArc->new([-2,0],[2,1],5,1,0,[0,-1],0.00001,1);
# delta theta is positive for this one - CCW in y-up frame
my $mp_arc_seg_b = Math::MPath::EllipticalArc->new([-2,0],[2,1],5,1,1,[0, 1],0.00001,1);

my @ys = (
         $mp_arc_seg_a->f(-1.5), 
         $mp_arc_seg_a->f( 1.5),
         $mp_arc_seg_b->f(-1.5),
         $mp_arc_seg_b->f( 1.5),
         );


ok(1);

1;