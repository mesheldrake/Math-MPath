#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 4;
use Test::Deep;

use Math::MPath;

# NOTE that though SVG uses a +Y goes down coordinate system, we generally
# concieve of our path math as operating in a +Y goes up coordinate system,
# like the one typically used in basic math and geometry illustrations.
# This means clockwise and counterclockwise, and left and right, as used here,
# would appear reversed if the geometry here were displayed in an SVG context.

# Wide ellipse, tilted 5 degrees, going about 3/4 the way around, CLOCKWISE (in a +Y goes up coordinate system)
# delta theta is negative for this one - CW in y-up frame
my $mp_arc_seg_a = Math::MPath::EllipticalArc->new([-2,0],[2,1],5,1,0,[0,-1],0.00001,1);
# Wide ellipse, tilted 5 degrees, going about 3/4 the way around, COUNTERCLOCKWISE (in a +Y goes up coordinate system)
# delta theta is positive for this one - CCW in y-up frame
my $mp_arc_seg_b = Math::MPath::EllipticalArc->new([-2,0],[2,1],5,1,1,[0, 1],0.00001,1);

my @ys = (
         $mp_arc_seg_a->f(-1.5), 
         $mp_arc_seg_a->f( 1.5),
         $mp_arc_seg_b->f(-1.5),
         $mp_arc_seg_b->f( 1.5),
         );
my @ys_old = (
         $mp_arc_seg_a->f_old(-1.5), 
         $mp_arc_seg_a->f_old( 1.5),
         $mp_arc_seg_b->f_old(-1.5),
         $mp_arc_seg_b->f_old( 1.5),
         );

is(scalar(@ys),6,'f(x) result count');
cmp_bag(
           [map {sprintf("%.10f",$_)} @ys],
           [map {sprintf("%.10f",$_)} @ys_old],
           "new f(x) output same as old f(x)"
         );

my @xs = (
         $mp_arc_seg_a->F(-0.5), 
         $mp_arc_seg_a->F( 0.5),
         $mp_arc_seg_b->F(-0.5),
         $mp_arc_seg_b->F( 0.5),
         );
my @xs_old = (
         $mp_arc_seg_a->F_old(-0.5), 
         $mp_arc_seg_a->F_old( 0.5),
         $mp_arc_seg_b->F_old(-0.5),
         $mp_arc_seg_b->F_old( 0.5),
         );

is(scalar(@xs),6,'F(y) result count');
cmp_bag(
           [map {sprintf("%.10f",$_)} @xs],
           [map {sprintf("%.10f",$_)} @xs_old],
           "new F(y) output same as old F(y)"
         );

1;