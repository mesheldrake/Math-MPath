#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 8;
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
my @ys_expect = (
    0.577897339560286,
    0.749036807344408,
    -0.546964473007376,
    -0.775174413677617,
    -0.551684003753151,
    0.74257556139109
);

is(scalar(@ys),6,'f(x) result count');
cmp_bag(
           [map {sprintf("%.10f",$_)} @ys],
           [map {sprintf("%.10f",$_)} @ys_expect],
           "f(x) output ok"
         );

my @xs = (
         $mp_arc_seg_a->F(-0.5), 
         $mp_arc_seg_a->F( 0.5),
         $mp_arc_seg_b->F(-0.5),
         $mp_arc_seg_b->F( 0.5),
         );
my @xs_expect = (
    1.56774310863924,
    -1.61557183829046,
    1.8254247910661,
    -1.87354355902701,
    1.57402621672775,
    1.81905910459074
);
is(scalar(@xs),6,'F(y) result count');
cmp_bag(
    [map {sprintf("%.10f",$_)} @xs],
    [map {sprintf("%.10f",$_)} @xs_expect],
    "F(y) output ok"
);

my @tangent_angles = map $mp_arc_seg_a->angleTangent(undef,undef,$_), (0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0);

is(scalar(grep {defined($_)} @tangent_angles),11,'enough tangent angles');

is_deeply([map {sprintf("%.10f",$_)} @tangent_angles],
    [map {sprintf("%.10f",$_)}
    (1.32388931581687,
    0.689157557918182,
    0.344472230674316,
    0.100751222746389,
    -0.138348898601983,
    -0.46364760900081,
    -1.05283401122026,
    -1.92408338506479,
    -2.5085749151812,
    -2.83197608730836,
    -3.07065713969832)
    ],
    "tangent angles CW"
);

@tangent_angles = map $mp_arc_seg_b->angleTangent(undef,undef,$_), (0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0);

is(scalar(grep {defined($_)} @tangent_angles),11,'enough tangent angles');

is_deeply([map {sprintf("%.10f",$_)} @tangent_angles],
    [map {sprintf("%.10f",$_)}
    (1.82801406372719,
    0.920131319441548,
    0.373211710344221,
    0.0651101133518759,
    -0.18079531154284,
    -0.46364760900081,
    -0.931091113117833,
    -1.7874542551126,
    -2.54297145662636,
    -2.93831708370588,
    3.08210465266726)
    ],
    "tangent angles CCW"
);

1;