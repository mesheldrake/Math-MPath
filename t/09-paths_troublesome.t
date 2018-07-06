#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 1;

use Math::MPath;

# use similar process to MPath::constructSegments()
# (copy-paste the latest of that in here and adapt)
# but go step by step so you can debug troublesome paths



my $trouble_path = "M -654.28571,-127.63782 C -653.33333,-123.82829 -652.38095,-120.01877 -651.42857,-116.20925 C -651.42857,-114.30448 -651.80213,-112.36273 -651.42857,-110.49496 C -650.83793,-107.54175 -648.99735,-104.90496 -648.57143,-101.92353 C -648.03268,-98.152296 -649.49537,-94.190741 -648.57143,-90.49496 C -647.53843,-86.362952 -644.08101,-83.145941 -642.85714,-79.066389 C -641.19251,-73.517598 -641.13612,-67.604141 -640,-61.923532 C -639.2299,-58.07302 -637.99469,-54.32822 -637.14286,-50.49496 C -636.15545,-46.051626 -633.32893,-29.848285 -631.42857,-24.780675 C -629.93307,-20.792683 -627.39206,-17.266905 -625.71429,-13.352103 C -621.48566,-3.4853033 -619.7326,7.7174886 -617.14286,18.076468 C -616.19048,21.885992 -615.2381,25.695516 -614.28571,29.50504 C -613.33333,33.314564 -611.78408,37.02297 -611.42857,40.933611 C -610.825,47.572899 -612.37138,54.333948 -611.42857,60.933611 C -603.40281,117.11394 -607.14298,45.217723 -602.85714,100.93361 C -602.41888,106.63107 -603.87935,112.45435 -602.85714,118.07647 C -601.93969,123.12247 -598.38676,127.38659 -597.14286,132.36218 C -596.72383,134.03829 -597.84292,148.10491 -597.14286,149.50504 C -593.61228,156.5662 -594.28571,148.51882 -594.28571,155.21933";

my $self={};
bless $self,'Math::MPath';
$self->{resolution} = 0.00001;
$self->{precision} = $self->{resolution}/1000;
$self->{isLite} = 1;



#diag('trouble bez');
#my $bez = Math::MPath::BezierCubicSegment->new([-617.14286,18.076468],[-616.19048,21.885992],[-615.2381,25.695516],[-614.28571,29.50504],$self->{precision},$self->{isLite});
#diag('trouble bez do X_offset');
#my $lutind = 0;
#diag " [$bez->{p1}->[0],$bez->{p1}->[1]],[$bez->{cp1}->[0],$bez->{cp1}->[1]],[$bez->{cp2}->[0],$bez->{cp2}->[1]],[$bez->{p2}->[0],$bez->{p2}->[1]], $bez->{XtoTLUT}->[$lutind]->[2]->[-1],\n";
#my $xoff = $bez->X_offset(0,undef,$bez->{XtoTLUT}->[$lutind]->[0]->[1],$bez->{XtoTLUT}->[$lutind]->[3]);
#diag('trouble bez DONE ',$xoff);
#exit;

$self->{pathspec} = $trouble_path;
$self->{pathSegmentSpecs} = [];
@{$self->{pathSegmentSpecs}} = $self->{pathspec} =~ /[MmZzLlHhVvCcSsQqTtAa][0-9, \-\.e]*/g;
$self->{pathSegments}=[];
my $lastM=$self->{pathSegmentSpecs}->[0];
my $lastSegSpec=$self->{pathSegmentSpecs}->[0];
for (my $i=1;$i<scalar(@{$self->{pathSegmentSpecs}});$i++) {
    my $thisSegSpec=$self->{pathSegmentSpecs}->[$i];
    if ($self->{pathSegmentSpecs}->[$i] =~ /^M/i) {$lastM=$self->{pathSegmentSpecs}->[$i];} #so we can start to be smart about paths with subpaths (multiple Ms)
    if ($self->{pathSegmentSpecs}->[$i] =~ /^(Z)/i) {$thisSegSpec=$1.substr($lastM,1);} # so we can treat it as a LineSegment, and start to be smart about paths with subpaths (multiple Ms). Ee should flag it as special case somehow, but we don't yet
    if ($lastSegSpec =~ /^Z/i && $thisSegSpec !~ /^M/i) {$lastSegSpec=$lastM;} #per SVG spec - if no new M following Z, but something else, that something else is from last M

    if ($thisSegSpec=~/^H/i) {
        $thisSegSpec=~s/\s//g; #/
        if ($i==1) {
            $lastM=~/M\s*?[0-9\-\.eE]+\s*?[\s,]\s*?([0-9\-\.eE]+)/;
            $thisSegSpec.=','.$1;
            }
        else {$thisSegSpec.=','.$self->{pathSegments}->[$i - 2]->{p2}->[1];}
        $thisSegSpec=~s/^H/L/;
        }
    if ($thisSegSpec=~/^V/i) {
        if ($i==1) {
            $lastM=~/M\s*?([0-9\-\.eE]+)\s*?[\s,]\s*?[0-9\-\.eE]+/;
            my $h = $1;
            $thisSegSpec=~s/^V(.*)$/'L'.$h.','.$1/ei;
            }
        else {
            $thisSegSpec=~s/^V(.*)$/'L'.$self->{pathSegments}->[$i - 2]->{p2}->[0].','.$1/ei;
            }
        }

    #diag("seg[",($i-1),"]: ",$thisSegSpec);
    $self->{pathSegments}->[$i - 1] = $self->constructSegment($thisSegSpec,$lastSegSpec);
    #diag("seg[",($i-1),"] DONE\n\n");

    $lastSegSpec=$thisSegSpec;
    }



ok(1);
