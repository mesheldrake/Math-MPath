#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 1;

use Math::MPath;

my $mp_d = 'M -5,-5 L 5,-10 L 0,0 C 200,80   -100,80   100,0 A10,10,0,0,1,100,10'; # TODO: MAKE THIS A MULTI-SEGMENT PATH WITH SMOOTH AND NON-SMOOTH SEGMENT MEETINGS. ALSO DO A PATH WITH MULTIPLE "M" COMMANDS.
my $mp = Math::MPath->newlite($mp_d,0.00001);

my @xys;
my @xys_point;
my @xys_XY;

#diag("f_offset");
#push(@xys,$mp_bez_seg->f_offset($_,1.5)) for (map {$_ * (100/500)} (0 .. 500));
#diag("point_offset");
push(@xys_point,$mp->point_offset($_,0.75)) for (map {$_ * (1/500)} (0 .. 500));
#diag("X_offset,Y_offset");
#push(@xys_XY,[$mp_bez_seg->X_offset($_,0.375),$mp_bez_seg->Y_offset($_,0.375)]) for (map {$_ * (1/500)} (0 .. 500));


my $seq = 0; # you could animate dot display?

my @pts_svg;

#push @pts_svg, map {
#              '<circle cx="'.$_->[0].'" cy="'.$_->[1].'" seq="'.$_->[2].'" r="0.1" class="bzpt"/>'
#              } map {
#                     my $x=$_->[0];
#                     my $y=$_->[1];
#                     #warn "$x,$y\n";
#                     [$x,$y,$seq++];
#                     #(map {[$x,$_,$seq++]} @{$_->[1]})
#                } @xys;

push @pts_svg, map {
              '<circle cx="'.$_->[0].'" cy="'.$_->[1].'" seq="'.$_->[2].'" r="0.1" class="bzpt2"/>'
              } map {
                     my $x=$_->[0];
                     my $y=$_->[1];
                     #warn "$x,$y\n";
                     [$x,$y,$seq++];
                     #(map {[$x,$_,$seq++]} @{$_->[1]})
                } @xys_point;

#push @pts_svg, map {
#              '<circle cx="'.$_->[0].'" cy="'.$_->[1].'" seq="'.$_->[2].'" r="0.1" class="bzpt3"/>'
#              } map {
#                     my $x=$_->[0];
#                     my $y=$_->[1];
#                     #warn "$x,$y\n";
#                     [$x,$y,$seq++];
#                     #(map {[$x,$_,$seq++]} @{$_->[1]})
#                } @xys_XY;



my $pts_svg = join("\n", @pts_svg);


my $svg = <<"EOSVG";
<svg xmlns="http://www.w3.org/2000/svg"
     width="100%"
     viewBox="0,0,100,100"
     preserveAspectRatio="xMidYMin slice"
>
<style>
/* <![CDATA[ */
.bzpt { fill:blue; stroke:none; }
.bzpt2 { fill:green; stroke:none; }
.bzpt3 { fill:red; stroke:none; }
/* ]]> */
</style>
<path d="$mp_d" stroke-width="0.03" stroke="black" fill="none" />
$pts_svg
</svg>
EOSVG

my $fn = $0.'.svg';

#diag("svg:\n$svg\n$fn\n");

# print this to an async network destination instead of just file here
# make that an "author test" that gets skipped for anyone else, however that works
open(my $tout,'>',$fn);
print $tout $svg;
close($tout);

# simple approach first
# just rsync that file to a good destination with a live display
# then trigger that destination to refresh
# if this is supposed to work on live site, display should be browser - browser polls for updated image, or have easy interface - foot switch? to reload it (to avoid all logging of each of those polls)
# put it in new "test series" "site" on local machine.
# have an index.cgi there that serves whatever images are in that destination folder - later it gets smart about pagination and formats.

my $rsync_cmd = 'rsync ./'.$fn.' /var/www/test_series/wwwroot/deviews/';

#my $rsync_result = `$rsync_cmd`;

#diag("\n$rsync_cmd :\n$rsync_result\n");

ok(1);
