#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 1;

use Math::MPath;

my $pi = atan2(1,1)*4;
sub tan {sin($_[0])/cos($_[0]);}
sub asin {
    if ($_[0] eq -1)    {return $pi/2 * -1} #eq 2
    elsif ($_[0] eq 1)  {return $pi/2     } #eq 4
    elsif ($_[0] eq 0)  {return 0;             } #eq 3
    else {                                       #eq 15
        if (abs($_[0])>1) {warn("giving something bigger than +/-1 to your asin:",$_[0],"\n");}
        return atan2($_[0],sqrt(1-$_[0]**2));
    }
}

my $mp1 = Math::MPath->newlite('M100,0 A10,10,0,0,1,100,10',0.00001);
my $mp2 = Math::MPath->newlite('M100,0 A10,10,0,1,1,100,10',0.00001);
my $mp3 = Math::MPath->newlite('M100,0 A10,10,0,1,0,100,10',0.00001);
my $mp4 = Math::MPath->newlite('M100,0 A10,10,0,0,0,100,10',0.00001);
my $mp5 = Math::MPath->newlite('M-25,0 A25,10,0,1,0,25,0',0.00001);

my @xys1 = map $mp1->point_offset($_,0.75), (map {$_ * (1/100)} (0 .. 100));
my @xys2 = map $mp2->point_offset($_,0.75), (map {$_ * (1/100)} (0 .. 100));
my @xys3 = map $mp3->point_offset($_,0.75), (map {$_ * (1/100)} (0 .. 100));
my @xys4 = map $mp4->point_offset($_,0.75), (map {$_ * (1/100)} (0 .. 100));


my $e5 = $mp5->{pathSegments}->[0];
my @xys5 = map $mp5->point($_), (map {$_ * (1/100)} (0 .. 100));
my @lines5 = map    [$mp5->point($_),[ cos($e5->theta_of_t($_)) * ($e5->{rx} - ($e5->{ry}**2/$e5->{rx})) ,0]], (map {$_ * (1/20)} (0 .. 20));

my $offdist=2;

my $initial_delta_x =             $e5->evalXofTheta(0 ) - cos($e5->theta_of_t(0 )) * ($e5->{rx} - ($e5->{ry}**2/$e5->{rx}));
my @lines51 = map [[$_ , 0],[$_,20]], (
        map {

            # this appears to be hitting the midpoint of all the tangents when offdist is 2
            # 1/4 point when 1.
            my $x0 = $e5->evalXofTheta($_);
            my $x_intercept = cos($e5->theta_of_t($_)) * ($e5->{rx} - ($e5->{ry}**2/$e5->{rx}));
            my $x_intercept_prime = -sin($e5->theta_of_t($_)) * ($e5->{rx} - ($e5->{ry}**2/$e5->{rx})); # useful?
            my $delta_x = $x0 - $x_intercept;
            my $factor = $delta_x / $initial_delta_x;
            
            my $factor2=$x_intercept/$e5->evalXofTheta(0 );

            #$x0 + $factor * $offdist;
my $theta = $e5->theta_of_t($_);
            my $xoff_of_theta = cos(atan2($e5->evalYofArcTheta($theta),$e5->evalXofArcTheta($theta) - $x_intercept)) * $offdist;

            my $theta_of_xoff = -1 * atan2($e5->{ry}, -($e5->{ry}**2/$e5->{rx}) * tan($pi/2 - asin($xoff_of_theta/$offdist)));

# This correction works sometimes, so theta_of_xoff() is probably basically correct, just needs some sign and phase fiddling.
# ... later :
# THIS IS experimentation I had to set aside in favor of a root finding approach
# but if you did really find some kind of theta(offsetx) function here, that would be nice.
# Keep this mess around a while in case you get a chance to re-figure out what this is
# and whether you can really get it to always work.
$theta_of_xoff = -$theta_of_xoff - $pi/2;

if ($e5->{delta_theta} > 0) {while ($theta_of_xoff < $e5->{theta1}) {$theta_of_xoff += 2*$pi;}}
if ($e5->{delta_theta} < 0) {while ($theta_of_xoff > $e5->{theta1}) {$theta_of_xoff -= 2*$pi;}}


#warn "eh: $theta_of_xoff eq $theta ??$e5->{theta1} $e5->{theta2}\n";
            $x0 + -$xoff_of_theta;


        }
        map {$_ * (1/20)} (0 .. 20)
);

my @xys52 = map $mp5->point_offset($_,-$offdist), (map {$_ * (1/20)} (0 .. 20));


my @xys53;
foreach my $point_offset (@xys52) {
    my $offsetx = $point_offset->[0];
    my $ry = $e5->{ry};
    my $rx = $e5->{rx};
#not working yet - go re-work it out.
    my $frac = $rx/($rx**2 - $ry**2);
    my $radical = sqrt(1 + (($ry*$frac)/($rx*$frac - 1))**2);
    my $X_cl_of_offsetx = ($offsetx - ($offdist/$radical))
                          /
                          (1 + (sqrt($rx**2 + $ry**2)/($radical*($rx - ($ry**2/$rx)))))
    ;
    push @xys53, [$X_cl_of_offsetx,0];
}

my $focus1 = $e5->{f1};
my $focus2 = $e5->{f2};

my $seq = 0; # you could animate dot display?

my @pts_svg;
my @lines_svg;

push @pts_svg, map {
              '<circle cx="'.$_->[0].'" cy="'.$_->[1].'" seq="'.$_->[3].'" r="0.1" fill="'.$_->[2].'" class=""/>'
              } map {
                     my $color=$_->[0];
                     my $pts=$_->[1];
                     map {my $x=$_->[0];
                          my $y=$_->[1];
                          [$x,$y,$color,$seq++];
                         } @$pts
                } (['#FF0000',\@xys1],
                   ['#00FF00',\@xys2],
                   ['#0000FF',\@xys3],
                   ['#FF00FF',\@xys4],
                   ['#000000',\@xys5],
                   ['#00FFFF',\@xys52],
                   ['#FFFF00',\@xys53],
                   ['#000000',[$focus1,$focus2]],
                  );

push @lines_svg, map {
              '<line x1="'.$_->[0]->[0].'" y1="'.$_->[0]->[1].'" x2="'.$_->[1]->[0].'" y2="'.$_->[1]->[1].'" seq="'.$_->[3].'" stroke-width="0.01" stroke="'.$_->[2].'" class=""/>'
              } map {
                     my $color=$_->[0];
                     my $lines=$_->[1];
                     map {my $p1=$_->[0];
                          my $p2=$_->[1];
                          [$p1,$p2,$color,$seq++];
                         } @$lines
                } (
                   ['#FF0000',\@lines5],
                   ['#000000',\@lines51],
                  );

my $pts_svg = join("\n", @pts_svg);
my $lines_svg = join("\n", @lines_svg);


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
<g transform="scale(1,-1)">
<path d="$mp1->{pathspec}" stroke-width="0.03" stroke="black" fill="none" />
<path d="$mp2->{pathspec}" stroke-width="0.03" stroke="black" fill="none" />
<path d="$mp3->{pathspec}" stroke-width="0.03" stroke="black" fill="none" />
<path d="$mp4->{pathspec}" stroke-width="0.03" stroke="black" fill="none" />
<path d="$mp5->{pathspec}" stroke-width="0.03" stroke="black" fill="none" />
$pts_svg
$lines_svg
</g>
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
