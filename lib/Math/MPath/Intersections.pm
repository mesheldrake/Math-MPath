################################################################################
###      Math::MPath::Intersections    #########################################
################################################################################
#
# Every possible SVG path segment type intersection with every other
# segment type.
#
# Letters refer to path segments as defined in SVG 1.1 spec.
#
# Need intersection functions for all these combinations:
#
# LL : LL, LH, LV, HV
# AL : A1L, A2L, A1H, A2H, A1V, A2V
# AA : A1A1, A1A2, A2A2
# QL : QL, QH, QV
# QA : QA1, QA2
# QQ
# CL : CL, CH, CV
# CA : CA1, CA2
# CQ
# CC
#
# (A1 is circular arc - special simpler case of elliptical arc)
# (A2 is general case elliptical arc)
# (Z (closePath) is equivalent to L for intersections)
#





package Math::MPath::Intersections;
{
use Math::MPath::Function::Root qw(BrentsMethod FalsePosition);
use Math::MPath::CubicFormula;
use Math::MPath::QuadraticFormula;

our $pi = 4 * atan2(1,1);

sub intersect_LL {
    my ($seg1, $seg2, $wantThetas) = @_;

    my @ret;

    my $segsegret;

    my $x1= $seg1->{p1}->[0];my $y1= $seg1->{p1}->[1];
    my $x2= $seg1->{p2}->[0];my $y2= $seg1->{p2}->[1];
    my $u1=$seg2->{p1}->[0];my $v1=$seg2->{p1}->[1];
    my $u2=$seg2->{p2}->[0];my $v2=$seg2->{p2}->[1];
    my $m1 = (($x2 - $x1)==0)?'Inf':($y2 - $y1)/($x2 - $x1);
    my $m2 = (($u2 - $u1)==0)?'Inf':($v2 - $v1)/($u2 - $u1);

    my $b1;
    my $b2;

    my  $xi;
    my $dm = $m1 - $m2;
    if    ($m1 eq 'Inf' && $m2 ne 'Inf') {$xi = $x1;$b2 = $v1 - ($m2 * $u1);}
    elsif ($m2 eq 'Inf' && $m1 ne 'Inf') {$xi = $u1;$b1 = $y1 - ($m1 * $x1);}
    elsif (abs($dm) > 0.000000000001) {
        $b1 = $y1 - ($m1 * $x1);
        $b2 = $v1 - ($m2 * $u1);
        $xi=($b2-$b1)/$dm;
    }
    my @lowhiu=($u2>$u1)?($u1,$u2):($u2,$u1);
    if ($m1 ne 'Inf') {
        my @lowhix=($x2>$x1)?($x1,$x2):($x2,$x1);
        if ($m2 eq 'Inf' &&   ($u2<$lowhix[0] || $u2>$lowhix[1]) ) {
            # NO INTERSECTION
        }
        elsif (
            ($xi || $xi eq 0) &&
            ($xi < $lowhix[1] || $xi eq $lowhix[1]) &&
            ($xi > $lowhix[0] || $xi eq $lowhix[0]) &&
            ($xi < $lowhiu[1] || $xi eq $lowhiu[1]) &&
            ($xi > $lowhiu[0] || $xi eq $lowhiu[0])
            ) {
            my $y=($m1*$xi)+$b1;
            my @lowhiv=($v2>$v1)?($v1,$v2):($v2,$v1);
            if ($m2 eq 'Inf' &&
                ($y<$lowhiv[0] || $y>$lowhiv[1])
                ) {
                # NO INTERSECTION
                # In this case we set $xi above even though there might not
                # be an intersection. If $y is not in range of the other
                # seg's y extremes, there is no intersection.
            }
            else {
                $segsegret = [$xi,$y];
            }
        }
    }
    elsif ($m2 ne 'Inf'
        && (
            ($x1 > $lowhiu[0] && $x1 < $lowhiu[1])
            ||
            ($x1 eq $lowhiu[0] && $x1 eq $lowhiu[1])
            )
        ) {
        my @lowhiy=($y2>$y1)?($y1,$y2):($y2,$y1);
        my @lowhiv=($v2>$v1)?($v1,$v2):($v2,$v1);
        my $yi = ($m2*$xi)+$b2;
        if (($yi || $yi eq 0) &&
            ($yi < $lowhiy[1] || $yi eq $lowhiy[1]) &&
            ($yi > $lowhiy[0] || $yi eq $lowhiy[0]) &&
            ($yi < $lowhiv[1] || $yi eq $lowhiv[1]) &&
            ($yi > $lowhiv[0] || $yi eq $lowhiv[0])
            ) {
            $segsegret=[$xi,$yi];
        }
    }

    if (defined($segsegret)) {
        if ($wantThetas) {push(@ret,($m1 eq 'Inf')?$seg1->solveYforTheta($segsegret->[1]):$seg1->solveXforTheta($segsegret->[0]));}
        else {push(@ret,$segsegret);}
    }
    return @ret;
}

sub intersect_AL {
    my ($arc, $line, $wantThetas, $lineIsSelf) = @_;
    my @ret;
    my @intersections;
    my $x1=$line->{p1}->[0];
    my $y1=$line->{p1}->[1];
    my $x2=$line->{p2}->[0];
    my $y2=$line->{p2}->[1];
    $x1-=$arc->{cx};
    $y1-=$arc->{cy};
    $x2-=$arc->{cx};
    $y2-=$arc->{cy};
    if (
        ($line->{maxx}>$arc->{minx}  || $line->{maxx} eq $arc->{minx})&&
        ($line->{minx}<$arc->{maxx}  || $line->{minx} eq $arc->{maxx})&&
        ($line->{maxy}>$arc->{miny} || $line->{maxy} eq $arc->{miny}) &&
        ($line->{miny}<$arc->{maxy}  || $line->{miny} eq $arc->{maxy})
        ) {
        my $rot_line_p1 = _rotate2d([0,0],[$x1,$y1],-1 * $arc->{phi_radians});
        my $rot_line_p2 = _rotate2d([0,0],[$x2,$y2],-1 * $arc->{phi_radians});

        if (abs(($rot_line_p2->[0]-$rot_line_p1->[0]))>0.000001) { # had this at > 0.1, but that seems weak. make it smaller and see if it's a problem
            my $rot_line_slope=($rot_line_p2->[1]-$rot_line_p1->[1])/($rot_line_p2->[0]-$rot_line_p1->[0]);
            my $a = (($rot_line_slope)**2/$arc->{ry}**2) + 1/$arc->{rx}**2;
            my $b = ( 2 * ($rot_line_slope) * ($rot_line_p1->[1] - ($rot_line_slope)*$rot_line_p1->[0]))/$arc->{ry}**2;
            my $c =(($rot_line_p1->[1] - ($rot_line_slope)*$rot_line_p1->[0])**2 / $arc->{ry}**2 ) - 1;
            my @xs = &quadraticformula($a,$b,$c,1);
            for (my $i=0;$i<@xs;$i++) {
                my $y=$rot_line_slope * $xs[$i] + ($rot_line_p1->[1] - $rot_line_slope * $rot_line_p1->[0]); #line formula
                push(@intersections,[$xs[$i],$y]);
            }
        }
        else { #vertical line - use ellipse formula to get points
            my $y=sqrt($arc->{ry}**2 * (1 - ($x1**2)/($arc->{rx}**2)));#vertical line. use ellipse formula to get the +/- y vals
            push(@intersections,[$x1,$y],[$x1,-$y]);
        }
        for (my $i=0;$i<@intersections;$i++) {
            # there was an interesting point of failure here due to floating point error in addition
            # when adding cx back to intersection's x, sometimes the result would get an extra 0.000000000000001 (or so),
            # when really it should have been equal to the vert line's minx/maxx
            # so the equality tests later would fail when they shouldn't have.
            # So here I'm adjusting significant digits for the addition with a sprintf(),
            # and then evaling the result to clip off any trailing zeros.
            # (So four simple lines of code become a jumbled ugly eight.)
            # (plus eight more to comment on it all)
            my $h=sqrt($intersections[$i]->[0]**2 + $intersections[$i]->[1]**2);
            $intersections[$i] = _rotate2d([0,0],$intersections[$i],$arc->{phi_radians});
            my $sigcnt1=length(($arc->{cx}             =~/\.([0-9]*)/g)[0]) + length(($arc->{cx}             =~/([0-9]*)[0-9]\./g)[0]);
            my $sigcnt2=length(($intersections[$i]->[0]=~/\.([0-9]*)/g)[0]) + length(($intersections[$i]->[0]=~/([0-9]*)[0-9]\./g)[0]);
            $intersections[$i]->[0]+=$arc->{cx};
            my $sigcnt3=length(($intersections[$i]->[0]=~/\.([0-9]*)/g)[0]);
            $intersections[$i]->[0]=0 + sprintf("%.".(sort {$b<=>$a} ($sigcnt1,$sigcnt2))[0]."f",$intersections[$i]->[0]);
            $intersections[$i]->[1]+=$arc->{cy};
        }

        #Now check to see of those intersections are within bounds and within sweep

            @intersections = grep {
            ($_->[0] < $line->{maxx} || abs($_->[0] - $line->{maxx}) < $line->{precision}) &&
            ($_->[0] > $line->{minx} || abs($_->[0] - $line->{minx}) < $line->{precision}) &&
            ($_->[1] < $line->{maxy} || abs($_->[1] - $line->{maxy}) < $line->{precision}) &&
            ($_->[1] > $line->{miny} || abs($_->[1] - $line->{miny}) < $line->{precision})
            } @intersections;

        my $leg1;
        my $leg2;
        if ($arc->{large_arc_flag}==0) {
            if ($arc->{sweep_flag} == 0) {
                $leg1=[[$arc->{cx},$arc->{cy}],[$arc->{p1}->[0],$arc->{p1}->[1]]];
                $leg2=[[$arc->{cx},$arc->{cy}],[$arc->{p2}->[0],$arc->{p2}->[1]]];
            }
            else {
                $leg1=[[$arc->{cx},$arc->{cy}],[$arc->{p2}->[0],$arc->{p2}->[1]]];
                $leg2=[[$arc->{cx},$arc->{cy}],[$arc->{p1}->[0],$arc->{p1}->[1]]];
            }
        }
        else {
            if ($arc->{sweep_flag} == 0) {
                $leg1=[[$arc->{cx},$arc->{cy}],[$arc->{p2}->[0],$arc->{p2}->[1]]];
                $leg2=[[$arc->{cx},$arc->{cy}],[$arc->{p1}->[0],$arc->{p1}->[1]]];
            }
            else {
                $leg1=[[$arc->{cx},$arc->{cy}],[$arc->{p1}->[0],$arc->{p1}->[1]]];
                $leg2=[[$arc->{cx},$arc->{cy}],[$arc->{p2}->[0],$arc->{p2}->[1]]];
            }
        }
        @intersections = grep {
                                  ( $arc->{large_arc_flag} && !$arc->isWithinSweep($_,$leg1,$leg2))
                               || (!$arc->{large_arc_flag} &&  $arc->isWithinSweep($_,$leg1,$leg2))
                         } @intersections;

        if ($wantThetas) {
            foreach my $int (@intersections) {
                if ($lineIsSelf) {
                    push(@ret,($line->{m} eq 'inf')?$line->solveYforTheta($int->[1]):$line->solveXforTheta($int->[0]));
                }
                else {
                    my @allArcThetas=$arc->solveXforTheta($int->[0]);
                    foreach my $t (@allArcThetas) {
                        my $tp=$arc->point($t);
                        if (abs($tp->[1] - $int->[1]) < 0.0000000001) {push(@ret,$t);}
                    }
                }
            }
        }
        else {
            push(@ret,@intersections);
        }
    }

    return @ret;
}

# optimized circle-circle intersection case
# don't know whether it's actually better than the general elliptical arc case
# at this point, but should revisit this math and improve with new features
# of elliptical arc code if possible
sub intersect_A1A1 {
    my ($arc1, $arc2, $wantThetas) = @_;

    my @ret;

    # Paul Bourke
    # http://paulbourke.net/geometry/circlesphere/

    my $center_dist = sqrt(($arc2->{cx}-$arc1->{cx})**2 + ($arc2->{cy}-$arc1->{cy})**2);
    my $radius_sum  = $arc1->{rx} + $arc2->{rx};
    next if ($center_dist > $radius_sum);
    my $radius_diff  = $arc1->{rx} - $arc2->{rx};
    next if ($center_dist < abs($radius_diff));

    if ($center_dist eq 0 && $arc1->{rx} eq $arc2->{rx}) {
        die "crap. come figure out infinite intersection solution for circle-circle arc intersection";
        }

    my $a = ($arc1->{rx}**2 - $arc2->{rx}**2 + $center_dist**2)/(2*$center_dist);
    my $h = sqrt($arc1->{rx}**2 - $a**2);
    my $x2 = $arc1->{cx} + $a*($arc2->{cx} - $arc1->{cx})/$center_dist;
    my $y2 = $arc1->{cy} + $a*($arc2->{cy} - $arc1->{cy})/$center_dist;

    my @intersections;
    push @intersections, [$x2 + $h * ($arc2->{cy} - $arc1->{cy}) / $center_dist,
                          $y2 - $h * ($arc2->{cx} - $arc1->{cx}) / $center_dist];
    # unless the at-one-point intersect case, figure the other intersection
    unless (   $center_dist eq $arc1->{rx} + $arc2->{rx}
            || $center_dist eq $arc1->{rx} - $arc2->{rx}
           ) {
    push @intersections, [$x2 - $h * ($arc2->{cy} - $arc1->{cy}) / $center_dist,
                          $y2 + $h * ($arc2->{cx} - $arc1->{cx}) / $center_dist];
        }

    # now, see if those are within arc sweep
    # hopefully just copy paste from arc-line code...

    if (scalar(@intersections) > 0) {
        my $leg1;
        my $leg2;

        # copy pasted modded twice, once for each arc in this case
        if ($arc1->{large_arc_flag}==0) {
            if ($arc1->{sweep_flag} == 0) {
                $leg1=[[$arc1->{cx},$arc1->{cy}],[$arc1->{p1}->[0],$arc1->{p1}->[1]]];
                $leg2=[[$arc1->{cx},$arc1->{cy}],[$arc1->{p2}->[0],$arc1->{p2}->[1]]];
                }
            else {
                $leg1=[[$arc1->{cx},$arc1->{cy}],[$arc1->{p2}->[0],$arc1->{p2}->[1]]];
                $leg2=[[$arc1->{cx},$arc1->{cy}],[$arc1->{p1}->[0],$arc1->{p1}->[1]]];
                }
            }
        else {
            if ($arc1->{sweep_flag} == 0) {
                $leg1=[[$arc1->{cx},$arc1->{cy}],[$arc1->{p2}->[0],$arc1->{p2}->[1]]];
                $leg2=[[$arc1->{cx},$arc1->{cy}],[$arc1->{p1}->[0],$arc1->{p1}->[1]]];
                }
            else {
                $leg1=[[$arc1->{cx},$arc1->{cy}],[$arc1->{p1}->[0],$arc1->{p1}->[1]]];
                $leg2=[[$arc1->{cx},$arc1->{cy}],[$arc1->{p2}->[0],$arc1->{p2}->[1]]];
                }
            }
        @intersections = grep {
            (     $arc1->{large_arc_flag} && !$arc1->isWithinSweep($_,$leg1,$leg2))
             || (!$arc1->{large_arc_flag} &&  $arc1->isWithinSweep($_,$leg1,$leg2))
            } @intersections;

        # now for other arc
        if ($arc2->{large_arc_flag}==0) {
            if ($arc2->{sweep_flag} == 0) {
                $leg1=[[$arc2->{cx},$arc2->{cy}],[$arc2->{p1}->[0],$arc2->{p1}->[1]]];
                $leg2=[[$arc2->{cx},$arc2->{cy}],[$arc2->{p2}->[0],$arc2->{p2}->[1]]];
                }
            else {
                $leg1=[[$arc2->{cx},$arc2->{cy}],[$arc2->{p2}->[0],$arc2->{p2}->[1]]];
                $leg2=[[$arc2->{cx},$arc2->{cy}],[$arc2->{p1}->[0],$arc2->{p1}->[1]]];
                }
            }
        else {
            if ($arc2->{sweep_flag} == 0) {
                $leg1=[[$arc2->{cx},$arc2->{cy}],[$arc2->{p2}->[0],$arc2->{p2}->[1]]];
                $leg2=[[$arc2->{cx},$arc2->{cy}],[$arc2->{p1}->[0],$arc2->{p1}->[1]]];
                }
            else {
                $leg1=[[$arc2->{cx},$arc2->{cy}],[$arc2->{p1}->[0],$arc2->{p1}->[1]]];
                $leg2=[[$arc2->{cx},$arc2->{cy}],[$arc2->{p2}->[0],$arc2->{p2}->[1]]];
                }
            }
        @intersections = grep {
            (     $arc2->{large_arc_flag} && !$arc2->isWithinSweep($_,$leg1,$leg2))
             || (!$arc2->{large_arc_flag} &&  $arc2->isWithinSweep($_,$leg1,$leg2))
            } @intersections;

        if ($wantThetas) {
            foreach my $int (@intersections) {
                my @allArcThetas=$arc1->solveXforTheta($int->[0]);
                foreach my $t (@allArcThetas) {
                    my $tp=$arc1->point($t);
                    if (abs($tp->[1] - $int->[1]) < 0.0000000001) {push(@ret,$t);}
                    }
                }
            }
        else {
            push(@ret,@intersections);
            }
        }

    return @ret;

}

sub intersect_CL {
    my ($curve, $line, $wantThetas, $lineIsSelf) = @_;
    my @ret;

    # t^3 + [(F-mB)/(E-mA)]t^2 + [(G-mC)/(E-mA)]t + [(H-mD-x0)/(E-mA)] = 0

    # You took the equation of a line y = mx + b
    # Used the Bezier's y(t) and x(t) for y and x, used the line's m and b,
    # then arranged it like this: y - mx - b = 0
    # and the left side comes out as a cubic polynomial in t that we can
    # use the good old cubic solver for.

    # First some special cases for vertical and horizontal lines.
    my @thetas;
    if ($line->{m} eq 'inf' || $line->{m} eq '-inf') {
        @thetas = $curve->solveXforTheta($line->{maxx});
        foreach my $t (@thetas) {
            my $y = $curve->bezierEvalYofT($t);
            if (($y < $line->{maxy} || $y eq $line->{maxy}) &&
                ($y > $line->{miny} || $y eq $line->{miny})) {
                if ($wantThetas) {
                    if (!$lineIsSelf) {push(@ret,$t);}
                    else {push(@ret,$line->solveYforTheta($y));}
                }
                else {push(@ret,[$line->{maxx},$y]);}
            }
        }
    }
    elsif ($line->{m} eq 0) {
        if ($wantThetas && !$lineIsSelf) {
            my @ths = $curve->solveYforTheta($line->{p1}->[1]);
            push(@ret, grep {my $p=$curve->point($_); ($p->[0] < $line->{maxx} || $p->[0] eq $line->{maxx}) && ($p->[0] > $line->{minx} || $p->[0] eq $line->{minx})} @ths);
        }
        else {
            # do F(y) to get possible x vals from bezier
            my @xs = $curve->F($line->{p1}->[1]);
            # then filter to what's actually within line segment bounds
            # also, the reason this zero slope thing is handled seperately
            # is that we make sure the resulting y values are exactly the horizontal line's y value
            if ($wantThetas) {
                push(@ret,map { $line->solveXforTheta($_)} @xs);
            }
            else {
                my @intersections = map {[$_,$line->{p1}->[1]]} grep { ($_ < $line->{maxx} || $_ eq $line->{maxx}) && ($_ > $line->{minx} || $_ eq $line->{minx})} @xs;
                push(@ret,@intersections);
            }
        }
    }
    else {

        @thetas = &cubicformula(
            ($curve->{F}-$line->{m}*$curve->{B})            / ($curve->{E}-$line->{m}*$curve->{A}),
            ($curve->{G}-$line->{m}*$curve->{C})            / ($curve->{E}-$line->{m}*$curve->{A}),
            ($curve->{H}-$line->{m}*$curve->{D}-$line->{b}) / ($curve->{E}-$line->{m}*$curve->{A}),
            1);

        @thetas = sort {$a<=>$b} grep {(1 > $_ || 1 eq $_) && ($_ > 0 || $_ eq 0)}  @thetas;

        foreach my $t (@thetas) {
            my $x = $curve->bezierEvalXofT($t);
            if (($x < $line->{maxx} || $x eq $line->{maxx}) &&
                ($x > $line->{minx} || $x eq $line->{minx})) {
                if ($wantThetas) {
                    if (!$lineIsSelf) {push(@ret,$t);}
                    else {push(@ret,$line->solveXforTheta($x));}
                }
                else {push(@ret,[$x,$curve->bezierEvalYofT($t)]);}
            }
        }
    }
    return @ret;
}
sub intersect_CC {

    # returns list of intersections and the corresponding Bezier parameters
    # for each input Bezier.

    # lutA and lutB optional - used to re-call this function in the intersection
    # pair case, where we prepare LUTs to isolate each of the pair

    my ($bezA, $bezB, $lutA, $lutB) = @_;

    my $XtoTLUT_A = $lutA ? $lutA : $bezA->{XtoTLUT};
    my $XtoTLUT_B = $lutB ? $lutB : $bezB->{XtoTLUT};

    my @ret;

    my %aseen;

    foreach my $spanA (@{$XtoTLUT_A}) {

        $aseen{$spanA}={} if !$aseen{$spanA};

        foreach my $spanB (@{$XtoTLUT_B}) {

        # avoid duplicate runs when looking for self intersection
        next if $spanA == $spanB;
        next if $aseen{$spanA}->{$spanB};
        $aseen{$spanB}->{$spanA}++;

        next if ($spanA->[1]->[0] > $spanB->[1]->[-1]);
        next if ($spanB->[1]->[0] > $spanA->[1]->[-1]);

        my $spanx = [$spanA->[1]->[0] > $spanB->[1]->[0] ? $spanA->[1]->[0] : $spanB->[1]->[0], $spanA->[1]->[-1] < $spanB->[1]->[-1] ? $spanA->[1]->[-1] : $spanB->[1]->[-1]];

        my $tsA = [$spanA->[0]->[0]->($spanx->[0]), $spanA->[0]->[0]->($spanx->[1])];
        my $tsB = [$spanB->[0]->[0]->($spanx->[0]), $spanB->[0]->[0]->($spanx->[1])];

        my $ysA = [$bezA->bezierEvalYofT($tsA->[0]),$bezA->bezierEvalYofT($tsA->[1])];
        my $ysB = [$bezB->bezierEvalYofT($tsB->[0]),$bezB->bezierEvalYofT($tsB->[1])];

        #warn "xspan:\n $spanx->[0],$spanx->[1]\n";
        #warn "ts:\n $tsA->[0],$tsA->[1]\n $tsB->[0],$tsB->[1]\n";
        #warn "ys:\n $ysA->[0],$ysA->[1]\n $ysB->[0],$ysB->[1]\n";

        # one intersection case, similar to how you'd test for crossing line segments
        if    (($ysA->[0] > $ysB->[0] && $ysA->[0] ne $ysB->[0] &&
                $ysA->[1] < $ysB->[1] && $ysA->[1] ne $ysB->[1]   ) ||
               ($ysA->[0] < $ysB->[0] && $ysA->[0] ne $ysB->[0] &&
                $ysA->[1] > $ysB->[1] && $ysA->[1] ne $ysB->[1]   )
              ) {

            my $findintbezbez = sub {
                # Y1(t1) - Y2(t2(X1(t1))) secret sauce ingredient
                return $bezA->bezierEvalYofT($_[0]) - $bezB->bezierEvalYofT( $spanB->[0]->[0]->( $bezA->bezierEvalXofT($_[0]) ) );
            };
            my $bounds_tA = ($tsA->[0] < $tsA->[1]) ? [$tsA->[0],$tsA->[1]] : [$tsA->[1],$tsA->[0]];
            my ($int_t_A,$msg)=BrentsMethod($findintbezbez,$bounds_tA,0.00001,undef,'subBez-subBez intersection finding');

            my $intersection_x = $bezA->bezierEvalXofT($int_t_A);
            my $intersection_y = $bezA->bezierEvalYofT($int_t_A);

            $int_t_B = $spanB->[0]->[0]->($intersection_x);

            push @ret, [$intersection_x,$intersection_y,$int_t_A,$int_t_B];
        }
        # catch any endpoint overlap "intersections"
        elsif (   $ysA->[0] eq $ysB->[0] || $ysA->[1] eq $ysB->[1]
               || $ysA->[1] eq $ysB->[0] || $ysA->[0] eq $ysB->[1]
              ) {
            # what is the policy on this?
            # when checking for self intersections, def don't want these
            # but for diff bezs, these could sometimes be legit intersections
            #push @ret, ["xoverlap","yoverlap"] if $bezA != $bezB;
        }
        # the zero or two intersection case
        else {

            warn "failed to find one intersection when a special LUT was provided\n" if defined($lutA) || defined($lutB);
            next if defined($lutA) || defined($lutB);

            # zero intersections if y ranges don't overlap
            my ($lowyA,$highyA) = ($ysA->[0]<$ysA->[1]) ? (@$ysA) : (reverse @$ysA);
            my ($lowyB,$highyB) = ($ysB->[0]<$ysB->[1]) ? (@$ysB) : (reverse @$ysB);
            next if ($lowyA > $highyB);
            next if ($lowyB > $highyA);

            # zero or two intersections

            # more secret sauce
            # if g(t1) = [Y1(t1) - Y2(t2(X1(t1))) ]'
            # crosses zero over the t1 parameter range [$tsA->[0],$tsA->[1]]
            # then there are two intersections in that range.
            # We can then root find to find the parameter for that zero crossing,
            # and use that to split the parameter range, guaranteeing one
            # intersection on each of the two new split parameter ranges.
            # We can then re-call this intersection subroutine, specifying those
            # new ranges to operate on.

            my $y_diff_prime = sub {
                return
                $bezA->bezierEvalYPrimeofT($_[0])
                -
                $bezB->bezierEvalYPrimeofT( $spanB->[0]->[0]->( $bezA->bezierEvalXofT($_[0]) ) )
                *
                $spanB->[0]->[1]->($bezA->bezierEvalXofT($_[0]))
                *
                $bezA->bezierEvalXPrimeofT($_[0]);
            };

            my $at_start = $y_diff_prime->($tsA->[0]);
            my $at_end   = $y_diff_prime->($tsA->[1]);

            #warn "ydiffprime start, end: $at_start, $at_end\n";

            if ($at_start > 0 && $at_end < 0 || $at_start < 0 && $at_end > 0) {

                my $bounds_tA = ($tsA->[0] < $tsA->[1]) ? [$tsA->[0],$tsA->[1]] : [$tsA->[1],$tsA->[0]];
                my ($split_t_A,$msg)=BrentsMethod($y_diff_prime,$bounds_tA,0.00001,undef,'subBez-subBez intersection finding - find pair split parameter');

                my $span_split_x = $bezA->bezierEvalXofT($split_t_A);

                $split_t_B = $spanB->[0]->[0]->($span_split_x);

                # set up new LUTs to pass to a re-call of this intersection sub
                # so re-call run will have proper x bounds to find each of
                # the intersection pair individually.

                my $sub_t_span_1_A = [$tsA->[0], $split_t_A];
                my $sub_t_span_2_A = [$split_t_A, $tsA->[1]];
                my $sub_t_span_1_B = [$tsB->[0], $split_t_B];
                my $sub_t_span_2_B = [$split_t_B, $tsB->[1]];

                my @int1 = intersect_CC($bezA,$bezB,
                                             [[$spanA->[0],[$spanx->[0]  , $span_split_x], $spanA->[3]?$sub_t_span_2_A:$sub_t_span_1_A, $spanA->[3]]],
                                             [[$spanB->[0],[$spanx->[0]  , $span_split_x], $spanB->[3]?$sub_t_span_2_B:$sub_t_span_1_B, $spanB->[3]]]
                                             );

                my @int2 = intersect_CC($bezA,$bezB,
                                             [[$spanA->[0],[$span_split_x, $spanx->[1]  ], $spanA->[3]?$sub_t_span_1_A:$sub_t_span_2_A, $spanA->[3]]],
                                             [[$spanB->[0],[$span_split_x, $spanx->[1]  ], $spanB->[3]?$sub_t_span_1_B:$sub_t_span_2_B, $spanB->[3]]]
                                             );

                push @ret, @int1, @int2;
            }
            else {
                #warn "no pair";
            }

        }

        }
    }

    return @ret;
}




################################################################################
### OFFSET INTERSECTIONS                                                     ###
################################################################################

sub intersect_LoLo {
    my ($lineA, $lineB, $offA, $offB, $wantThetas) = @_;

    my $lineA = Math::MPath::LineSegment->new(
        [ $lineA->{p1}->[0] + $offA * cos($lineA->{angleNormal}),
          $lineA->{p1}->[1] + $offA * sin($lineA->{angleNormal}) ],
        [ $lineA->{p2}->[0] + $offA * cos($lineA->{angleNormal}),
          $lineA->{p2}->[1] + $offA * sin($lineA->{angleNormal}) ]
    );
    my $lineB = Math::MPath::LineSegment->new(
        [ $lineB->{p1}->[0] + $offB * cos($lineB->{angleNormal}),
          $lineB->{p1}->[1] + $offB * sin($lineB->{angleNormal}) ],
        [ $lineB->{p2}->[0] + $offB * cos($lineB->{angleNormal}),
          $lineB->{p2}->[1] + $offB * sin($lineB->{angleNormal}) ]
    );

    return intersect_LL($lineA, $lineB, $wantThetas);
}

sub intersect_A1oLo {
    my ($arc, $line, $offCircle, $offLine, $wantThetas, $lineIsSelf) = @_;

    my $arc_off_radius = $arc->{rx} + $offCircle;
    my $arc_off = Math::MPath::EllipticalArc->new(
        [$arc->{p1}->[0], $arc->{p1}->[1]],
        [$arc_off_radius,$arc_off_radius],
        $arc->{phi},
        $arc->{large_arc_flag},
        $arc->{sweep_flag},
        [$arc->{p2}->[0], $arc->{p2}->[1]],
        $arc->{precision},
        $arc->{isLite}
    );

    my $line_off = Math::MPath::LineSegment->new(
        [ $line->{p1}->[0] + $offLine * cos($line->{angleNormal}),
          $line->{p1}->[1] + $offLine * sin($line->{angleNormal}) ],
        [ $line->{p2}->[0] + $offLine * cos($line->{angleNormal}),
          $line->{p2}->[1] + $offLine * sin($line->{angleNormal}) ]
    );

    return intersect_AL($arc_off, $line_off, $wantThetas, $lineIsSelf);
}

sub intersect_AoAo {

    # Adaptation of CoCo (offset Bezier curve intersections), 
    # for offset elliptical arcs.
    
    # lutA and lutB optional - used to re-call this function in the
    # two intersection case

    my ($arcA, $arcB, $offA, $offB, $lutA, $lutB) = @_;

    my $XtoTLUT_A = $lutA ? $lutA : [map [$_->[0],[$_->[1]->[0],$_->[1]->[-1]],[$_->[2]->[0],$_->[2]->[-1]],$_->[3]], @{$arcA->{XtoTLUT}}];
    my $XtoTLUT_B = $lutB ? $lutB : [map [$_->[0],[$_->[1]->[0],$_->[1]->[-1]],[$_->[2]->[0],$_->[2]->[-1]],$_->[3]], @{$arcB->{XtoTLUT}}];

    # assume any provided LUTs already contain offset Xs in their x range entries
    # otherwise we need to do that now for the LUT copies we made
    if (!defined($lutA)) {
        foreach my $span (@{$XtoTLUT_A}) {
            $span->[1]->[0]  = $arcA->X_offset($span->[2]->[0] ,$offA,$span->[0]->[1],$span->[3]);
            $span->[1]->[-1] = $arcA->X_offset($span->[2]->[-1],$offA,$span->[0]->[1],$span->[3]);
        }
    }

    if (!defined($lutB)) {
        foreach my $span (@{$XtoTLUT_B}) {
            $span->[1]->[0]  = $arcB->X_offset($span->[2]->[0] ,$offB,$span->[0]->[1],$span->[3]);
            $span->[1]->[-1] = $arcB->X_offset($span->[2]->[-1],$offB,$span->[0]->[1],$span->[3]);
        }
    }

    my @ret;

    #my %aseen;

    my $offintloopcnt=0;

    foreach my $spanA (@{$XtoTLUT_A}) {

        #$aseen{$spanA}={} if !$aseen{$spanA};

        foreach my $spanB (@{$XtoTLUT_B}) {

            # avoid duplicate runs
            # TODO redo this logic

            # skip if x spans don't overlap
            next if ($spanA->[1]->[0] > $spanB->[1]->[-1]);
            next if ($spanB->[1]->[0] > $spanA->[1]->[-1]);

            # x clipping to minimum x span that's valid for both lut entry spans

            my $spanx = [];
            my $tsA = [];
            my $tsB = [];

            # the x clipping is more involved with these offset curves

            if ($spanA->[1]->[0] > $spanB->[1]->[0]) {
                $spanx->[0] = $spanA->[1]->[0];
                $tsA->[0] = $spanA->[2]->[0];
                $tsB->[0] = $arcB->t_from_xoff($spanx->[0],$offB,[$spanB->[2]->[0],$spanB->[2]->[-1]],$spanB->[0]->[1]);

            }
            else {
                $spanx->[0] = $spanB->[1]->[0];
                $tsB->[0] = $spanB->[2]->[0];
                $tsA->[0] = $arcA->t_from_xoff($spanx->[0],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1]);

            }
            if ($spanA->[1]->[-1] < $spanB->[1]->[-1]) {
                $spanx->[1] = $spanA->[1]->[-1];
                $tsA->[1] = $spanA->[2]->[-1];
                $tsB->[1] = $arcB->t_from_xoff($spanx->[1],$offB,[$spanB->[2]->[0],$spanB->[2]->[-1]],$spanB->[0]->[1]);

            }
            else {
                $spanx->[1] = $spanB->[1]->[-1];
                $tsB->[1] = $spanB->[2]->[-1];
                $tsA->[1] = $arcA->t_from_xoff($spanx->[1],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1]);

            }

            my $ysA = [$arcA->Y_offset($tsA->[0],$offA,$spanA->[0]->[1],$spanA->[3]),$arcA->Y_offset($tsA->[1],$offA,$spanA->[0]->[1],$spanA->[3])];
            my $ysB = [$arcB->Y_offset($tsB->[0],$offB,$spanB->[0]->[1],$spanB->[3]),$arcB->Y_offset($tsB->[1],$offB,$spanB->[0]->[1],$spanB->[3])];

            #warn "\nspanx:\n $spanx->[0],$spanx->[1]\n";
            #warn "ts:\n $tsA->[0],$tsA->[1]\n $tsB->[0],$tsB->[1]\n";
            #warn "ys:\n $ysA->[0],$ysA->[1]\n $ysB->[0],$ysB->[1]\n";

            # one intersection case, similar to how you'd test for crossing line segments
            if    (($ysA->[0] > $ysB->[0] && $ysA->[0] ne $ysB->[0] &&
                    $ysA->[1] < $ysB->[1] && $ysA->[1] ne $ysB->[1]   ) ||
                   ($ysA->[0] < $ysB->[0] && $ysA->[0] ne $ysB->[0] &&
                    $ysA->[1] > $ysB->[1] && $ysA->[1] ne $ysB->[1]   )
                  ) {

                my $findintarcarc = sub {
                    # need Yoff1(t1(xoff)) - Yoff2(t2(xoff))
                    my $t1 = $arcA->t_from_xoff($_[0],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1]);
                    my $t2 = $arcB->t_from_xoff($_[0],$offB,[$spanB->[2]->[0],$spanB->[2]->[-1]],$spanB->[0]->[1]);

                    my $yoff1 = $arcA->Y_offset($t1,$offA,$spanA->[0]->[1],$spanA->[3]);
                    my $yoff2 = $arcB->Y_offset($t2,$offB,$spanB->[0]->[1],$spanB->[3]);

                    return $yoff2 - $yoff1;
                };

                my $bounds_xoff = [$spanx->[0],$spanx->[1]];
                my ($int_xoff,$msg)=BrentsMethod($findintarcarc,$bounds_xoff,0.00001,undef,'subArcOff-subArcOff intersection finding');

                if ($msg) { warn "arcoffint message: $msg\n"; }
                else {
                    my $t1 = $arcA->t_from_xoff($int_xoff,$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1]);
                    my $intersection_x = $arcA->X_offset($t1,$offA,$spanA->[0]->[1],$spanA->[3]);
                    my $intersection_y = $arcA->Y_offset($t1,$offA,$spanA->[0]->[1],$spanA->[3]);
                    my $x = $arcA->evalXofTheta($t1);
                    $int_t_A = $spanA->[0]->[0]->($x);
                    $int_t_B = $spanB->[0]->[0]->($x);
                    #warn "GOT THE ONE INTERSECTION: [$intersection_x,$intersection_y,$int_t_A,$int_t_B]\n";
                    push @ret, [$intersection_x,$intersection_y,$int_t_A,$int_t_B];
                }

            }
            # catch any endpoint overlap "intersections"
            elsif (    $ysA->[0] eq $ysB->[0] || $ysA->[1] eq $ysB->[1]
                    || $ysA->[1] eq $ysB->[0] || $ysA->[0] eq $ysB->[1]
                  ) {
                # what is the policy on this?
                # when checking for self intersections, def don't want these
                # but for diff bezs, these could sometimes be legit intersections
                #push @ret, ["xoverlap","yoverlap"] if $bezA != $bezB;
            }
            # the zero or two intersection case
            else {

                warn "failed to find one offset intersection when a special LUT was provided\n" if defined($lutA) || defined($lutB);
                next if defined($lutA) || defined($lutB);

                # zero intersections if y ranges don't overlap
                my ($lowyA,$highyA) = ($ysA->[0]<$ysA->[1]) ? (@$ysA) : (reverse @$ysA);
                my ($lowyB,$highyB) = ($ysB->[0]<$ysB->[1]) ? (@$ysB) : (reverse @$ysB);
                next if ($lowyA > $highyB);
                next if ($lowyB > $highyA);

                # zero or two intersections

                # Similar to offset Bezier approach
                #
                # root finding function to test for zero or two, then find the split point for two
                #
                # g(t1) = [offsetY2(t2_of_xoff(xoff)) - offsetY1(t1_of_xoff(xoff))]'

                my $y_diff_prime = sub {

                    # need to adapt this for offset elliptical arcs

                    my $tA = $arcA->t_from_xoff($_[0],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1]);
                    my $tB = $arcB->t_from_xoff($_[0],$offB,[$spanB->[2]->[0],$spanB->[2]->[-1]],$spanB->[0]->[1]);

die "CURRENTLY WORKING ON AoAo ZERO OR TWO INTERSECTION CASE";

#                    my $nA = $arcA->F_prime(undef,$tA);
#                    my $nB = $arcB->F_prime(undef,$tB);

#                    my $nA_prime = $arcA->F_2prime(undef,$tA);
#                    my $nB_prime = $arcB->F_2prime(undef,$tB);

#                    my $YPrimeA = $arcA->bezierEvalYPrimeofT($tA);
#                    my $YPrimeB = $arcB->bezierEvalYPrimeofT($tB);

#                    my $yoffset_primeA = -($offA/2.0) * 1/(sqrt($nA**2 + 1)**3) * 2*$nA * $nA_prime;
#                    my $yoffset_primeB = -($offB/2.0) * 1/(sqrt($nB**2 + 1)**3) * 2*$nB * $nB_prime;

#                    my $YoffPrimeA = ($YPrimeA + $yoffset_primeA);
#                    my $YoffPrimeB = ($YPrimeB + $yoffset_primeB);

                    my $ret = $YoffPrimeB - $YoffPrimeA;

                    return $ret;

                };

                my $at_start = $y_diff_prime->($spanx->[0]);
                my $at_end   = $y_diff_prime->($spanx->[1]);

                #warn "ydiffprime start, end: [$spanx->[0]] $at_start, [$spanx->[1]] $at_end\n";

                if ($at_start > 0 && $at_end < 0 || $at_start < 0 && $at_end > 0) {

                    my $bounds_xoffA = [$spanx->[0],$spanx->[1]];
                    my ($split_xoff_A,$msg)=FalsePosition($y_diff_prime,$bounds_xoffA,0.00001,(($bounds_xoffA->[1]-$bounds_xoffA->[0])/2),'subArc-subArc intersection finding - find pair split parameter');

                    warn "split find fail msg: $msg\n" if $msg;

                    my $span_split_x = $split_xoff_A;

                    #warn "split xoff: $span_split_x\n";

                    $split_t_A = $arcA->t_from_xoff($span_split_x,$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1]);
                    $split_t_B = $arcB->t_from_xoff($span_split_x,$offB,[$spanB->[2]->[0],$spanB->[2]->[-1]],$spanB->[0]->[1]);

                    # set up new LUTs to pass to a re-call of this intersection sub
                    # so re-call run will have proper x bounds to find each of
                    # the intersection pair individually.

                    my $sub_t_span_1_A = [$tsA->[0], $split_t_A];
                    my $sub_t_span_2_A = [$split_t_A, $tsA->[1]];
                    my $sub_t_span_1_B = [$tsB->[0], $split_t_B];
                    my $sub_t_span_2_B = [$split_t_B, $tsB->[1]];

                    #warn "re-call 1\n";
                    #warn "  [[$spanA->[0],[$spanx->[0]  , $span_split_x], $spanA->[3]?$sub_t_span_2_A:$sub_t_span_1_A, $spanA->[3]]]\n";
                    #warn "  [[$spanB->[0],[$spanx->[0]  , $span_split_x], $spanB->[3]?$sub_t_span_2_B:$sub_t_span_1_B, $spanB->[3]]]\n";
                    my @int1 = intersect_AoAo($arcA,$arcB,$offA,$offB,
                                                 [[$spanA->[0],[$spanx->[0]  , $span_split_x], $spanA->[3]?$sub_t_span_2_A:$sub_t_span_1_A, $spanA->[3]]],
                                                 [[$spanB->[0],[$spanx->[0]  , $span_split_x], $spanB->[3]?$sub_t_span_2_B:$sub_t_span_1_B, $spanB->[3]]]
                                                 );
                    #warn "re-call 2\n";
                    #warn "  [[$spanA->[0],[$span_split_x, $spanx->[1]  ], $spanA->[3]?$sub_t_span_1_A:$sub_t_span_2_A, $spanA->[3]]]\n";
                    #warn "  [[$spanB->[0],[$span_split_x, $spanx->[1]  ], $spanB->[3]?$sub_t_span_1_B:$sub_t_span_2_B, $spanB->[3]]]\n";

                    my @int2 = intersect_AoAo($arcA,$arcB,$offA,$offB,
                                                 [[$spanA->[0],[$span_split_x, $spanx->[1]  ], $spanA->[3]?$sub_t_span_1_A:$sub_t_span_2_A, $spanA->[3]]],
                                                 [[$spanB->[0],[$span_split_x, $spanx->[1]  ], $spanB->[3]?$sub_t_span_1_B:$sub_t_span_2_B, $spanB->[3]]]
                                                 );

                    push @ret, @int1, @int2;
                }
                else {
                    #warn "no pair";
                }

            }
        }
    }
    return @ret;
}

# This is the first-worked-out reference for how to do offset curve intersections
# when you have curves decomposed into monotonic subsections, with one-to-one
# t(x) functions for each subsection.
sub intersect_CoCo {

    # The first attempt at adapting bez_bez_intersect() for offset beziers

    # returns list of intersections and the corresponding Bezier parameters
    # for each input Bezier offset curve.

    # lutA and lutB optional - used to re-call this function in the
    # two intersection case

    my ($bezA, $bezB, $offA, $offB, $lutA, $lutB) = @_;

    my $XtoTLUT_A = $lutA ? $lutA : [map [$_->[0],[$_->[1]->[0],$_->[1]->[-1]],[$_->[2]->[0],$_->[2]->[-1]],$_->[3]], @{$bezA->{XtoTLUT}}];
    my $XtoTLUT_B = $lutB ? $lutB : [map [$_->[0],[$_->[1]->[0],$_->[1]->[-1]],[$_->[2]->[0],$_->[2]->[-1]],$_->[3]], @{$bezB->{XtoTLUT}}];

    # assume any provided LUTs already contain offset Xs in their x range entries
    # otherwise we need to do that now for the LUT copies we made
    if (!defined($lutA)) {
        foreach my $span (@{$XtoTLUT_A}) {
            $span->[1]->[0]  = $bezA->X_offset($span->[2]->[0] ,$offA,$span->[0]->[1],$span->[3]);
            $span->[1]->[-1] = $bezA->X_offset($span->[2]->[-1],$offA,$span->[0]->[1],$span->[3]);
        }
    }

    if (!defined($lutB)) {
        foreach my $span (@{$XtoTLUT_B}) {
            $span->[1]->[0]  = $bezB->X_offset($span->[2]->[0] ,$offB,$span->[0]->[1],$span->[3]);
            $span->[1]->[-1] = $bezB->X_offset($span->[2]->[-1],$offB,$span->[0]->[1],$span->[3]);
        }
    }

    my @ret;

    my %aseen;

    my $offintloopcnt=0;

    foreach my $spanA (@{$XtoTLUT_A}) {

        $aseen{$spanA}={} if !$aseen{$spanA};

        foreach my $spanB (@{$XtoTLUT_B}) {

            # avoid duplicate runs when looking for self intersection
            next if $spanA == $spanB && ($offA==$offB);
            next if $aseen{$spanA}->{$spanB} && ($offA==$offB);
            $aseen{$spanB}->{$spanA}++;

            # skip if x spans don't overlap
            next if ($spanA->[1]->[0] > $spanB->[1]->[-1]);
            next if ($spanB->[1]->[0] > $spanA->[1]->[-1]);

            # x clipping to minimum x span that's valid for both lut entry spans

            my $spanx = [];
            my $tsA = [];
            my $tsB = [];

            # the x clipping is more involved with these offset curves

            if ($spanA->[1]->[0] > $spanB->[1]->[0]) {
                $spanx->[0] = $spanA->[1]->[0];
                $tsA->[0] = $spanA->[2]->[0];
                $tsB->[0] = $bezB->t_from_xoff($spanx->[0],$offB,[$spanB->[2]->[0],$spanB->[2]->[-1]],$spanB->[0]->[1],$spanB->[3]);

            }
            else {
                $spanx->[0] = $spanB->[1]->[0];
                $tsB->[0] = $spanB->[2]->[0];
                $tsA->[0] = $bezA->t_from_xoff($spanx->[0],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);

            }
            if ($spanA->[1]->[-1] < $spanB->[1]->[-1]) {
                $spanx->[1] = $spanA->[1]->[-1];
                $tsA->[1] = $spanA->[2]->[-1];
                $tsB->[1] = $bezB->t_from_xoff($spanx->[1],$offB,[$spanB->[2]->[0],$spanB->[2]->[-1]],$spanB->[0]->[1],$spanB->[3]);

            }
            else {
                $spanx->[1] = $spanB->[1]->[-1];
                $tsB->[1] = $spanB->[2]->[-1];
                $tsA->[1] = $bezA->t_from_xoff($spanx->[1],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);

            }

            my $ysA = [$bezA->Y_offset($tsA->[0],$offA,$spanA->[0]->[1],$spanA->[3]),$bezA->Y_offset($tsA->[1],$offA,$spanA->[0]->[1],$spanA->[3])];
            my $ysB = [$bezB->Y_offset($tsB->[0],$offB,$spanB->[0]->[1],$spanB->[3]),$bezB->Y_offset($tsB->[1],$offB,$spanB->[0]->[1],$spanB->[3])];

            # warn "spanx:\n $spanx->[0],$spanx->[1]\n";
            # warn "ts:\n $tsA->[0],$tsA->[1]\n $tsB->[0],$tsB->[1]\n";
            # warn "ys:\n $ysA->[0],$ysA->[1]\n $ysB->[0],$ysB->[1]\n";

            # one intersection case, similar to how you'd test for crossing line segments
            if    (($ysA->[0] > $ysB->[0] && $ysA->[0] ne $ysB->[0] &&
                    $ysA->[1] < $ysB->[1] && $ysA->[1] ne $ysB->[1]   ) ||
                   ($ysA->[0] < $ysB->[0] && $ysA->[0] ne $ysB->[0] &&
                    $ysA->[1] > $ysB->[1] && $ysA->[1] ne $ysB->[1]   )
                  ) {

                my $findintbezbez = sub {
                    # need Yoff1(t1(xoff)) - Yoff2(t2(xoff))
                    my $t1 = $bezA->t_from_xoff($_[0],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);
                    my $t2 = $bezB->t_from_xoff($_[0],$offB,[$spanB->[2]->[0],$spanB->[2]->[-1]],$spanB->[0]->[1],$spanB->[3]);
                    my $yoff1 = $bezA->Y_offset($t1,$offA,$spanA->[0]->[1],$spanA->[3]);
                    my $yoff2 = $bezB->Y_offset($t2,$offB,$spanB->[0]->[1],$spanB->[3]);
                    return $yoff2 - $yoff1;
                };

                my $bounds_xoff = [$spanx->[0],$spanx->[1]];
                my ($int_xoff,$msg)=BrentsMethod($findintbezbez,$bounds_xoff,0.00001,undef,'subBezOff-subBezOff intersection finding');

                if ($msg) { warn "bezoffint message: $msg\n"; }
                else {
                    my $t1 = $bezA->t_from_xoff($int_xoff,$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);
                    my $intersection_x = $bezA->X_offset($t1,$offA,$spanA->[0]->[1],$spanA->[3]);
                    my $intersection_y = $bezA->Y_offset($t1,$offA,$spanA->[0]->[1],$spanA->[3]);
                    my $x = $bezA->bezierEvalXofT($t1);
                    $int_t_A = $spanA->[0]->[0]->($x);
                    $int_t_B = $spanB->[0]->[0]->($x);
                    #warn "GOT THE ONE INTERSECTION: [$intersection_x,$intersection_y,$int_t_A,$int_t_B]\n";
                    push @ret, [$intersection_x,$intersection_y,$int_t_A,$int_t_B];
                }

            }
            # catch any endpoint overlap "intersections"
            elsif (    $ysA->[0] eq $ysB->[0] || $ysA->[1] eq $ysB->[1]
                    || $ysA->[1] eq $ysB->[0] || $ysA->[0] eq $ysB->[1]
                  ) {
                # what is the policy on this?
                # when checking for self intersections, def don't want these
                # but for diff bezs, these could sometimes be legit intersections
                #push @ret, ["xoverlap","yoverlap"] if $bezA != $bezB;
            }
            # the zero or two intersection case
            else {

                warn "failed to find one offset intersection when a special LUT was provided\n" if defined($lutA) || defined($lutB);
                next if defined($lutA) || defined($lutB);

                # zero intersections if y ranges don't overlap
                my ($lowyA,$highyA) = ($ysA->[0]<$ysA->[1]) ? (@$ysA) : (reverse @$ysA);
                my ($lowyB,$highyB) = ($ysB->[0]<$ysB->[1]) ? (@$ysB) : (reverse @$ysB);
                next if ($lowyA > $highyB);
                next if ($lowyB > $highyA);

                # zero or two intersections

                # more secret sauce
                # for non-offset version it was this:
                # g(t1) = [Y1(t1) - Y2(t2(X1(t1))) ]'

                #die "YOU ARE HERE"; # hey, first tries with properly (hopefully) worked out y_diff_prime sub are now working! might have licked this
                # how about for offset version?
                # g(t1) = [offsetY2(t2_of_xoff(xoff)) - offsetY1(t1_of_xoff(xoff))]'
                #
                # where the form of offsetY(t) is
                #
                # Y(t) + y_offset(t) ,
                #
                # Y(t) being a normal Bezier eval
                #
                # y_offset(t) being something like
                #
                # offset / ( sqrt( (dx/dy)**2 + 1 ) )
                #
                # and [Y(t) + y_offset(t)]' being differentiable using the nice facilities
                # found here in MPath :)

                my $y_diff_prime = sub {

                    # worked it out on paper on 4/28/2018
                    # how I wish I had that paper here now 6/10/2018

                    my $tA = $bezA->t_from_xoff($_[0],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);
                    my $tB = $bezB->t_from_xoff($_[0],$offB,[$spanB->[2]->[0],$spanB->[2]->[-1]],$spanB->[0]->[1],$spanB->[3]);

                    my $nA = $bezA->F_prime(undef,$tA);
                    my $nB = $bezB->F_prime(undef,$tB);

                    my $nA_prime = $bezA->F_2prime(undef,$tA);
                    my $nB_prime = $bezB->F_2prime(undef,$tB);

                    my $YPrimeA = $bezA->bezierEvalYPrimeofT($tA);
                    my $YPrimeB = $bezB->bezierEvalYPrimeofT($tB);

                    my $yoffset_primeA = -($offA/2.0) * 1/(sqrt($nA**2 + 1)**3) * 2*$nA * $nA_prime;
                    my $yoffset_primeB = -($offB/2.0) * 1/(sqrt($nB**2 + 1)**3) * 2*$nB * $nB_prime;

                    my $YoffPrimeA = ($YPrimeA + $yoffset_primeA);
                    my $YoffPrimeB = ($YPrimeB + $yoffset_primeB);

                    my $ret = $YoffPrimeB - $YoffPrimeA;

                    return $ret;

                };

                my $at_start = $y_diff_prime->($spanx->[0]);
                my $at_end   = $y_diff_prime->($spanx->[1]);

                #warn "ydiffprime start, end: [$spanx->[0]] $at_start, [$spanx->[1]] $at_end\n";

                if ($at_start > 0 && $at_end < 0 || $at_start < 0 && $at_end > 0) {

                    my $bounds_xoffA = [$spanx->[0],$spanx->[1]];
                    my ($split_xoff_A,$msg)=FalsePosition($y_diff_prime,$bounds_xoffA,0.00001,(($bounds_xoffA->[1]-$bounds_xoffA->[0])/2),'subBez-subBez intersection finding - find pair split parameter');

                    warn "split find fail msg: $msg\n" if $msg;

                    my $span_split_x = $split_xoff_A;

                    #warn "split xoff: $span_split_x\n";

                    $split_t_A = $bezA->t_from_xoff($span_split_x,$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);
                    $split_t_B = $bezB->t_from_xoff($span_split_x,$offB,[$spanB->[2]->[0],$spanB->[2]->[-1]],$spanB->[0]->[1],$spanB->[3]);

                    # set up new LUTs to pass to a re-call of this intersection sub
                    # so re-call run will have proper x bounds to find each of
                    # the intersection pair individually.

                    my $sub_t_span_1_A = [$tsA->[0], $split_t_A];
                    my $sub_t_span_2_A = [$split_t_A, $tsA->[1]];
                    my $sub_t_span_1_B = [$tsB->[0], $split_t_B];
                    my $sub_t_span_2_B = [$split_t_B, $tsB->[1]];

                    #warn "re-call 1\n";
                    #warn "  [[$spanA->[0],[$spanx->[0]  , $span_split_x], $spanA->[3]?$sub_t_span_2_A:$sub_t_span_1_A, $spanA->[3]]]\n";
                    #warn "  [[$spanB->[0],[$spanx->[0]  , $span_split_x], $spanB->[3]?$sub_t_span_2_B:$sub_t_span_1_B, $spanB->[3]]]\n";
                    my @int1 = intersect_CoCo($bezA,$bezB,$offA,$offB,
                                                 [[$spanA->[0],[$spanx->[0]  , $span_split_x], $spanA->[3]?$sub_t_span_2_A:$sub_t_span_1_A, $spanA->[3]]],
                                                 [[$spanB->[0],[$spanx->[0]  , $span_split_x], $spanB->[3]?$sub_t_span_2_B:$sub_t_span_1_B, $spanB->[3]]]
                                                 );
                    #warn "re-call 2\n";
                    #warn "  [[$spanA->[0],[$span_split_x, $spanx->[1]  ], $spanA->[3]?$sub_t_span_1_A:$sub_t_span_2_A, $spanA->[3]]]\n";
                    #warn "  [[$spanB->[0],[$span_split_x, $spanx->[1]  ], $spanB->[3]?$sub_t_span_1_B:$sub_t_span_2_B, $spanB->[3]]]\n";

                    my @int2 = intersect_CoCo($bezA,$bezB,$offA,$offB,
                                                 [[$spanA->[0],[$span_split_x, $spanx->[1]  ], $spanA->[3]?$sub_t_span_1_A:$sub_t_span_2_A, $spanA->[3]]],
                                                 [[$spanB->[0],[$span_split_x, $spanx->[1]  ], $spanB->[3]?$sub_t_span_1_B:$sub_t_span_2_B, $spanB->[3]]]
                                                 );

                    push @ret, @int1, @int2;
                }
                else {
                    #warn "no pair";
                }

            }
        }
    }
    return @ret;
}

sub intersect_CoLo {
    # This is similar to the approach of intersect_CoCo(), simplified for the
    # line part of the math.

    my ($bezA, $lineB, $offA, $offB, $lutA, $lutB) = @_;

    if (!$lutB && !exists $lineB->{XtoTLUT}) {
        # too complex probably for line segments, but hack it in here while adapting this bez-bez algo to work with bez-line
        $lineB->{XtoTLUT} = [
                             [
                              [sub {$lineB->solveXforTheta($_[0])}], # bez would have 3 or 6 sub refs here. line only needs 1 or 2. 1 for our purposes here so far
                              [$lineB->{minx},$lineB->{maxx}],
                              [0,1],
                              ($lineB->{p2}->[0] >= $lineB->{p1}->[0] ? 0 : 1)
                             ]
                            ];
    }

    my $XtoTLUT_A = $lutA ? $lutA : [map [$_->[0],[$_->[1]->[0],$_->[1]->[-1]],[$_->[2]->[0],$_->[2]->[-1]],$_->[3]], @{$bezA->{XtoTLUT}}];
    my $XtoTLUT_B = $lutB ? $lutB : [map [$_->[0],[$_->[1]->[0],$_->[1]->[-1]],[$_->[2]->[0],$_->[2]->[-1]],$_->[3]], @{$lineB->{XtoTLUT}}];

    # assume any provided LUTs already contain offset Xs in their x range entries
    # otherwise we need to do that now for the LUT copies we made
    if (!defined($lutA)) {
        foreach my $span (@{$XtoTLUT_A}) {
            $span->[1]->[0]  = $bezA->X_offset($span->[2]->[0] ,$offA,$span->[0]->[1],$span->[3]);
            $span->[1]->[-1] = $bezA->X_offset($span->[2]->[-1],$offA,$span->[0]->[1],$span->[3]);
        }
    }

    if (!defined($lutB) && (0 + $offB) ne '0') {
        foreach my $span (@{$XtoTLUT_B}) {
            $span->[1]->[0]  = $lineB->X_offset($span->[2]->[0] ,$offB);
            $span->[1]->[-1] = $lineB->X_offset($span->[2]->[-1],$offB);
        }
    }

    my @ret;

    my %aseen;

    my $offintloopcnt=0;

    foreach my $spanA (@{$XtoTLUT_A}) {

        $aseen{$spanA}={} if !$aseen{$spanA};

        #foreach my $spanB (@{$XtoTLUT_B}) { # don't need this loop for line - already deleted closing bracket for this and unindented old body following here
        $spanB = $XtoTLUT_B->[0]; # for the line, we only expect one lut entry, and only because we're hacking in a lut for the line - prob don't even want that - simplify later

        # skip if x spans don't overlap
        next if ($spanA->[1]->[0] > $spanB->[1]->[-1]);
        next if ($spanB->[1]->[0] > $spanA->[1]->[-1]);

        # x clipping to minimum x span that's valid for both lut entry spans

        my $spanx = [];
        my $tsA = [];
        my $tsB = [];

        # the x clipping is more involved with these offset segments

        if ($spanA->[1]->[0] > $spanB->[1]->[0]) {
            $spanx->[0] = $spanA->[1]->[0];
            $tsA->[0] = $spanA->[2]->[0];
            $tsB->[0] = $lineB->t_from_xoff($spanx->[0],$offB);

        }
        else {
            $spanx->[0] = $spanB->[1]->[0];
            $tsB->[0] = $spanB->[2]->[0];
            $tsA->[0] = $bezA->t_from_xoff($spanx->[0],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);

        }
        if ($spanA->[1]->[-1] < $spanB->[1]->[-1]) {
            $spanx->[1] = $spanA->[1]->[-1];
            $tsA->[1] = $spanA->[2]->[-1];
            $tsB->[1] = $lineB->t_from_xoff($spanx->[1],$offB);

        }
        else {
            $spanx->[1] = $spanB->[1]->[-1];
            $tsB->[1] = $spanB->[2]->[-1];
            $tsA->[1] = $bezA->t_from_xoff($spanx->[1],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);

        }

        my $ysA = [$bezA->Y_offset($tsA->[0],$offA,$spanA->[0]->[1],$spanA->[3]),$bezA->Y_offset($tsA->[1],$offA,$spanA->[0]->[1],$spanA->[3])];
        my $ysB = [$lineB->Y_offset($tsB->[0],$offB),$lineB->Y_offset($tsB->[1],$offB)];

        # warn "spanx:\n $spanx->[0],$spanx->[1]\n";
        # warn "ts:\n $tsA->[0],$tsA->[1]\n $tsB->[0],$tsB->[1]\n";
        # warn "ys:\n $ysA->[0],$ysA->[1]\n $ysB->[0],$ysB->[1]\n";

        # one intersection case, similar to how you'd test for crossing line segments
        if    (($ysA->[0] > $ysB->[0] && $ysA->[0] ne $ysB->[0] &&
                $ysA->[1] < $ysB->[1] && $ysA->[1] ne $ysB->[1]   ) ||
               ($ysA->[0] < $ysB->[0] && $ysA->[0] ne $ysB->[0] &&
                $ysA->[1] > $ysB->[1] && $ysA->[1] ne $ysB->[1]   )
              ) {

            my $findintbezline = sub {
                # need Yoff1(t1(xoff)) - Yoff2(t2(xoff))
                my $t1 = $bezA->t_from_xoff($_[0],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);
                my $t2 = $lineB->t_from_xoff($_[0],$offB);
                my $yoff1 = $bezA->Y_offset($t1,$offA,$spanA->[0]->[1],$spanA->[3]);
                my $yoff2 = $lineB->Y_offset($t2,$offB);
                return $yoff2 - $yoff1;
            };

            my $bounds_xoff = [$spanx->[0],$spanx->[1]];
            my ($int_xoff,$msg)=BrentsMethod($findintbezline,$bounds_xoff,0.00001,undef,'subBezOff-lineOff intersection finding');

            if ($msg) { warn "subBezOff-lineOff message: $msg\n"; }
            else {
                my $t1 = $bezA->t_from_xoff($int_xoff,$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);
                my $intersection_x = $bezA->X_offset($t1,$offA,$spanA->[0]->[1],$spanA->[3]);
                my $intersection_y = $bezA->Y_offset($t1,$offA,$spanA->[0]->[1],$spanA->[3]);
                my $x = $bezA->bezierEvalXofT($t1);
                $int_t_A = $spanA->[0]->[0]->($x);
                $int_t_B = $spanB->[0]->[0]->($x);
                #warn "GOT THE ONE INTERSECTION: [$intersection_x,$intersection_y,$int_t_A,$int_t_B]\n";
                push @ret, [$intersection_x,$intersection_y,$int_t_A,$int_t_B];
            }

        }
        # catch any endpoint overlap "intersections"
        elsif (    $ysA->[0] eq $ysB->[0] || $ysA->[1] eq $ysB->[1]
                || $ysA->[1] eq $ysB->[0] || $ysA->[0] eq $ysB->[1]
              ) {
            # what is the policy on this?
            # when checking for self intersections, def don't want these
            # but for diff bezs, these could sometimes be legit intersections
            #push @ret, ["xoverlap","yoverlap"] if $bezA != $bezB;
        }
        # the zero or two intersection case
        else {

            warn "failed to find one offset intersection when a special LUT was provided\n" if defined($lutA) || defined($lutB);
            next if defined($lutA) || defined($lutB);

            # zero intersections if y ranges don't overlap
            my ($lowyA,$highyA) = ($ysA->[0]<$ysA->[1]) ? (@$ysA) : (reverse @$ysA);
            my ($lowyB,$highyB) = ($ysB->[0]<$ysB->[1]) ? (@$ysB) : (reverse @$ysB);
            next if ($lowyA > $highyB);
            next if ($lowyB > $highyA);

            # zero or two intersections

            # See intersect_CoCo() for note on what this is.
            # Here it's pared down for CoLo case.
            # There we were looking for this:
            # g(t1) = offsetY2(t2_of_xoff(xoff))' - offsetY1(t1_of_xoff(xoff))'
            # and had to use xoff as the common thing between the Beziers, to
            # get corresponding t1 and t2.
            # But for a line, offsetY(t)' (same as Y(t)') is always the same -
            # for any t, and for any xoff. And it's dy/dt, which here is the
            # same as the paramaterized line's delta_y / delta_t.
            # And delta_t == 1.
            # So offsetY(t)' = line->{maxy} - line->{miny}
            # Seems right on first tests.

            my $y_diff_prime = sub {

                my $tA = $bezA->t_from_xoff($_[0],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);
                my $nA = $bezA->F_prime(undef,$tA);
                my $nA_prime = $bezA->F_2prime(undef,$tA);
                my $YPrimeA = $bezA->bezierEvalYPrimeofT($tA);
                my $yoffset_primeA = -($offA/2.0) * 1/(sqrt($nA**2 + 1)**3) * 2*$nA * $nA_prime;
                my $YoffPrimeA = ($YPrimeA + $yoffset_primeA);

                my $YoffPrimeB = $lineB->{dy};

                my $ret = $YoffPrimeB - $YoffPrimeA;

                return $ret;

            };

            my $at_start = $y_diff_prime->($spanx->[0]);
            my $at_end   = $y_diff_prime->($spanx->[1]);

            #warn "ydiffprime start, end: [$spanx->[0]] $at_start, [$spanx->[1]] $at_end\n";

            if ($at_start > 0 && $at_end < 0 || $at_start < 0 && $at_end > 0) {

                my $bounds_xoffA = [$spanx->[0],$spanx->[1]];
                my ($split_xoff_A,$msg)=FalsePosition($y_diff_prime,$bounds_xoffA,0.00001,(($bounds_xoffA->[1]-$bounds_xoffA->[0])/2),'subBez-line intersection finding - find pair split parameter');

                warn "split find fail msg: $msg\n" if $msg;

                my $span_split_x = $split_xoff_A;

                #warn "split xoff: $span_split_x\n";

                $split_t_A = $bezA->t_from_xoff($span_split_x,$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);
                $split_t_B = $lineB->t_from_xoff($span_split_x,$offB);

                # set up new LUTs to pass to a re-call of this intersection sub
                # so re-call run will have proper x bounds to find each of
                # the intersection pair individually.

                my $sub_t_span_1_A = [$tsA->[0], $split_t_A];
                my $sub_t_span_2_A = [$split_t_A, $tsA->[1]];
                my $sub_t_span_1_B = [$tsB->[0], $split_t_B];
                my $sub_t_span_2_B = [$split_t_B, $tsB->[1]];

                #warn "re-call 1\n";
                #warn "  [[$spanA->[0],[$spanx->[0]  , $span_split_x], $spanA->[3]?$sub_t_span_2_A:$sub_t_span_1_A, $spanA->[3]]]\n";
                #warn "  [[$spanB->[0],[$spanx->[0]  , $span_split_x], $spanB->[3]?$sub_t_span_2_B:$sub_t_span_1_B, $spanB->[3]]]\n";
                my @int1 = intersect_CoLo($bezA,$lineB,$offA,$offB,
                                             [[$spanA->[0],[$spanx->[0]  , $span_split_x], $spanA->[3]?$sub_t_span_2_A:$sub_t_span_1_A, $spanA->[3]]],
                                             [[$spanB->[0],[$spanx->[0]  , $span_split_x], $spanB->[3]?$sub_t_span_2_B:$sub_t_span_1_B, $spanB->[3]]]
                                             );
                #warn "re-call 2\n";
                #warn "  [[$spanA->[0],[$span_split_x, $spanx->[1]  ], $spanA->[3]?$sub_t_span_1_A:$sub_t_span_2_A, $spanA->[3]]]\n";
                #warn "  [[$spanB->[0],[$span_split_x, $spanx->[1]  ], $spanB->[3]?$sub_t_span_1_B:$sub_t_span_2_B, $spanB->[3]]]\n";

                my @int2 = intersect_CoLo($bezA,$lineB,$offA,$offB,
                                             [[$spanA->[0],[$span_split_x, $spanx->[1]  ], $spanA->[3]?$sub_t_span_1_A:$sub_t_span_2_A, $spanA->[3]]],
                                             [[$spanB->[0],[$span_split_x, $spanx->[1]  ], $spanB->[3]?$sub_t_span_1_B:$sub_t_span_2_B, $spanB->[3]]]
                                             );

                push @ret, @int1, @int2;
            }
            else {
                #warn "no pair";
            }

        }

    }
    return @ret;

}
sub intersect_CoL {
    my ($bezA, $lineB, $offA) = @_;
    return intersect_CoLo($bezA, $lineB, $offA, 0);
}

sub _rotate2d {
    my ($origin,$point,$angle) = @_;
    my $dx=($point->[0]-$origin->[0]);
    my $dy=($point->[1]-$origin->[1]);
    #{a c-b d, a d+b c}
    return [$origin->[0] + ($dx*cos($angle) - $dy*sin($angle)),$origin->[1] + ($dx*sin($angle) + $dy*cos($angle))];
    }

} # end package
1;

