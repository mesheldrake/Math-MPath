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
# AL : A1L, AL, A1H, AH, A1V, AV
# AA : A1A1, A1A, AA
# QL : QL, QH, QV
# QA : QA1, QA
# QQ
# CL : CL, CH, CV
# CA : CA1, CA
# CQ
# CC
#
# (A1 is circular arc - special simpler case of elliptical arc)
# (A is general case elliptical arc)
# (Z (closePath) is equivalent to L for intersections)
#





package Math::MPath::Intersections;
{
use Math::MPath::Function::Root qw(BrentsMethod FalsePosition);
use Math::MPath::CubicFormula;
use Math::MPath::QuadraticFormula;

our $pi = 4 * atan2(1,1);

sub intersect_LL {
    my ($seg1, $seg2) = @_;

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
        my $theta1 = ($m1 eq 'Inf')?$seg1->solveYforTheta($segsegret->[1]):$seg1->solveXforTheta($segsegret->[0]);
        my $theta2 = ($m2 eq 'Inf')?$seg2->solveYforTheta($segsegret->[1]):$seg2->solveXforTheta($segsegret->[0]);
        push @$segsegret, $theta1, $theta2;
        push(@ret,$segsegret);
    }
    return @ret;
}

# This one needs a rewrite. There's a better way to do this.
# Look at your t(x) or t(y) code for EllipticalArcs, where you combine
# the equation of a line with the Ellipse equations and solve for ellipse theta
# which you would then convert to your normalized ellipse seg parameter.
sub intersect_AL {
    my ($arc, $line) = @_;
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
        ($line->{maxy}>$arc->{miny}  || $line->{maxy} eq $arc->{miny})&&
        ($line->{miny}<$arc->{maxy}  || $line->{miny} eq $arc->{maxy})
        ) {
        my $rot_line_p1 = [ ($x1*cos(-$arc->{phi_radians}) - $y1*sin(-$arc->{phi_radians})),
                            ($x1*sin(-$arc->{phi_radians}) + $y1*cos(-$arc->{phi_radians})) ];
        my $rot_line_p2 = [ ($x2*cos(-$arc->{phi_radians}) - $y2*sin(-$arc->{phi_radians})),
                            ($x2*sin(-$arc->{phi_radians}) + $y2*cos(-$arc->{phi_radians})) ];

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
            $intersections[$i] = [ ($intersections[$i]->[0]*cos($arc->{phi_radians}) - $intersections[$i]->[1]*sin($arc->{phi_radians})),
                                   ($intersections[$i]->[0]*sin($arc->{phi_radians}) + $intersections[$i]->[1]*cos($arc->{phi_radians})) ];
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


        foreach my $int (@intersections) {
            # this is a clunky old approach
            my @allArcThetas=$arc->solveXforTheta($int->[0]);
            foreach my $t (@allArcThetas) {
                my $tp=$arc->point($t);
                if (abs($tp->[1] - $int->[1]) < 0.0000000001) {
                    push(@$int,$t);

                }
            }
            die "missed capturing an arc segment parameter in Arc Line intersect" if scalar(@$int) != 3;
            push(@$int,($line->{m} eq 'inf')?$line->solveYforTheta($int->[1]):$line->solveXforTheta($int->[0]));
        }
    

        push(@ret,@intersections);

    }

    return @ret;
}

# "optimized" circle-circle intersection case
# don't know whether it's actually better than the general elliptical arc case
# at this point, but should revisit this math and improve with new features
# of elliptical arc code if possible
sub intersect_A1A1 {
    my ($arc1, $arc2) = @_;

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


        # a reason to revisit this - this is a bad approach for looking up
        # the arc parameters. Whatever your approach is above, it should probably
        # have calculated these already, and more reliably.
        foreach my $int (@intersections) {
            my @allArcThetas=$arc1->solveXforTheta($int->[0]);
            foreach my $t (@allArcThetas) {
                my $tp=$arc1->point($t);
                if (abs($tp->[1] - $int->[1]) < 0.0000000001) {push(@$int,$t);}
            }
            die "missed capturing an arc segment parameter in Arc Line intersect (1)" if scalar(@$int) != 3;
            @allArcThetas=$arc2->solveXforTheta($int->[0]);
            foreach my $t (@allArcThetas) {
                my $tp=$arc2->point($t);
                if (abs($tp->[1] - $int->[1]) < 0.0000000001) {push(@$int,$t);}
            }
            die "missed capturing an arc segment parameter in Arc Line intersect (2)" if scalar(@$int) != 4;
        }

        push(@ret,@intersections);

        }

    return @ret;

}

sub intersect_AA {

    # just intersect_CC, adapted for arcs

    my ($arcA, $arcB, $lutA, $lutB) = @_;

    my $XtoTLUT_A = $lutA ? $lutA : $arcA->{XtoTLUT};
    my $XtoTLUT_B = $lutB ? $lutB : $arcB->{XtoTLUT};

    my @ret;

    my %aseen;

    foreach my $spanA (@{$XtoTLUT_A}) {

        #$aseen{$spanA}={} if !$aseen{$spanA};

        foreach my $spanB (@{$XtoTLUT_B}) {

        # avoid duplicate runs when looking for self intersection
        #next if $spanA == $spanB;
        #next if $aseen{$spanA}->{$spanB};
        #$aseen{$spanB}->{$spanA}++;

        next if ($spanA->[1]->[0] > $spanB->[1]->[-1]);
        next if ($spanB->[1]->[0] > $spanA->[1]->[-1]);

        my $spanx=[];
        my $tsA=[];
        my $tsB=[];

        if ($spanA->[1]->[0] eq $spanB->[1]->[0]) {
            $spanx->[0] = $spanA->[1]->[0];
            $tsA->[0] = $spanA->[2]->[0];
            $tsB->[0] = $spanB->[2]->[0];
        }
        elsif ($spanA->[1]->[0] > $spanB->[1]->[0]) {
            $spanx->[0] = $spanA->[1]->[0];
            $tsA->[0] = $spanA->[2]->[0];
            #$tsB->[0] = $arcB->t_from_xoff($spanx->[0],$offB,[$spanB->[2]->[0],$spanB->[2]->[-1]],$spanB->[0]->[1]);
            $tsB->[0] = $spanB->[0]->[0]->($spanx->[0]);
        }
        else {
            $spanx->[0] = $spanB->[1]->[0];
            $tsB->[0] = $spanB->[2]->[0];
            #$tsA->[0] = $arcA->t_from_xoff($spanx->[0],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1]);
            $tsA->[0] = $spanA->[0]->[0]->($spanx->[0]);
        }
        if ($spanA->[1]->[-1] eq $spanB->[1]->[-1]) {
            $spanx->[1] = $spanA->[1]->[-1];
            $tsA->[1] = $spanA->[2]->[-1];
            $tsB->[1] = $spanB->[2]->[-1];
        }
        elsif ($spanA->[1]->[-1] < $spanB->[1]->[-1]) {
            $spanx->[1] = $spanA->[1]->[-1];
            $tsA->[1] = $spanA->[2]->[-1];
            #$tsB->[1] = $arcB->t_from_xoff($spanx->[1],$offB,[$spanB->[2]->[0],$spanB->[2]->[-1]],$spanB->[0]->[1]);
            $tsB->[1] = $spanB->[0]->[0]->($spanx->[1]);
        }
        else {
            $spanx->[1] = $spanB->[1]->[-1];
            $tsB->[1] = $spanB->[2]->[-1];
            #$tsA->[1] = $arcA->t_from_xoff($spanx->[1],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1]);
            $tsA->[1] = $spanA->[0]->[0]->($spanx->[1]);
        }

        my $ysA = [$arcA->evalYofTheta($tsA->[0]),$arcA->evalYofTheta($tsA->[1])];
        my $ysB = [$arcB->evalYofTheta($tsB->[0]),$arcB->evalYofTheta($tsB->[1])];

        #warn "xspan:\n $spanx->[0],$spanx->[1]\n";
        #warn "ts:\n $tsA->[0],$tsA->[1]\n $tsB->[0],$tsB->[1]\n";
        #warn "ys:\n $ysA->[0],$ysA->[1]\n $ysB->[0],$ysB->[1]\n";

        # one intersection case, similar to how you'd test for crossing line segments
        if    (($ysA->[0] > $ysB->[0] && $ysA->[0] ne $ysB->[0] &&
                $ysA->[1] < $ysB->[1] && $ysA->[1] ne $ysB->[1]   ) ||
               ($ysA->[0] < $ysB->[0] && $ysA->[0] ne $ysB->[0] &&
                $ysA->[1] > $ysB->[1] && $ysA->[1] ne $ysB->[1]   )
              ) {

            my $findintarcarc = sub {
                # Y1(t1) - Y2(t2(X1(t1))) secret sauce ingredient
                return $arcA->evalYofTheta($_[0]) - $arcB->evalYofTheta( $spanB->[0]->[0]->( $arcA->evalXofTheta($_[0]) ) );
            };
            my $bounds_tA = ($tsA->[0] < $tsA->[1]) ? [$tsA->[0],$tsA->[1]] : [$tsA->[1],$tsA->[0]];
            my ($int_t_A,$msg)=BrentsMethod($findintarcarc,$bounds_tA,0.00001,undef,'subArc-subArc intersection finding');

            my $intersection_x = $arcA->evalXofTheta($int_t_A);
            my $intersection_y = $arcA->evalYofTheta($int_t_A);

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
            #push @ret, ["xoverlap","yoverlap"] if $arcA != $arcB;
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
                my $ret =
                $arcA->evalYPrimeofTheta($_[0])
                -
                $arcB->evalYPrimeofTheta( $spanB->[0]->[0]->( $arcA->evalXofTheta($_[0]) ) )
                *
                $spanB->[0]->[1]->($arcA->evalXofTheta($_[0]))
                *
                $arcA->evalXPrimeofTheta($_[0]);
                return $ret;
            };

            my $at_start = $y_diff_prime->($tsA->[0]);
            my $at_end   = $y_diff_prime->($tsA->[1]);

            #warn "ydiffprime start, end: $at_start, $at_end\n";

            if ($at_start > 0 && $at_end < 0 || $at_start < 0 && $at_end > 0) {

                my $bounds_tA = ($tsA->[0] < $tsA->[1]) ? [$tsA->[0],$tsA->[1]] : [$tsA->[1],$tsA->[0]];
                my ($split_t_A,$msg)=BrentsMethod($y_diff_prime,$bounds_tA,0.00001,undef,'subArc-subArc intersection finding - find pair split parameter');

                my $span_split_x = $arcA->evalXofTheta($split_t_A);
                
                #warn "split: $span_split_x\n";

                $split_t_B = $spanB->[0]->[0]->($span_split_x);

                # set up new LUTs to pass to a re-call of this intersection sub
                # so re-call run will have proper x bounds to find each of
                # the intersection pair individually.

                my $sub_t_span_1_A = [$tsA->[0], $split_t_A];
                my $sub_t_span_2_A = [$split_t_A, $tsA->[1]];
                my $sub_t_span_1_B = [$tsB->[0], $split_t_B];
                my $sub_t_span_2_B = [$split_t_B, $tsB->[1]];

                #warn "[[$spanA->[0],[$spanx->[0]  , $span_split_x], $spanA->[3]?@$sub_t_span_2_A:@$sub_t_span_1_A, $spanA->[3]]]\n";
                #warn "[[$spanB->[0],[$spanx->[0]  , $span_split_x], $spanB->[3]?@$sub_t_span_2_B:@$sub_t_span_1_B, $spanB->[3]]]\n";
                my @int1 = intersect_AA($arcA,$arcB,
                                             [[$spanA->[0],[$spanx->[0]  , $span_split_x], $spanA->[3]?$sub_t_span_2_A:$sub_t_span_1_A, $spanA->[3]]],
                                             [[$spanB->[0],[$spanx->[0]  , $span_split_x], $spanB->[3]?$sub_t_span_2_B:$sub_t_span_1_B, $spanB->[3]]]
                                             );

                #warn "[[$spanA->[0],[$span_split_x, $spanx->[1]  ], $spanA->[3]?@$sub_t_span_1_A:@$sub_t_span_2_A, $spanA->[3]]]\n";
                #warn "[[$spanB->[0],[$span_split_x, $spanx->[1]  ], $spanB->[3]?@$sub_t_span_1_B:@$sub_t_span_2_B, $spanB->[3]]]\n";
                my @int2 = intersect_AA($arcA,$arcB,
                                             [[$spanA->[0],[$span_split_x, $spanx->[1]  ], $spanA->[3]?$sub_t_span_1_A:$sub_t_span_2_A, $spanA->[3]]],
                                             [[$spanB->[0],[$span_split_x, $spanx->[1]  ], $spanB->[3]?$sub_t_span_1_B:$sub_t_span_2_B, $spanB->[3]]]
                                             );

                warn "expecting two intersections but got [",(scalar(@int1)," + ",scalar(@int2)),"]" if (scalar(@int1) != 1 || scalar(@int2) != 1);

                push @ret, $int1[0]->[2] <= $int2[0]->[2]
                           ? ($int1[0],$int2[0])
                           : ($int2[0],$int1[0]);
            }
            else {
                #warn "no pair";
            }

        }

        }
    }

    return @ret;
}

sub intersect_CL {
    my ($curve, $line) = @_;
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
        @ts = $curve->solveXforTheta($line->{maxx});
        foreach my $t (@ts) {
            my $y = $curve->bezierEvalYofT($t);
            if (($y < $line->{maxy} || $y eq $line->{maxy}) &&
                ($y > $line->{miny} || $y eq $line->{miny})) {
                push(@ret,[$line->{maxx},$y,$t,$line->solveYforTheta($y)]);
            }
        }
    }
    elsif ($line->{m} eq 0) {
        @ts = $curve->solveYforTheta($line->{p1}->[1]);
        foreach my $t (@ts) {
            my $x = $curve->bezierEvalXofT($t);
            if (($x < $line->{maxx} || $x eq $line->{maxx}) &&
                ($x > $line->{minx} || $x eq $line->{minx})) {
                push(@ret,[$x,$line->{p1}->[1],$t,$line->solveXforTheta($x)]);
            }
        }
    }
    else {

        @ts = &cubicformula(
            ($curve->{F}-$line->{m}*$curve->{B})            / ($curve->{E}-$line->{m}*$curve->{A}),
            ($curve->{G}-$line->{m}*$curve->{C})            / ($curve->{E}-$line->{m}*$curve->{A}),
            ($curve->{H}-$line->{m}*$curve->{D}-$line->{b}) / ($curve->{E}-$line->{m}*$curve->{A}),
            1);

        @ts = sort {$a<=>$b} grep {(1 > $_ || 1 eq $_) && ($_ > 0 || $_ eq 0)}  @ts;

        foreach my $t (@ts) {
            my $x = $curve->bezierEvalXofT($t);
            if (($x < $line->{maxx} || $x eq $line->{maxx}) &&
                ($x > $line->{minx} || $x eq $line->{minx})) {
                push(@ret,[$x,$curve->bezierEvalYofT($t),$t,$line->solveXforTheta($x)]);
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

                warn "expecting two intersections but got [",(scalar(@int1)," + ",scalar(@int2)),"]" if (scalar(@int1) != 1 || scalar(@int2) != 1);

                push @ret, $int1[0]->[2] <= $int2[0]->[2]
                           ? ($int1[0],$int2[0])
                           : ($int2[0],$int1[0]);
            }
            else {
                #warn "no pair";
            }

        }

        }
    }

    return @ret;
}

sub intersect_CA {

    # like intersect_CC, with second bezier stuff replaced with arc stuff

    my ($bezA, $arcB, $lutA, $lutB) = @_;

    my $XtoTLUT_A = $lutA ? $lutA : $bezA->{XtoTLUT};
    my $XtoTLUT_B = $lutB ? $lutB : $arcB->{XtoTLUT};

    my @ret;

    #my %aseen;

    foreach my $spanA (@{$XtoTLUT_A}) {

        #$aseen{$spanA}={} if !$aseen{$spanA};

        foreach my $spanB (@{$XtoTLUT_B}) {

        # avoid duplicate runs when looking for self intersection
        #next if $spanA == $spanB;
        #next if $aseen{$spanA}->{$spanB};
        #$aseen{$spanB}->{$spanA}++;

        next if ($spanA->[1]->[0] > $spanB->[1]->[-1]);
        next if ($spanB->[1]->[0] > $spanA->[1]->[-1]);

        my $spanx = [$spanA->[1]->[0] > $spanB->[1]->[0] ? $spanA->[1]->[0] : $spanB->[1]->[0], $spanA->[1]->[-1] < $spanB->[1]->[-1] ? $spanA->[1]->[-1] : $spanB->[1]->[-1]];

        my $tsA = [$spanA->[0]->[0]->($spanx->[0]), $spanA->[0]->[0]->($spanx->[1])];
        my $tsB = [$spanB->[0]->[0]->($spanx->[0]), $spanB->[0]->[0]->($spanx->[1])];

        my $ysA = [$bezA->bezierEvalYofT($tsA->[0]),$bezA->bezierEvalYofT($tsA->[1])];
        my $ysB = [$arcB->evalYofTheta($tsB->[0]),$arcB->evalYofTheta($tsB->[1])];

        #warn "xspan:\n $spanx->[0],$spanx->[1]\n";
        #warn "ts:\n $tsA->[0],$tsA->[1]\n $tsB->[0],$tsB->[1]\n";
        #warn "ys:\n $ysA->[0],$ysA->[1]\n $ysB->[0],$ysB->[1]\n";

        # one intersection case, similar to how you'd test for crossing line segments
        if    (($ysA->[0] > $ysB->[0] && $ysA->[0] ne $ysB->[0] &&
                $ysA->[1] < $ysB->[1] && $ysA->[1] ne $ysB->[1]   ) ||
               ($ysA->[0] < $ysB->[0] && $ysA->[0] ne $ysB->[0] &&
                $ysA->[1] > $ysB->[1] && $ysA->[1] ne $ysB->[1]   )
              ) {

            my $findintbezarc = sub {
                # Y1(t1) - Y2(t2(X1(t1))) secret sauce ingredient
                return $bezA->bezierEvalYofT($_[0]) - $arcB->evalYofTheta( $spanB->[0]->[0]->( $bezA->bezierEvalXofT($_[0]) ) );
            };
            my $bounds_tA = ($tsA->[0] < $tsA->[1]) ? [$tsA->[0],$tsA->[1]] : [$tsA->[1],$tsA->[0]];
            my ($int_t_A,$msg)=BrentsMethod($findintbezarc,$bounds_tA,0.00001,undef,'subBez-subArc intersection finding');

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
            #push @ret, ["xoverlap","yoverlap"] if $bezA != $arcB;
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
                my $bez_x_of_t = $bezA->bezierEvalXofT($_[0]);
                return
                $bezA->bezierEvalYPrimeofT($_[0])
                -
                $arcB->evalYPrimeofTheta( $spanB->[0]->[0]->( $bez_x_of_t ) )
                *
                $spanB->[0]->[1]->( $bez_x_of_t )
                *
                $bezA->bezierEvalXPrimeofT( $_[0] );
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

                my @int1 = intersect_CA($bezA,$arcB,
                                             [[$spanA->[0],[$spanx->[0]  , $span_split_x], $spanA->[3]?$sub_t_span_2_A:$sub_t_span_1_A, $spanA->[3]]],
                                             [[$spanB->[0],[$spanx->[0]  , $span_split_x], $spanB->[3]?$sub_t_span_2_B:$sub_t_span_1_B, $spanB->[3]]]
                                             );

                my @int2 = intersect_CA($bezA,$arcB,
                                             [[$spanA->[0],[$span_split_x, $spanx->[1]  ], $spanA->[3]?$sub_t_span_1_A:$sub_t_span_2_A, $spanA->[3]]],
                                             [[$spanB->[0],[$span_split_x, $spanx->[1]  ], $spanB->[3]?$sub_t_span_1_B:$sub_t_span_2_B, $spanB->[3]]]
                                             );

                warn "expecting two intersections but got [",(scalar(@int1)," + ",scalar(@int2)),"]" if (scalar(@int1) != 1 || scalar(@int2) != 1);

                push @ret, $int1[0]->[2] <= $int2[0]->[2]
                           ? ($int1[0],$int2[0])
                           : ($int2[0],$int1[0]);
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
    my ($lineA, $lineB, $offA, $offB) = @_;

    $offA //= (exists($lineA->{offset})?$lineA->{offset}:0);
    $offB //= (exists($lineB->{offset})?$lineB->{offset}:0);
    return intersect_LL($lineA, $lineB) if (!$offA && !$offB);

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

    return intersect_LL($lineA, $lineB);
}

sub intersect_A1oLo {
    my ($arc, $line, $offCircle, $offLine) = @_;

    $offCircle //= (exists($arc->{offset})?$arc->{offset}:0);
    $offLine //= (exists($line->{offset})?$line->{offset}:0);
    return intersect_AL($arc, $line) if (!$offCircle && !$offLine);

    # will the parameters at intersections for these temporary offset segs
    # really correspond to the parameter you would use with the original
    # seg, with seg->point(t), to get that same intersection point?
    # Sounds like you need a test for that.

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

    return intersect_AL($arc_off, $line_off);
}

sub intersect_A1oA1o {
    my ($arcA, $arcB, $offA, $offB) = @_;

    $offA //= (exists($arcA->{offset})?$arcA->{offset}:0);
    $offB //= (exists($arcB->{offset})?$arcB->{offset}:0);
    return intersect_A1A1($arcA, $arcB) if (!$offA && !$offB);

    # see parameter correspondence question in intersect_A1oLo()

    my $arcA_off_radius = $arcA->{rx} + $offA;
    my $arcA_off = Math::MPath::EllipticalArc->new(
        [$arcA->{p1}->[0], $arcA->{p1}->[1]],
        [$arcA_off_radius,$arcA_off_radius],
        $arcA->{phi},
        $arcA->{large_arc_flag},
        $arcA->{sweep_flag},
        [$arcA->{p2}->[0], $arcA->{p2}->[1]],
        $arcA->{precision},
        $arcA->{isLite}
    );
    my $arcB_off_radius = $arcB->{rx} + $offB;
    my $arcB_off = Math::MPath::EllipticalArc->new(
        [$arcB->{p1}->[0], $arcB->{p1}->[1]],
        [$arcB_off_radius,$arcB_off_radius],
        $arcB->{phi},
        $arcB->{large_arc_flag},
        $arcB->{sweep_flag},
        [$arcB->{p2}->[0], $arcB->{p2}->[1]],
        $arcB->{precision},
        $arcB->{isLite}
    );

    return intersect_A1A1($arcA, $arcB, $wantThetas);
}

sub intersect_AoLo {
    # This is similar to the approach of intersect_CoLo(),
    # just adapting all the Bezier stuff to be arc stuff.

    my ($arcA, $lineB, $offA, $offB, $lutA, $lutB) = @_;

    $offA //= (exists($arcA->{offset})?$arcA->{offset}:0);
    $offB //= (exists($lineB->{offset})?$lineB->{offset}:0);
    return intersect_AL($arcA, $lineB, $lutA, $lutB) if (!$offA && !$offB);

    if (!$lutB && !exists $lineB->{XtoTLUT}) {
        # too complex probably for line segments, but hack it in here while adapting this bez-bez algo to work with bez-line
        my $line_is_reversed = $lineB->{p2}->[0] >= $lineB->{p1}->[0] ? 0 : 1;
        $lineB->{XtoTLUT} = [
                             [
                              [sub {$lineB->solveXforTheta($_[0])}], # bez would have 3 or 6 sub refs here. line only needs 1 or 2. 1 for our purposes here so far
                              [$lineB->{minx},$lineB->{maxx}],
                              [$line_is_reversed?(1,0):(0,1)],
                              $line_is_reversed
                             ]
                            ];
    }

    my $XtoTLUT_A = $lutA ? $lutA : [map [$_->[0],[$_->[1]->[0],$_->[1]->[-1]],[$_->[2]->[0],$_->[2]->[-1]],$_->[3]], @{$arcA->{XtoTLUT}}];
    my $XtoTLUT_B = $lutB ? $lutB : [map [$_->[0],[$_->[1]->[0],$_->[1]->[-1]],[$_->[2]->[0],$_->[2]->[-1]],$_->[3]], @{$lineB->{XtoTLUT}}];

    # assume any provided LUTs already contain offset Xs in their x range entries
    # otherwise we need to do that now for the LUT copies we made
    if (!defined($lutA)) {
        foreach my $span (@{$XtoTLUT_A}) {
            $span->[1]->[0]  = $arcA->X_offset($span->[2]->[0] ,$offA,$span->[0]->[1],$span->[3]);
            $span->[1]->[-1] = $arcA->X_offset($span->[2]->[-1],$offA,$span->[0]->[1],$span->[3]);
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
            $tsA->[0] = $arcA->t_from_xoff($spanx->[0],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);

        }
        if ($spanA->[1]->[-1] < $spanB->[1]->[-1]) {
            $spanx->[1] = $spanA->[1]->[-1];
            $tsA->[1] = $spanA->[2]->[-1];
            $tsB->[1] = $lineB->t_from_xoff($spanx->[1],$offB);

        }
        else {
            $spanx->[1] = $spanB->[1]->[-1];
            $tsB->[1] = $spanB->[2]->[-1];
            $tsA->[1] = $arcA->t_from_xoff($spanx->[1],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);

        }

        my $ysA = [$arcA->Y_offset($tsA->[0],$offA,$spanA->[0]->[1],$spanA->[3]),$arcA->Y_offset($tsA->[1],$offA,$spanA->[0]->[1],$spanA->[3])];
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

            my $findintarcline = sub {
                # need Yoff1(t1(xoff)) - Yoff2(t2(xoff))
                my $t1 = $arcA->t_from_xoff($_[0],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);
                my $t2 = $lineB->t_from_xoff($_[0],$offB);
                my $yoff1 = $arcA->Y_offset($t1,$offA,$spanA->[0]->[1],$spanA->[3]);
                #warn "### E($t1) $offA : ",($yoff1 - $arcA->evalYofTheta($t1)),"\n";
                my $yoff2 = $lineB->Y_offset($t2,$offB);
                #warn "### L($t2) $offB : ",($yoff2 - $lineB->point($t2)->[1]),"\n";
                return $yoff1 - $yoff2;
            };

            my $bounds_xoff = [$spanx->[0],$spanx->[1]];
            my ($int_xoff,$msg)=BrentsMethod($findintarcline,$bounds_xoff,0.00001,undef,'subArcOff-lineOff intersection finding');

            if ($msg) { warn "subArcOff-lineOff message: $msg\n"; }
            else {
                my $t1 = $arcA->t_from_xoff($int_xoff,$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);
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

            my $y_diff_prime = sub {

                # for the arc
                my $tA = $arcA->t_from_xoff($_[0],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1]);
                my $xA = $arcA->evalXofTheta($tA);
                my $tA_prime  = $spanA->[0]->[1]->($xA);
                my $tA_2prime = $spanA->[0]->[2]->($xA);
                # arc->f_prime()
                my $mA = $arcA->evalYPrimeofTheta($tA) * $tA_prime;
                # arc->f_2prime()
                my $mA_prime = $arcA->evalYDoublePrimeofTheta($tA) * $tA_prime**2 + $arcA->evalYPrimeofTheta($tB) * $tA_2prime;
                my $YPrimeA = $arcA->evalYPrimeofTheta($tA);
                my $yoffset_primeA = -($offA/2.0) * 1/(sqrt($mA**2 + 1)**3) * 2*$mA * $mA_prime;
                $YPrimeA *= -1 if $spanA->[3];
                $yoffset_primeA *= -1 if $spanA->[3];
                my $YoffPrimeA = ($YPrimeA + $yoffset_primeA);

                # for the Line
                my $YoffPrimeB = $lineB->{dy}; # dy/dt which is (y2-y1)/(1-0)
                $YoffPrimeB *= -1 if $spanB->[3];

                my $ret = $YoffPrimeA - $YoffPrimeB;

                return $ret;
            };

            my $at_start = $y_diff_prime->($spanx->[0]);
            my $at_end   = $y_diff_prime->($spanx->[1]);

            #warn "ydiffprime start, end: [$spanx->[0]] $at_start, [$spanx->[1]] $at_end\n";

            if ($at_start > 0 && $at_end < 0 || $at_start < 0 && $at_end > 0) {

                my $bounds_xoffA = [$spanx->[0],$spanx->[1]];
                my ($split_xoff_A,$msg)=FalsePosition($y_diff_prime,$bounds_xoffA,0.00001,(($bounds_xoffA->[1]-$bounds_xoffA->[0])/2),'subArc-line intersection finding - find pair split parameter');

                warn "split find fail msg: $msg\n" if $msg;

                my $span_split_x = $split_xoff_A;

                #warn "split xoff: $span_split_x\n";

                $split_t_A = $arcA->t_from_xoff($span_split_x,$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);
                $split_t_B = $lineB->t_from_xoff($span_split_x,$offB);

                # set up new LUTs to pass to a re-call of this intersection sub
                # so re-call run will have proper x bounds to find each of
                # the intersection pair individually.

                my $sub_t_span_1_A = [$tsA->[$spanA->[3]?1:0], $split_t_A];
                my $sub_t_span_2_A = [$split_t_A, $tsA->[$spanA->[3]?0:1]];
                my $sub_t_span_1_B = [$tsB->[$spanB->[3]?1:0], $split_t_B];
                my $sub_t_span_2_B = [$split_t_B, $tsB->[$spanB->[3]?0:1]];

                #warn "re-call 1\n";
                #warn "  [[$spanA->[0],[$spanx->[0]  , $span_split_x], $spanA->[3]?$sub_t_span_2_A:$sub_t_span_1_A, $spanA->[3]]]\n";
                #warn "  [[$spanB->[0],[$spanx->[0]  , $span_split_x], $spanB->[3]?$sub_t_span_2_B:$sub_t_span_1_B, $spanB->[3]]]\n";
                my @int1 = intersect_AoLo($arcA,$lineB,$offA,$offB,
                                             [[$spanA->[0],[$spanx->[0]  , $span_split_x], $spanA->[3]?$sub_t_span_2_A:$sub_t_span_1_A, $spanA->[3]]],
                                             [[$spanB->[0],[$spanx->[0]  , $span_split_x], $spanB->[3]?$sub_t_span_2_B:$sub_t_span_1_B, $spanB->[3]]]
                                             );

                #warn "re-call 2\n";
                #warn "  [[$spanA->[0],[$span_split_x, $spanx->[1]  ], $spanA->[3]?@$sub_t_span_1_A:@$sub_t_span_2_A, $spanA->[3]]]\n";
                #warn "  [[$spanB->[0],[$span_split_x, $spanx->[1]  ], $spanB->[3]?@$sub_t_span_1_B:@$sub_t_span_2_B, $spanB->[3]]]\n";

                my @int2 = intersect_AoLo($arcA,$lineB,$offA,$offB,
                                             [[$spanA->[0],[$span_split_x, $spanx->[1]  ], $spanA->[3]?$sub_t_span_1_A:$sub_t_span_2_A, $spanA->[3]]],
                                             [[$spanB->[0],[$span_split_x, $spanx->[1]  ], $spanB->[3]?$sub_t_span_1_B:$sub_t_span_2_B, $spanB->[3]]]
                                             );

                warn "expecting two intersections but got [",(scalar(@int1)," + ",scalar(@int2)),"]offs:[$offA,$offB]" if (scalar(@int1) != 1 || scalar(@int2) != 1);

                push @ret, $int1[0]->[2] <= $int2[0]->[2]
                           ? ($int1[0],$int2[0])
                           : ($int2[0],$int1[0]);
            }
            else {
                #warn "no pair";
            }

        }

    }
    return @ret;

}

sub intersect_AoAo {

    # Adaptation of CoCo (offset Bezier curve intersections), 
    # for offset elliptical arcs.
    
    # lutA and lutB optional - used to re-call this function in the
    # two intersection case

    my ($arcA, $arcB, $offA, $offB, $lutA, $lutB) = @_;

    $offA //= (exists($arcA->{offset})?$arcA->{offset}:0);
    $offB //= (exists($arcB->{offset})?$arcB->{offset}:0);
    return intersect_AA($arcA, $arcB, $lutA, $lutB) if (!$offA && !$offB);

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

                    # for the first arc
                    my $tA = $arcA->t_from_xoff($_[0],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1]);
                    my $xA = $arcA->evalXofTheta($tA);
                    my $tA_prime  = $spanA->[0]->[1]->($xA);
                    my $tA_2prime = $spanA->[0]->[2]->($xA);
                    # arc->f_prime()
                    my $mA = $arcA->evalYPrimeofTheta($tA) * $tA_prime;
                    # arc->f_2prime()
                    my $mA_prime = $arcA->evalYDoublePrimeofTheta($tA) * $tA_prime**2 + $arcA->evalYPrimeofTheta($tB) * $tA_2prime;
                    my $YPrimeA = $arcA->evalYPrimeofTheta($tA);
                    my $yoffset_primeA = -($offA/2.0) * 1/(sqrt($mA**2 + 1)**3) * 2*$mA * $mA_prime;
                    $YPrimeA *= -1 if $spanA->[3];
                    $yoffset_primeA *= -1 if $spanA->[3];
                    my $YoffPrimeA = ($YPrimeA + $yoffset_primeA);

                    # for the second arc
                    my $tB = $arcB->t_from_xoff($_[0],$offB,[$spanB->[2]->[0],$spanB->[2]->[-1]],$spanB->[0]->[1]);
                    my $xB = $arcB->evalXofTheta($tB);
                    my $tB_prime  = $spanB->[0]->[1]->($xB);
                    my $tB_2prime = $spanB->[0]->[2]->($xB);
                    # arc->f_prime()
                    my $mB = $arcB->evalYPrimeofTheta($tB) * $tB_prime;
                    # arc->f_2prime()
                    my $mB_prime = $arcB->evalYDoublePrimeofTheta($tB) * $tB_prime**2 + $arcB->evalYPrimeofTheta($tB) * $tB_2prime;
                    my $YPrimeB = $arcB->evalYPrimeofTheta($tB);
                    my $yoffset_primeB = -($offB/2.0) * 1/(sqrt($mB**2 + 1)**3) * 2*$mB * $mB_prime;
                    $YPrimeB *= -1 if $spanB->[3];
                    $yoffset_primeB *= -1 if $spanB->[3];
                    my $YoffPrimeB = ($YPrimeB + $yoffset_primeB);

                    my $ret = $YoffPrimeA - $YoffPrimeB;

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

                    my $sub_t_span_1_A = [$tsA->[$spanA->[3]?1:0], $split_t_A];
                    my $sub_t_span_2_A = [$split_t_A, $tsA->[$spanA->[3]?0:1]];
                    my $sub_t_span_1_B = [$tsB->[$spanB->[3]?1:0], $split_t_B];
                    my $sub_t_span_2_B = [$split_t_B, $tsB->[$spanB->[3]?0:1]];

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

                warn "expecting two intersections but got [",(scalar(@int1)," + ",scalar(@int2)),"]" if (scalar(@int1) != 1 || scalar(@int2) != 1);

                push @ret, $int1[0]->[2] <= $int2[0]->[2]
                           ? ($int1[0],$int2[0])
                           : ($int2[0],$int1[0]);
                }
                else {
                    #warn "no pair";
                }

            }
        }
    }
    return @ret;
}

sub intersect_CoAo {

    # Adaptation of CoCo (offset Bezier curve intersections), 
    # for offset elliptical arcs.
    
    # lutA and lutB optional - used to re-call this function in the
    # two intersection case

    my ($bezA, $arcB, $offA, $offB, $lutA, $lutB) = @_;

    $offA //= (exists($bezA->{offset})?$bezA->{offset}:0);
    $offB //= (exists($arcB->{offset})?$arcB->{offset}:0);
    return intersect_CA($bezA, $arcB, $lutA, $lutB) if (!$offA && !$offB);

    my $XtoTLUT_A = $lutA ? $lutA : [map [$_->[0],[$_->[1]->[0],$_->[1]->[-1]],[$_->[2]->[0],$_->[2]->[-1]],$_->[3]], @{$bezA->{XtoTLUT}}];
    my $XtoTLUT_B = $lutB ? $lutB : [map [$_->[0],[$_->[1]->[0],$_->[1]->[-1]],[$_->[2]->[0],$_->[2]->[-1]],$_->[3]], @{$arcB->{XtoTLUT}}];

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
                $tsA->[0] = $bezA->t_from_xoff($spanx->[0],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1]);

            }
            if ($spanA->[1]->[-1] < $spanB->[1]->[-1]) {
                $spanx->[1] = $spanA->[1]->[-1];
                $tsA->[1] = $spanA->[2]->[-1];
                $tsB->[1] = $arcB->t_from_xoff($spanx->[1],$offB,[$spanB->[2]->[0],$spanB->[2]->[-1]],$spanB->[0]->[1]);

            }
            else {
                $spanx->[1] = $spanB->[1]->[-1];
                $tsB->[1] = $spanB->[2]->[-1];
                $tsA->[1] = $bezA->t_from_xoff($spanx->[1],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1]);

            }

            my $ysA = [$bezA->Y_offset($tsA->[0],$offA,$spanA->[0]->[1],$spanA->[3]),$bezA->Y_offset($tsA->[1],$offA,$spanA->[0]->[1],$spanA->[3])];
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

                my $findintbezarc = sub {
                    # need Yoff1(t1(xoff)) - Yoff2(t2(xoff))
                    my $t1 = $bezA->t_from_xoff($_[0],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1]);
                    my $t2 = $arcB->t_from_xoff($_[0],$offB,[$spanB->[2]->[0],$spanB->[2]->[-1]],$spanB->[0]->[1]);

                    my $yoff1 = $bezA->Y_offset($t1,$offA,$spanA->[0]->[1],$spanA->[3]);
                    my $yoff2 = $arcB->Y_offset($t2,$offB,$spanB->[0]->[1],$spanB->[3]);

                    return $yoff2 - $yoff1;
                };

                my $bounds_xoff = [$spanx->[0],$spanx->[1]];
                my ($int_xoff,$msg)=BrentsMethod($findintbezarc,$bounds_xoff,0.00001,undef,'subBezOff-subArcOff intersection finding');

                if ($msg) { warn "bezarcoffint message: $msg\n"; }
                else {
                    my $t1 = $bezA->t_from_xoff($int_xoff,$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1]);
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
                #push @ret, ["xoverlap","yoverlap"] if $bezA != $arcB;
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

                    # for the bezier
                    my $tA = $bezA->t_from_xoff($_[0],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1]);
                    my $xA = $bezA->bezierEvalXofT($tA);
                    my $tA_prime  = $spanA->[0]->[1]->($xA);
                    my $tA_2prime = $spanA->[0]->[2]->($xA);
                    #bez->f_prime()
                    my $mA = $bezA->{Em3} * $tA**2 * $tA_prime  +  $bezA->{Fm2} * $tA * $tA_prime + $bezA->{G} * $tA_prime;
                    #bez->f_2prime()
                    my $mA_prime = $bezA->{Em3} * ($tA**2 * $tA_2prime + 2*$tA*$tA_prime**2) + $bezA->{Fm2} * ($tA * $tA_2prime + $tA_prime**2) + $bezA->{G} * $tA_2prime;
                    my $YPrimeA = $bezA->bezierEvalYPrimeofT($tA);
                    my $yoffset_primeA = -($offA/2.0) * (1/(sqrt($mA**2 + 1)**3)) * 2*$mA * $mA_prime;
                    $YPrimeA *= -1 if $spanA->[3];
                    $yoffset_primeA *= -1 if $spanA->[3];
                    my $YoffPrimeA = ($YPrimeA + $yoffset_primeA);

                    # for the arc
                    my $tB = $arcB->t_from_xoff($_[0],$offB,[$spanB->[2]->[0],$spanB->[2]->[-1]],$spanB->[0]->[1]);
                    my $xB = $arcB->evalXofTheta($tB);
                    my $tB_prime  = $spanB->[0]->[1]->($xB);
                    my $tB_2prime = $spanB->[0]->[2]->($xB);
                    # arc->f_prime()
                    my $mB = $arcB->evalYPrimeofTheta($tB) * $tB_prime;
                    # arc->f_2prime()
                    my $mB_prime = $arcB->evalYDoublePrimeofTheta($tB) * $tB_prime**2 + $arcB->evalYPrimeofTheta($tB) * $tB_2prime;
                    my $YPrimeB = $arcB->evalYPrimeofTheta($tB);
                    my $yoffset_primeB = -($offB/2.0) * 1/(sqrt($mB**2 + 1)**3) * 2*$mB * $mB_prime;
                    $YPrimeB *= -1 if $spanB->[3];
                    $yoffset_primeB *= -1 if $spanB->[3];
                    my $YoffPrimeB = ($YPrimeB + $yoffset_primeB);

                    my $ret = $YoffPrimeA - $YoffPrimeB;

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

                    $split_t_A = $bezA->t_from_xoff($span_split_x,$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1]);
                    $split_t_B = $arcB->t_from_xoff($span_split_x,$offB,[$spanB->[2]->[0],$spanB->[2]->[-1]],$spanB->[0]->[1]);

                    # set up new LUTs to pass to a re-call of this intersection sub
                    # so re-call run will have proper x bounds to find each of
                    # the intersection pair individually.

                    my $sub_t_span_1_A = [$tsA->[$spanA->[3]?1:0], $split_t_A];
                    my $sub_t_span_2_A = [$split_t_A, $tsA->[$spanA->[3]?0:1]];
                    my $sub_t_span_1_B = [$tsB->[$spanB->[3]?1:0], $split_t_B];
                    my $sub_t_span_2_B = [$split_t_B, $tsB->[$spanB->[3]?0:1]];

                    #warn "re-call 1\n";
                    #warn "  [[$spanA->[0],[$spanx->[0]  , $span_split_x], $spanA->[3]?$sub_t_span_2_A:$sub_t_span_1_A, $spanA->[3]]]\n";
                    #warn "  [[$spanB->[0],[$spanx->[0]  , $span_split_x], $spanB->[3]?$sub_t_span_2_B:$sub_t_span_1_B, $spanB->[3]]]\n";
                    my @int1 = intersect_CoAo($bezA,$arcB,$offA,$offB,
                                                 [[$spanA->[0],[$spanx->[0]  , $span_split_x], $spanA->[3]?$sub_t_span_2_A:$sub_t_span_1_A, $spanA->[3]]],
                                                 [[$spanB->[0],[$spanx->[0]  , $span_split_x], $spanB->[3]?$sub_t_span_2_B:$sub_t_span_1_B, $spanB->[3]]]
                                                 );
                    #warn "re-call 2\n";
                    #warn "  [[$spanA->[0],[$span_split_x, $spanx->[1]  ], $spanA->[3]?$sub_t_span_1_A:$sub_t_span_2_A, $spanA->[3]]]\n";
                    #warn "  [[$spanB->[0],[$span_split_x, $spanx->[1]  ], $spanB->[3]?$sub_t_span_1_B:$sub_t_span_2_B, $spanB->[3]]]\n";

                    my @int2 = intersect_CoAo($bezA,$arcB,$offA,$offB,
                                                 [[$spanA->[0],[$span_split_x, $spanx->[1]  ], $spanA->[3]?$sub_t_span_1_A:$sub_t_span_2_A, $spanA->[3]]],
                                                 [[$spanB->[0],[$span_split_x, $spanx->[1]  ], $spanB->[3]?$sub_t_span_1_B:$sub_t_span_2_B, $spanB->[3]]]
                                                 );

                warn "expecting two intersections but got [",(scalar(@int1)," + ",scalar(@int2)),"]" if (scalar(@int1) != 1 || scalar(@int2) != 1);

                push @ret, $int1[0]->[2] <= $int2[0]->[2]
                           ? ($int1[0],$int2[0])
                           : ($int2[0],$int1[0]);
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

    $offA //= (exists($bezA->{offset})?$bezA->{offset}:0);
    $offB //= (exists($bezB->{offset})?$bezB->{offset}:0);
    return intersect_CC($bezA, $bezB, $lutA, $lutB) if (!$offA && !$offB);

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


                # TODO
                # Need to confirm this in test cases that push the limits of this math
                # - need to work out the math for where those limits are, to make those test cases.

                my $y_diff_prime = sub {
                    my $tA = $bezA->t_from_xoff($_[0],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);
                    my $tB = $bezB->t_from_xoff($_[0],$offB,[$spanB->[2]->[0],$spanB->[2]->[-1]],$spanB->[0]->[1],$spanB->[3]);

                    my $xA = $bezA->bezierEvalXofT($tA);
                    my $tA_prime  = $spanA->[0]->[1]->($xA);
                    my $tA_2prime = $spanA->[0]->[2]->($xA);
                    my $xB = $bezB->bezierEvalXofT($tB);
                    my $tB_prime  = $spanB->[0]->[1]->($xB);
                    my $tB_2prime = $spanB->[0]->[2]->($xB);

                    #bez->f_prime()
                    my $mA = $bezA->{Em3} * $tA**2 * $tA_prime  +  $bezA->{Fm2} * $tA * $tA_prime + $bezA->{G} * $tA_prime;
                    my $mB = $bezB->{Em3} * $tB**2 * $tB_prime  +  $bezB->{Fm2} * $tB * $tB_prime + $bezB->{G} * $tB_prime;

                    #bez->f_2prime()
                    my $mA_prime = $bezA->{Em3} * ($tA**2 * $tA_2prime + 2*$tA*$tA_prime**2) + $bezA->{Fm2} * ($tA * $tA_2prime + $tA_prime**2) + $bezA->{G} * $tA_2prime;
                    my $mB_prime = $bezB->{Em3} * ($tB**2 * $tB_2prime + 2*$tB*$tB_prime**2) + $bezB->{Fm2} * ($tB * $tB_2prime + $tB_prime**2) + $bezB->{G} * $tB_2prime;

                    my $YPrimeA = $bezA->bezierEvalYPrimeofT($tA);
                    my $YPrimeB = $bezB->bezierEvalYPrimeofT($tB);

                    my $yoffset_primeA = -($offA/2.0) * (1/(sqrt($mA**2 + 1)**3)) * 2*$mA * $mA_prime;
                    my $yoffset_primeB = -($offB/2.0) * (1/(sqrt($mB**2 + 1)**3)) * 2*$mB * $mB_prime;
 
                    # if isReversed flag true, (increasing x gives decreasing t)
                    # dY/dt and dYoffset/dt need to be negated,
                    # to correspond to dY/dx
                    $YPrimeA *= -1 if $spanA->[3];
                    $YPrimeB *= -1 if $spanB->[3];
                    $yoffset_primeA *= -1 if $spanA->[3];
                    $yoffset_primeB *= -1 if $spanB->[3];

                    my $YoffPrimeA = ($YPrimeA + $yoffset_primeA);
                    my $YoffPrimeB = ($YPrimeB + $yoffset_primeB);

                    my $ret = $YoffPrimeB - $YoffPrimeA;

                    return $ret;

                };


                my ($at_start,$at_end);

                $at_start = $y_diff_prime->($spanx->[0]);
                $at_end   = $y_diff_prime->($spanx->[1]);

                #warn "\nydiffprime start, end: [$tsA->[0]] $at_start, [$tsA->[1]] $at_end\n\n";

                if ($at_start > 0 && $at_end < 0 || $at_start < 0 && $at_end > 0) {

                    my ($span_split_x,$split_t_A,$msg);

                    my $bounds_xoffA = [$spanx->[0],$spanx->[1]];
                    ($span_split_x,$msg)=FalsePosition($y_diff_prime,$bounds_xoffA,0.00001,(($bounds_xoffA->[1]+$bounds_xoffA->[0])/2),'subBezOffset-subBezOffset intersection finding - find pair split parameter');


                    die "split find fail msg: $msg\n" if $msg;

                    #warn "split xoff: $span_split_x\n";

                    $split_t_A = $bezA->t_from_xoff($span_split_x,$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);
                    my $split_t_B = $bezB->t_from_xoff($span_split_x,$offB,[$spanB->[2]->[0],$spanB->[2]->[-1]],$spanB->[0]->[1],$spanB->[3]);

                    # set up new LUTs to pass to a re-call of this intersection sub
                    # so re-call run will have proper x bounds to find each of
                    # the intersection pair individually.

                    my $sub_t_span_1_A = [$tsA->[$spanA->[3]?1:0], $split_t_A];
                    my $sub_t_span_2_A = [$split_t_A, $tsA->[$spanA->[3]?0:1]];
                    my $sub_t_span_1_B = [$tsB->[$spanB->[3]?1:0], $split_t_B];
                    my $sub_t_span_2_B = [$split_t_B, $tsB->[$spanB->[3]?0:1]];

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

                warn "expecting two intersections but got [",(scalar(@int1)," + ",scalar(@int2)),"]" if (scalar(@int1) != 1 || scalar(@int2) != 1);

                push @ret, $int1[0]->[2] <= $int2[0]->[2]
                           ? ($int1[0],$int2[0])
                           : ($int2[0],$int1[0]);
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
    
    $offA //= (exists($bezA->{offset})?$bezA->{offset}:0);
    $offB //= (exists($lineB->{offset})?$lineB->{offset}:0);
    return intersect_CL($bezA, $lineB, $lutA, $lutB) if (!$offA && !$offB);

    if (!$lutB && !exists $lineB->{XtoTLUT}) {
        # too complex probably for line segments, but hack it in here while adapting this bez-bez algo to work with bez-line
        my $line_is_reversed = $lineB->{p2}->[0] >= $lineB->{p1}->[0] ? 0 : 1;
        $lineB->{XtoTLUT} = [
                             [
                              [sub {$lineB->solveXforTheta($_[0])}], # bez would have 3 or 6 sub refs here. line only needs 1 or 2. 1 for our purposes here so far
                              [$lineB->{minx},$lineB->{maxx}],
                              [$line_is_reversed?(1,0):(0,1)],
                              $line_is_reversed
                             ]
                            ];
    }

    my $XtoTLUT_A = $lutA ? $lutA : [map [$_->[0],[$_->[1]->[0],$_->[1]->[-1]],[$_->[2]->[0],$_->[2]->[-1]],$_->[3]], @{$bezA->{XtoTLUT}}];
    my $XtoTLUT_B = $lutB ? $lutB : [map [$_->[0],[$_->[1]->[0],$_->[1]->[-1]],[$_->[2]->[0],$_->[2]->[-1]],$_->[3]], @{$lineB->{XtoTLUT}}];

    # assume any provided LUTs already contain offset Xs in their x range entries
    # otherwise we need to do that now for the LUT copies we made
    #my $lcnt=0;
    if (!defined($lutA)) {
        foreach my $span (@{$XtoTLUT_A}) {
            #warn "lcnt: $lcnt\n";
            #$lcnt++;
            #warn "hello colo1 [$bezA->{p1}->[0],$bezA->{p1}->[1]],[$bezA->{cp1}->[0],$bezA->{cp1}->[1]],[$bezA->{cp2}->[0],$bezA->{cp2}->[1]],[$bezA->{p2}->[0],$bezA->{p2}->[1]], $span->[2]->[0],$offA\n";
            $span->[1]->[0]  = $bezA->X_offset($span->[2]->[0] ,$offA,$span->[0]->[1],$span->[3]);
            #warn "hello colo2 [$bezA->{p1}->[0],$bezA->{p1}->[1]],[$bezA->{cp1}->[0],$bezA->{cp1}->[1]],[$bezA->{cp2}->[0],$bezA->{cp2}->[1]],[$bezA->{p2}->[0],$bezA->{p2}->[1]], $span->[2]->[-1],$offA\n";
            $span->[1]->[-1] = $bezA->X_offset($span->[2]->[-1],$offA,$span->[0]->[1],$span->[3]);
            #warn "hello colo3\n";
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

                # for the Bezier
                my $tA = $bezA->t_from_xoff($_[0],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);
                my $xA = $bezA->bezierEvalXofT($tA);
                my $tA_prime  = $spanA->[0]->[1]->($xA);
                my $tA_2prime = $spanA->[0]->[2]->($xA);
                # bez->f_prime()
                my $mA = $bezA->{Em3} * $tA**2 * $tA_prime  +  $bezA->{Fm2} * $tA * $tA_prime + $bezA->{G} * $tA_prime;
                # bez->f_2prime()
                my $mA_prime = $bezA->{Em3} * ($tA**2 * $tA_2prime + 2*$tA*$tA_prime**2) + $bezA->{Fm2} * ($tA * $tA_2prime + $tA_prime**2) + $bezA->{G} * $tA_2prime;

                my $YPrimeA = $bezA->bezierEvalYPrimeofT($tA);
                my $yoffset_primeA = -($offA/2.0) * (1/(sqrt($mA**2 + 1)**3)) * 2*$mA * $mA_prime;

                $YPrimeA *= -1 if $spanA->[3];
                $yoffset_primeA *= -1 if $spanA->[3];

                my $YoffPrimeA = ($YPrimeA + $yoffset_primeA);

                # for the Line
                my $YoffPrimeB = $lineB->{dy}; # dy/dt which is (y2-y1)/(1-0)
                $YoffPrimeB *= -1 if $spanB->[3];

                my $ret = $YoffPrimeA - $YoffPrimeB;

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

                my $sub_t_span_1_A = [$tsA->[$spanA->[3]?1:0], $split_t_A];
                my $sub_t_span_2_A = [$split_t_A, $tsA->[$spanA->[3]?0:1]];
                my $sub_t_span_1_B = [$tsB->[$spanB->[3]?1:0], $split_t_B];
                my $sub_t_span_2_B = [$split_t_B, $tsB->[$spanB->[3]?0:1]];

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

                warn "expecting two intersections but got [",(scalar(@int1)," + ",scalar(@int2)),"]" if (scalar(@int1) != 1 || scalar(@int2) != 1);

                push @ret, $int1[0]->[2] <= $int2[0]->[2]
                           ? ($int1[0],$int2[0])
                           : ($int2[0],$int1[0]);
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

} # end package
1;

