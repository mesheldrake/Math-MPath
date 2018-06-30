####################################################################################
###      Math::MPath::EllipticalArc           ###################################
# http://www.w3.org/TR/SVG/implnote.html#ArcImplementationNotes
package Math::MPath::EllipticalArc;
{
use Math::MPath::QuadraticFormula;
use Math::MPath::Function::Root qw(BrentsMethod FalsePosition);
our $pi = 4 * atan2(1,1);

sub new {
    my $class = shift;
    my $self={};
    bless $self,$class;
    $self->{p1} = shift;
    $self->{r} = shift;
    $self->{rx} = abs($self->{r}->[0]);
    $self->{ry} = abs($self->{r}->[1]);
    $self->{phi} = shift;
    # angle starts as degrees, mod 360, per svg 1.1 implementation notes
    while ($self->{phi} >  360) {$self->{phi} -= 360;}
    while ($self->{phi} < -360) {$self->{phi} += 360;}
    $self->{phi_radians} = $self->{phi} * ($pi/180);
    $self->{large_arc_flag} = shift;
    $self->{sweep_flag} = shift;
    if ($self->{large_arc_flag} ne 0) {$self->{large_arc_flag}=1;}
    if ($self->{sweep_flag} ne 0) {$self->{sweep_flag}=1;}
    $self->{p2} = shift;
    $self->{precision} = shift;
    $self->{isLite} = shift;
    $self->{maxdiaglength} = sqrt(($self->{maxx} - $self->{minx})**2 + ($self->{maxy} - $self->{miny})**2);

    #calculate the center of the ellipse
    #step 1: Compute (x1', y1')
    my $x1prime=(($self->{p1}->[0] - $self->{p2}->[0])/2) *  cos($self->{phi_radians})  +  (($self->{p1}->[1] - $self->{p2}->[1])/2) * sin($self->{phi_radians});
    my $y1prime=(($self->{p1}->[0] - $self->{p2}->[0])/2) * -sin($self->{phi_radians})  +  (($self->{p1}->[1] - $self->{p2}->[1])/2) * cos($self->{phi_radians});
    #make sure they are large enough to make an ellipse that will reach the destination point
    if    ($self->{rx}==0 || $self->{rx}==0) {warn "elliptical arc with zero x or y radius. x,y:",$self->{p2}->[0],", ",$self->{p2}->[1],", rx,ry: $self->{rx}, $self->{ry}\n";}
    # avoid divide by zero and force $lam to be too big if rx or ry == 0
    # ... but then you're just kicking the problem down to later divide by zero cases.
    # TODO: the issue is, what do you do with zero-radius arcs? let them exist? If so, should probably skip over all this attempt to calc with 0.
    #       what about zero length line segs, zero length beziers?
    #       But in those cases and this, the endpoints should give that away. Like, zero radius must be secondary to p1 and p2 being the same.
    #       Diff here would be with possible full arc with largearcflag=1.
    #       Hmm.
    my $lam = ($self->{rx} == 0 || $self->{ry} == 0 ) ? 2 : ($x1prime)**2/($self->{rx}**2) + ($y1prime)**2/($self->{ry}**2);
    if ($lam > 1) {
        my $sqrtlam = sqrt($lam);
        $self->{rx} *= $sqrtlam + 0.000000001;
        $self->{ry} *= $sqrtlam + 0.000000001;
    }
    #step 2: Compute (cX ', cY  ')

    my $to_sqrt_for_hairy_radical;
    if ($x1prime == 0 && $y1prime == 0) {
        $to_sqrt_for_hairy_radical=0; # otherwise we get divide by zero
    }
    else {
        $to_sqrt_for_hairy_radical = ( ($self->{rx}**2 * $self->{ry}**2) - ($self->{rx}**2 * $y1prime**2) - ($self->{ry}**2 * $x1prime**2) )/( ($self->{rx}**2 * $y1prime**2) + ($self->{ry}**2 * $x1prime**2) );
    }

    if ($to_sqrt_for_hairy_radical<0 && abs($to_sqrt_for_hairy_radical)<10**-15) {$to_sqrt_for_hairy_radical=0;} #had some really small negative - like 10**-19 - choking sqrt. Maybe should find how that happened. But snapping to zero for now to get work done.
    #my $hairy_radical = (($self->{large_arc_flag} eq $self->{sweep_flag})?-1:1) * sqrt(( ($self->{rx}**2 * $self->{ry}**2) - ($self->{rx}**2 * $y1prime**2) - ($self->{ry}**2 * $x1prime**2) )/( ($self->{rx}**2 * $y1prime**2) + ($self->{ry}**2 * $x1prime**2) ));
    my $hairy_radical = (($self->{large_arc_flag} eq $self->{sweep_flag})?-1:1) * sqrt($to_sqrt_for_hairy_radical);
    my $cxprime = $hairy_radical *      (($self->{rx} * $y1prime)/$self->{ry});
    my $cyprime = $hairy_radical * -1 * (($self->{ry} * $x1prime)/$self->{rx});
    #step 3: Compute (cX, cY) from (cX ', cY  ')
    $self->{cx} = $cxprime * cos($self->{phi_radians}) + $cyprime * -sin($self->{phi_radians}) + ($self->{p1}->[0] + $self->{p2}->[0])/2;
    $self->{cy} = $cxprime * sin($self->{phi_radians}) + $cyprime *  cos($self->{phi_radians}) + ($self->{p1}->[1] + $self->{p2}->[1])/2;
    #Step 4: Compute theta1     and     delta-theta
    #my $theta1_arccos_arg = (($y1prime - $cyprime)/$self->{ry})/sqrt((($y1prime - $cyprime)/$self->{ry})**2 + (($x1prime - $cxprime)/$self->{rx})**2);
    my $theta1_arccos_arg;
    if ($x1prime == 0 && $cxprime == 0) {
        $theta1_arccos_arg =0;
    }
    else {
        $theta1_arccos_arg =  (($x1prime - $cxprime)/$self->{rx})/sqrt((($x1prime - $cxprime)/$self->{rx})**2 + (($y1prime - $cyprime)/$self->{ry})**2);
    }
    my $theta1sign = (($y1prime-$cyprime)/$self->{ry}) < 0 ? -1:1; # - ((0)*(...))   ux*vy-uy*vx;
    $self->{theta1} = $theta1sign * ($pi/2 - asin($theta1_arccos_arg));

    my $delta_theta_arccos_arg;

    if ( $x1prime == 0 && $cxprime == 0 && $y1prime == 0 && $cyprime == 0) {
        $delta_theta_arccos_arg =0;
    }
    else {
        $delta_theta_arccos_arg =  ((($x1prime - $cxprime)/$self->{rx}) * ((-$x1prime - $cxprime)/$self->{rx}) + (($y1prime - $cyprime)/$self->{ry}) * ((-$y1prime - $cyprime)/$self->{ry}))/sqrt((($x1prime - $cxprime)/$self->{rx})**2 + (($y1prime - $cyprime)/$self->{ry})**2);
    }

    my $delta_theta_sign=(
        ( (($x1prime - $cxprime)/$self->{rx}) * ((-$y1prime - $cyprime)/$self->{ry}) )
            -
        ( (($y1prime - $cyprime)/$self->{ry}) * ((-$x1prime - $cxprime)/$self->{rx}) )

        < 0
    ) ? -1 : 1;

    $self->{delta_theta} = $delta_theta_sign * (($pi/2) - asin($delta_theta_arccos_arg));

    # Make all these theta mods sensible, 
    # in accordance with what's supposed to happen, 
    # overcoming any angle flattenning that comes from atan sqrt usage above.

    #mod(360deg) for delta_theta
    if (abs($self->{delta_theta}) > 2*$pi) {
        my $div=$self->{delta_theta}/(2*$pi);
        my $rem=$div - int($div);
        $self->{delta_theta} = $rem * (2*$pi);
    }

    if ($self->{sweep_flag}) {
        if ($self->{delta_theta} < 0) {
            if (!$self->{large_arc_flag}) {$self->{delta_theta} *= -1 ;}
            else {$self->{delta_theta} += 2 * $pi;}
            }
    }
    else {
        if ($self->{delta_theta} > 0) {$self->{delta_theta} -= 2 * $pi;}
        if ($self->{large_arc_flag} && $self->{delta_theta} > -$pi) {$self->{delta_theta} = -(2*$pi+$self->{delta_theta});}
    }

    $self->{theta2}=$self->{theta1} + $self->{delta_theta};

    #calculate the foci of the ellipse
    $self->{f1}= $self->{rx} > $self->{ry} ? [sqrt($self->{rx}**2 - $self->{ry}**2),0] : [0,sqrt($self->{ry}**2 - $self->{rx}**2)];
    $self->{f2}= $self->{rx} > $self->{ry} ? [-$self->{f1}->[0],$self->{f1}->[1]] : [$self->{f1}->[0],-$self->{f1}->[1]];
    #now is a good time to calculate eccentricity, too - used in circumference and arc calculations
    $self->{eccentricity}=$self->{f1}/(($self->{rx}>$self->{ry})?$self->{rx}:$self->{ry});
    $self->{f1} = [($self->{f1}->[0]*cos($self->{phi_radians}) - $self->{f1}->[1]*sin($self->{phi_radians})),
                   ($self->{f1}->[0]*sin($self->{phi_radians}) + $self->{f1}->[1]*cos($self->{phi_radians}))];
    $self->{f1}->[0]+=$self->{cx};
    $self->{f1}->[1]+=$self->{cy};
    $self->{f2} = [($self->{f2}->[0]*cos($self->{phi_radians}) - $self->{f2}->[1]*sin($self->{phi_radians})),
                   ($self->{f2}->[0]*sin($self->{phi_radians}) + $self->{f2}->[1]*cos($self->{phi_radians}))];

    $self->{f2}->[0]+=$self->{cx};
    $self->{f2}->[1]+=$self->{cy};

    if (!$self->{isLite}) {
        $self->initBigs();
    }

    my @extremexs_is=(0,((!$self->{isLite})?$self->solveXPrimeforThetaBig(Math::BigFloat->bzero()):$self->solveXPrimeforTheta(0)),1);
    my @extremexs = map {(ref($_) && !$self->{isLite}) ? $self->evalXofThetaBig($_):$self->evalXofTheta($_)} @extremexs_is;
    my @extremexs_sorted = sort {$a<=>$b} @extremexs;

    my @extremeys_is=(0,((!$self->{isLite})?$self->solveYPrimeforThetaBig(Math::BigFloat->bzero()):$self->solveYPrimeforTheta(0)),1);
    my @extremeys = map {(ref($_) && !$self->{isLite}) ? $self->evalYofThetaBig($_):$self->evalYofTheta($_)} @extremeys_is;
    my @extremeys_sorted = sort {$a<=>$b} @extremeys;

    $self->{extremexs_is}=[(@extremexs_is)];
    $self->{extremeys_is}=[(@extremeys_is)];
    $self->{extremexs}=[(@extremexs)];
    $self->{extremeys}=[(@extremeys)];
    $self->{minx} = $extremexs_sorted[0];
    $self->{maxx} = $extremexs_sorted[$#extremexs];
    $self->{miny} = $extremeys_sorted[0];
    $self->{maxy} = $extremeys_sorted[$#extremeys];


    if (!$self->{isLite}) {
        $self->initDangerRanges();
    }

    # At some point you stopped making these
    # but they're still refered to in f() and F().
    # And there are similar/same ranges you calculate
    # when initDangerRanges() is run, so maybe you meant to
    # switch to using those?
    $self->{Fydangerranges} = [];
    $self->{fxdangerranges} = [];

    # Following not ideal but workable approach developed for cubic Bezier,
    # make a LUT of x spans where the elliptical arc is monotonic, where
    # we can assign definite t(x) and t(y) one-to-one functions, to facilitate
    # intersections and offset intersections with other curves, using the
    # same approach worked out for intersections of two cubic Beziers
    # in Intersections.pm

    my @div_angles;
    my $n = int($self->{theta1}/($pi/2)); # multiple of pi/2 that gives the axis angle _before_ theta1, when sweeping away from 0 radians in same direction as delta_theta; [-1,0,1]
    for (0 .. abs(int($self->{delta_theta}/($pi/2)))) { # should loop 1 to 4 times
        $n += $self->{delta_theta} > 0 ? 1 : -1; # now a multiple of pi/2 that gives an axis angle _after_ theta1, when sweeping away from 0 radians in same direction as delta_theta; can be as high/low as +/-5 if n started at +/-1
        my $axis_angle = $n*($pi/2);
        my $axis_delta_angle = $axis_angle - $self->{theta1};
        # this is only potentially false on the last loop
        if (abs($axis_delta_angle) < abs($self->{delta_theta})) {
            push @div_angles, $axis_angle;
        }
    }

    my @div_ts = (0,1);

    push @div_ts, map $self->t_of_theta($_), @div_angles;

    if ($self->{phi} ne 0) {
        push @div_ts, ((!$self->{isLite})?$self->solveXPrimeforThetaBig(Math::BigFloat->bzero()):$self->solveXPrimeforTheta(0));
        push @div_ts, ((!$self->{isLite})?$self->solveYPrimeforThetaBig(Math::BigFloat->bzero()):$self->solveYPrimeforTheta(0));
    }

    @div_ts = sort {$a<=>$b} @div_ts;

    # Create the same kind of sort of clunky LUT we worked up in cubic Bezier case
    # so we can work with these elliptical arc sub segments in the same way -
    # especially in subsegment intersections. Later we'll clean and simplify how
    # all that gets set up there and here.

    my @XtoTLUT; # [ [t(x),t'(x),t''(x),t(y),t'(y),t''(y)],
                 #   [xlow,xhigh], # in left-to-right order (as opposed to t order, which may or may not be the same)
                 #   [t corresponding to xlow,t corresponding to xhigh],
                 #   [isReversed - true if t0 corresponds to xhigh]
                 #   [ylow,yhigh] # in numerical lower to higher order (not visual low to high, because y-axis orientation might confuse the issue)
                 # ]

    my $t_of_x_eqn_1and2 = sub {
        my ($x,$which_theta) = @_;
        my $x_unshift = $x - $self->{cx};
        my $x_u = $x_unshift * cos(-$self->{phi_radians});
        my $y_u = $x_unshift * sin(-$self->{phi_radians});
        my $m = sin($pi/2 - $self->{phi_radians})/cos($pi/2 - $self->{phi_radians});
        my $y1_min_mx1 = $y_u-$m*$x_u;
        my $a =        $self->{ry}  / $y1_min_mx1;
        my $b = -($m * $self->{rx}) / $y1_min_mx1;
        my $phase_angle = atan2( $b , $a );
        my $theta  = $which_theta
                   ? asin(1/ sqrt($a**2 + $b**2)) - $phase_angle
                   : asin(1/-sqrt($a**2 + $b**2)) - $phase_angle + $pi # + pi by educated trial and error
                   ;
        if ($self->{sweep_flag}) {while ($theta < $self->{theta1}) {$theta += 2*$pi;}}
        else                     {while ($theta > $self->{theta1}) {$theta -= 2*$pi;}}
        my $t = $self->t_of_theta($theta);
        return $t;
    };
    my $t_of_y_eqn_1and2 = sub {
        my ($y,$which_theta) = @_;
        my $y_unshift = $y - $self->{cy};
        my $x_u = $y_unshift * -sin(-$self->{phi_radians});
        my $y_u = $y_unshift *  cos(-$self->{phi_radians});
        my $m = sin(-$self->{phi_radians})/cos(-$self->{phi_radians});
        my $y1_min_mx1 = $y_u-$m*$x_u;
        my $a =        $self->{ry}  / $y1_min_mx1;
        my $b = -($m * $self->{rx}) / $y1_min_mx1;
        my $phase_angle = atan2( $b , $a );
        my $theta  = $which_theta
                   ? asin(1/ sqrt($a**2 + $b**2)) - $phase_angle
                   : asin(1/-sqrt($a**2 + $b**2)) - $phase_angle + $pi # + pi by educated trial and error
                   ;
        if ($self->{sweep_flag}) {while ($theta < $self->{theta1}) {$theta += 2*$pi;}}
        else                     {while ($theta > $self->{theta1}) {$theta -= 2*$pi;}}
        my $t = $self->t_of_theta($theta);
        return $t;
    };

    my $t_prime_of_x_eqn_1and2 = sub {
        my ($x,$which) = @_;
        my $x_unshift = $x - $self->{cx};
        my $x_u = $x_unshift * cos(-$self->{phi_radians});
        my $y_u = $x_unshift * sin(-$self->{phi_radians});
        my $m = sin($pi/2 - $self->{phi_radians})/cos($pi/2 - $self->{phi_radians});
        my $z = (($self->{ry}**2+($m*$self->{rx})**2)/((1+$m**2)*sin(-$self->{phi_radians}+atan2(-$m,1))**2));
        my $theta_prime = ($which?1:-1) / sqrt($z - $x_unshift**2);
        my $t_prime = $theta_prime / $self->{delta_theta};
        return $t_prime;
    };

    # just copy-paste and tried to adapt from x case
    # needs testing:
    my $t_prime_of_y_eqn_1and2 = sub {
        my ($y,$which) = @_;
        my $y_unshift = $y - $self->{cy};
        my $x_u = $y_unshift * -sin(-$self->{phi_radians}); # cos(pi/2 - phi)
        my $y_u = $y_unshift *  cos(-$self->{phi_radians}); # sin(pi/2 - phi)
        my $m = sin(-$self->{phi_radians})/cos(-$self->{phi_radians});
        my $z = (($self->{ry}**2+($m*$self->{rx})**2)/((1+$m**2)*sin( ( ($pi/2) - $self->{phi_radians}) + atan2(-$m,1))**2));
        my $theta_prime = ($which?1:-1) / sqrt($z - $y_unshift**2);
        my $t_prime = $theta_prime / $self->{delta_theta};
        return $t_prime;
    };

    # needs testing
    my $t_2prime_of_x_eqn_1and2 = sub {
        my ($self,$x,$which) = @_;
        my $x_unshift = $x - $self->{cx};
        my $x_u = $x_unshift * cos(-$self->{phi_radians});
        my $y_u = $x_unshift * sin(-$self->{phi_radians});
        my $m = sin($pi/2 - $self->{phi_radians})/cos($pi/2 - $self->{phi_radians});
        my $z = (($self->{ry}**2+($m*$self->{rx})**2)/((1+$m**2)*sin(-$self->{phi_radians}+atan2(-$m,1))**2));
        my $theta_2prime = ($which?1:-1) * $x_unshift / (sqrt($z - $x_unshift**2)**3);
        my $t_2prime = $theta_2prime / $self->{delta_theta};
        return $t_2prime;
    };

    # needs testing
    my $t_2prime_of_y_eqn_1and2 = sub {
        my ($y,$which) = @_;
        my $y_unshift = $y - $self->{cy};
        my $x_u = $y_unshift * -sin(-$self->{phi_radians}); # cos(pi/2 - phi)
        my $y_u = $y_unshift *  cos(-$self->{phi_radians}); # sin(pi/2 - phi)
        my $m = sin(-$self->{phi_radians})/cos(-$self->{phi_radians});
        my $z = (($self->{ry}**2+($m*$self->{rx})**2)/((1+$m**2)*sin( ( ($pi/2) - $self->{phi_radians}) + atan2(-$m,1))**2));
        my $theta_2prime = ($which?1:-1) * $y_unshift / (sqrt($z - $y_unshift**2)**3);
        my $t_2prime = $theta_2prime / $self->{delta_theta};
        return $t_2prime;
    };

    for (my $i = 1; $i < @div_ts; $i++) {
        my $ta = $div_ts[$i-1];
        my $tb = $div_ts[$i];
        my $tmid = ($ta + $tb) / 2;

        my $xa = $self->evalXofTheta($ta);
        my $xb = $self->evalXofTheta($tb);
        my $ya = $self->evalYofTheta($ta);
        my $yb = $self->evalYofTheta($tb);
        my $isReversed = $xa > $xb;

        my $angle_mid = $self->theta_of_t($tmid);

        while ($angle_mid >  $pi) {$angle_mid -= 2*$pi;}
        while ($angle_mid < -$pi) {$angle_mid += 2*$pi;}

        my @tx_eqns;
        my @ty_eqns;

        # The "which" flag in these calls to t_of_x_eqn_1and2() and
        # t_of_y_eqn_1and2() was figured by trial and error.

        if ($angle_mid > 0) {
            if ($angle_mid < $pi/2) { # quadrant 1
                push @tx_eqns, (
                    sub {#t(x)
                         #warn "(eq a)\n";
                         return $t_of_x_eqn_1and2->($_[0],    0    );
                    },
                    sub {#t_prime(x)
                        #warn "(eq a)\n";
                        return $t_prime_of_x_eqn_1and2->($_[0],    0    );
                    },
                    sub {#t_2prime(x)
                        #warn "(eq a)\n";
                        return $t_2prime_of_x_eqn_1and2->($_[0],    0    );
                    },
                    sub {#t(y)
                         #warn "(eq a)\n";
                         return $t_of_y_eqn_1and2->($_[0],    1    );
                    },
                    sub {#t'(y)
                         #warn "(eq a)\n";
                        return $t_prime_of_y_eqn_1and2->($_[0],    1    );
                    },
                    sub {#t_2prime(y)
                        #warn "(eq a)\n";
                        return $t_2prime_of_y_eqn_1and2->($_[0],    1    );
                    }
                );
            }
            elsif ($angle_mid < $pi) { # quadrant 2
                push @tx_eqns, (
                    sub {
                         #warn "(eq b)\n";
                         return $t_of_x_eqn_1and2->($_[0],    1    );
                    },
                    sub {#t_prime(x)
                        #warn "(eq b)\n";
                        return $t_prime_of_x_eqn_1and2->($_[0],    1    );
                    },
                    sub {#t_2prime(x)
                        #warn "(eq b)\n";
                        return $t_2prime_of_x_eqn_1and2->($_[0],    1    );
                    },
                    sub {
                         #warn "(eq b)\n";
                         return $t_of_y_eqn_1and2->($_[0],    0    );
                    },
                    sub {#t'(y)
                         #warn "(eq a)\n";
                        return $t_prime_of_y_eqn_1and2->($_[0],    0    );
                    },
                    sub {#t_2prime(y)
                        #warn "(eq b)\n";
                        return $t_2prime_of_y_eqn_1and2->($_[0],    0    );
                    }
                );
            }
            else { warn "out of bounds angle [$angle_mid]";}
        }
        elsif ($angle_mid < 0) {
            if ($angle_mid > -$pi/2) { # quadrant 4
                push @tx_eqns, (
                    sub {
                         #warn "(eq c)\n";
                         return $t_of_x_eqn_1and2->($_[0],    1    );
                    },
                    sub {#t_prime(x)
                        #warn "(eq c)\n";
                        return $t_prime_of_x_eqn_1and2->($_[0],    1    );
                    },
                    sub {#t_2prime(x)
                        #warn "(eq c)\n";
                        return $t_2prime_of_x_eqn_1and2->($_[0],    1    );
                    },
                    sub {
                         #warn "(eq c)\n";
                         return $t_of_y_eqn_1and2->($_[0],    0    );
                    },
                    sub {#t'(y)
                         #warn "(eq a)\n";
                        return $t_prime_of_y_eqn_1and2->($_[0],    0    );
                    },
                    sub {#t_2prime(y)
                        #warn "(eq c)\n";
                        return $t_2prime_of_y_eqn_1and2->($_[0],    0    );
                    }
                );
            }
            elsif ($angle_mid > -$pi) { # quadrant 3
                push @tx_eqns, (
                    sub {
                         #warn "(eq d)\n";
                         return $t_of_x_eqn_1and2->($_[0],    0    );
                    },
                    sub {#t_prime(x)
                        #warn "(eq d)\n";
                        return $t_prime_of_x_eqn_1and2->($_[0],    0    );
                    },
                    sub {#t_2prime(x)
                        #warn "(eq d)\n";
                        return $t_2prime_of_x_eqn_1and2->($_[0],    0    );
                    },
                    sub {
                         #warn "(eq d)\n";
                         return $t_of_y_eqn_1and2->($_[0],    1    );
                    },
                    sub {#t'(y)
                         #warn "(eq a)\n";
                        return $t_prime_of_y_eqn_1and2->($_[0],    1    );
                    },
                    sub {#t_2prime(y)
                        #warn "(eq d)\n";
                        return $t_2prime_of_y_eqn_1and2->($_[0],    1    );
                    }

                );
            }
            else { warn "out of bounds angle [$angle_mid]";}
        }

        push @XtoTLUT, [[@tx_eqns,@ty_eqns],
                        [$isReversed ? ($xb,$xa):($xa,$xb)],
                        [$isReversed ? ($tb,$ta):($ta,$tb)],
                        $isReversed,
                        [$ya > $yb ? ($yb,$ya):($ya,$yb)] # hack this in here quick to get F(y) working
                       ];

    }

    # sort in ascending t order
    @XtoTLUT = sort {$a->[2]->[$a->[3]?0:1] <=> $b->[2]->[$b->[3]?0:1]} @XtoTLUT;

    $self->{XtoTLUT} = \@XtoTLUT;

    return $self;
}

sub theta_of_x {
    my ($self,$x,$which) = @_;
    my $x_unshift = $x - $self->{cx};
    my $x_u = $x_unshift * cos(-$self->{phi_radians});
    my $y_u = $x_unshift * sin(-$self->{phi_radians});
    my $m = sin($pi/2 - $self->{phi_radians})/cos($pi/2 - $self->{phi_radians});
    my $y1_min_mx1 = $y_u-$m*$x_u;
    my $a =        $self->{ry}  / $y1_min_mx1;
    my $b = -($m * $self->{rx}) / $y1_min_mx1;
    my $phase_angle = atan2( $b , $a );
    my $theta  = $which
               ? asin(1/ sqrt($a**2 + $b**2)) - $phase_angle
               : asin(1/-sqrt($a**2 + $b**2)) - $phase_angle + $pi # + pi by educated trial and error
               ;
    if ($self->{sweep_flag}) {while ($theta < $self->{theta1}) {$theta += 2*$pi;}}
    else                     {while ($theta > $self->{theta1}) {$theta -= 2*$pi;}}
    return $theta;
}
sub theta_prime_of_x {
    my ($self,$x,$which) = @_;
    my $x_unshift = $x - $self->{cx};
    my $x_u = $x_unshift * cos(-$self->{phi_radians});
    my $y_u = $x_unshift * sin(-$self->{phi_radians});
    my $m = sin($pi/2 - $self->{phi_radians})/cos($pi/2 - $self->{phi_radians});
    my $z = (($self->{ry}**2+($m*$self->{rx})**2)/((1+$m**2)*sin(-$self->{phi_radians}+atan2(-$m,1))**2));
    my $theta_prime = ($which?1:-1) / sqrt($z - $x_unshift**2);
    return $theta_prime;
}
sub theta_2prime_of_x {
    my ($self,$x,$which) = @_;
    my $x_unshift = $x - $self->{cx};
    my $x_u = $x_unshift * cos(-$self->{phi_radians});
    my $y_u = $x_unshift * sin(-$self->{phi_radians});
    my $m = sin($pi/2 - $self->{phi_radians})/cos($pi/2 - $self->{phi_radians});
    my $z = (($self->{ry}**2+($m*$self->{rx})**2)/((1+$m**2)*sin(-$self->{phi_radians}+atan2(-$m,1))**2));
    my $theta_2prime = ($which?1:-1) * $x_unshift / (sqrt($z - $x_unshift**2)**3);
    return $theta_prime;
}

sub t_of_theta {
    my ($self, $arcTheta) = @_;
    my $num = ($arcTheta - $self->{theta1});
    return 0 if $num == 0; # avoid divide by zero if delta_theta == 0
    return $num / $self->{delta_theta};
}

sub t_prime_of_theta { return 1 / $self->{delta_theta}; }

sub theta_of_t {
    my ($self, $t) = @_;
    return $self->{theta1} + $self->{delta_theta} * $t;
}

sub theta_prime_of_t { return $self->{delta_theta}; }

sub t_prime_of_x {
    my ($self, $x, $which) = @_;
    # t_prime_of_theta * theta_prime_of_x
    my $ret = $self->theta_prime_of_x($x,$which) / $self->{delta_theta};
}

sub t_2prime_of_x {
    my ($self, $x, $which) = @_;
    # d_t/d_theta * d_theta/d_x
    my $ret = $self->theta_2prime_of_x($x,$which) / $self->{delta_theta}
            ;
}


sub solveXPrimeforThetaBig {return $_[0]->solveXPrimeforTheta($_[1])}
sub solveYPrimeforThetaBig {return $_[0]->solveYPrimeforTheta($_[1])}
sub evalXofThetaBig {return $_[0]->evalXofTheta($_[1])}
sub evalYofThetaBig {return $_[0]->evalYofTheta($_[1])}

sub initDangerRanges {
    my $self = shift;

    #redo these in all big - this might muck stuff up, so might delete this section later
    my @extremexs_is=(0,($self->solveXPrimeforThetaBig(Math::BigFloat->bzero())),1);
    my @extremexs = map {ref($_)?$self->evalXofThetaBig($_):$self->evalXofTheta($_)} @extremexs_is;
    my @extremexs_sorted = sort {$a<=>$b} @extremexs;

    my @extremeys_is=(0,($self->solveYPrimeforThetaBig(Math::BigFloat->bzero())),1);
    my @extremeys = map {ref($_)?$self->evalYofThetaBig($_):$self->evalYofTheta($_)} @extremeys_is;
    my @extremeys_sorted = sort {$a<=>$b} @extremeys;

    $self->{extremexs_is}=[(@extremexs_is)];
    $self->{extremeys_is}=[(@extremeys_is)];
    $self->{extremexs}=[(@extremexs)];
    $self->{extremeys}=[(@extremeys)];
    $self->{minx} = $extremexs_sorted[0];
    $self->{maxx} = $extremexs_sorted[$#extremexs];
    $self->{miny} = $extremeys_sorted[0];
    $self->{maxy} = $extremeys_sorted[$#extremeys];
    #end redo

    my $perlprecision = 0.00000000000001;
    my $maxdim=(sort {$b<=>$a} map {abs($_)} ($self->{maxx},$self->{maxy},$self->{minx},$self->{miny}))[0];
    $maxdim=~/^([0-9]+)\./;
    $self->{curveprecision}=$perlprecision * (10**(length($1))) ;
    my @xdangeris = sort {$a<=>$b} (
        ($self->solveXPrimeforTheta($perlprecision/$self->{curveprecision}))  , # near zero slope, positive
        ($self->solveXPrimeforTheta(0)),                  # zero slope
        ($self->solveXPrimeforTheta(-$perlprecision/$self->{curveprecision})) , # near zero slope, negative
        ($self->solveXPrimeforTheta($self->{curveprecision}/$perlprecision)) ,  # toward infinite slope, positive
        ($self->solveXPrimeforTheta(-$self->{curveprecision}/$perlprecision)) , # toward infinite slope, negative
    );
    my @ydangeris = sort {$a<=>$b} (
        ($self->solveYPrimeforTheta($perlprecision/$self->{curveprecision}))  ,
        ($self->solveYPrimeforTheta(0)) ,
        ($self->solveYPrimeforTheta(-$perlprecision/$self->{curveprecision})) ,
        ($self->solveYPrimeforTheta($self->{curveprecision}/$perlprecision)) ,
        ($self->solveYPrimeforTheta(-$self->{curveprecision}/$perlprecision)) ,
    );

    my @mxidangerranges;
    for (my $i=0;$i<@xdangeris;$i++) {
        my $p=$self->evalXPrimeofTheta($xdangeris[$i]);
        my $pp=$self->evalXDoublePrimeofTheta($xdangeris[$i]);
        if (($p < 0 && $pp < 0) || (($p > 0 || $p eq 0) && ($pp > 0 || $pp eq 0))) {
            #is end of range
            if ($i eq 0)          {push(@mxidangerranges,[0,(sort {$a<=>$b} (1,($xdangeris[$i] + $self->{curveprecision})))[0]]);}
            else                  {push(@mxidangerranges,[(sort {$b<=>$a} (0,($xdangeris[$i-1] - $self->{curveprecision})))[0],(sort {$a<=>$b} (1,($xdangeris[$i] + $self->{curveprecision})))[0]]);}
        }
        elsif ($i eq $#xdangeris) {push(@mxidangerranges,[(sort {$b<=>$a} (0,($xdangeris[$i] - $self->{curveprecision})))[0],1]);}
    }
    $self->{mxidangerranges} = \@mxidangerranges;

    $self->{xofidangerranges} = [];
    $self->{xofidangerranges}->[scalar(@{$self->{xofidangerranges}})]=[$self->{p1}->[0],$self->{p1}->[0]]; #zero-length "range" to
    $self->{xofidangerranges}->[scalar(@{$self->{xofidangerranges}})]=[$self->{p2}->[0],$self->{p2}->[0]]; #try to be more exact at endpoints

    foreach (@mxidangerranges) {
        $self->{xofidangerranges}->[scalar(@{$self->{xofidangerranges}})]=[$self->evalXofTheta($_->[0]),$self->evalXofThetaBig($_->[1])];
    }

    my @myidangerranges;
    $self->{yofidangerranges} = [];
    $self->{yofidangerranges}->[scalar(@{$self->{yofidangerranges}})]=[$self->{p1}->[1],$self->{p1}->[1]];
    $self->{yofidangerranges}->[scalar(@{$self->{yofidangerranges}})]=[$self->{p2}->[1],$self->{p2}->[1]]; #try to be more exact at endpoints
    for (my $i=0;$i<@ydangeris;$i++) {
        my $p=$self->evalYPrimeofTheta($ydangeris[$i]);
        my $pp=$self->evalYDoublePrimeofTheta($ydangeris[$i]);
        if (($p < 0 && $pp < 0) || (($p > 0 || $p eq 0) && ($pp > 0 || $pp eq 0))) {
            #is end of range
            if ($i == 0)          {push(@myidangerranges,[0,(sort {$a<=>$b} (1,$ydangeris[$i] + $self->{curveprecision}))[0]]);}
            else                  {push(@myidangerranges,[(sort {$b<=>$a} (0,$ydangeris[$i-1] - $self->{curveprecision}))[0],(sort {$a<=>$b} (1,$ydangeris[$i] + $self->{curveprecision}))[0]]);}
        }
        elsif ($i == $#ydangeris) {push(@myidangerranges,[(sort {$b<=>$a} (0,$ydangeris[$i] - $self->{curveprecision}))[0],1]);}
    }
    $self->{myidangerranges} = \@myidangerranges;
    foreach (@myidangerranges) {
        $self->{yofidangerranges}->[scalar(@{$self->{yofidangerranges}})]=[$self->evalYofThetaBig($_->[0]),$self->evalYofThetaBig($_->[1])];
    }
}
sub initBigs {
    my $self=shift;
    #nothing here for arc, right?
    $self->{bigInitted}=1;
}
sub precision {
    my $self=shift;
    if (defined($_[0])) {
        $self->{precision}=$_[0];
    }
    return $self->{precision};
}
sub getRange {
    my $self = shift;
    return ($self->{minx},$self->{miny},$self->{maxx},$self->{maxy});
}
sub inRange {
    my $self = shift;
    my $coords = shift;
    my $xok=0;
    my $yok=0;
    # should get rid of all these eval()s
    # but inRange is a potential source of bugginess, as experienced with the
    # bezier version, so don't mess unless you're ready to test
    if (defined($coords->[0]) && (eval($self->{minx}) < eval($coords->[0]) || eval($self->{minx}) eq eval($coords->[0])) && (eval($self->{maxx}) > eval($coords->[0]) || eval($self->{maxx}) eq eval($coords->[0]))) {$xok=1;}
    if (defined($coords->[1]) && (eval($self->{miny}) < eval($coords->[1]) || eval($self->{miny}) eq eval($coords->[1])) && (eval($self->{maxy}) > eval($coords->[1]) || eval($self->{maxy}) eq eval($coords->[1]))) {$yok=1;}
    return $xok,$yok;
}

sub isWithinThetaRange {
    my $self = shift;
    my $theta = shift;

    if ($self->{large_arc_flag}==0) {
        if ($self->{sweep_flag} == 0) {
            return (($theta < $self->{theta1} || $theta eq $self->{theta1}) && ($theta > $self->{theta2} || $theta eq $self->{theta2})) ? 1:0;
        }
        else {
            return (($theta > $self->{theta1} || $theta eq $self->{theta1}) && ($theta < $self->{theta2} || $theta eq $self->{theta2})) ? 1:0;
        }
    }
    else {
        if ($self->{sweep_flag} == 0) {
            return (($theta < $self->{theta1} || $theta eq $self->{theta1}) && ($theta > $self->{theta2} || $theta eq $self->{theta2})) ? 1:0;
        }
        else {
            return (($theta > $self->{theta1} || $theta eq $self->{theta1}) && ($theta < $self->{theta2} || $theta eq $self->{theta2})) ? 1:0;
        }
    }
}

sub f {
    my ($self, $x) = @_;
    my @ys = map  {$self->evalYofTheta($_->[0]->[0]->($x))}  # Y(t(x))
             grep {$x >= $_->[1]->[0] && $x < $_->[1]->[-1]} # x is within sub seg's x span
             @{$self->{XtoTLUT}};                            # in sorted t order
    return wantarray ? @ys : (scalar(@ys) ? $ys[0] : undef);
}

sub F {
    my ($self, $y) = @_;
    my @xs = map  {$self->evalXofTheta($_->[0]->[3]->($y))}  # X(t(y))
             grep {$y >= $_->[4]->[0] && $y < $_->[4]->[-1]} # y is within sub seg's y span
             @{$self->{XtoTLUT}};                            # in sorted t order
    return wantarray ? @xs : (scalar(@xs) ? $xs[0] : undef);
}

sub point {
    my $self = shift;
    my $theta = shift;

    my $arc_theta = $self->theta_of_t($theta);
    return if !$self->isWithinThetaRange($arc_theta);
    for (my $i=0;$i<scalar(@{$self->{mxidangerranges}});$i++) {
        if (($theta > $self->{mxidangerranges}->[$i]->[0] || $theta eq $self->{mxidangerranges}->[$i]->[0]) && ($theta < $self->{mxidangerranges}->[$i]->[1] || $theta eq $self->{mxidangerranges}->[$i]->[1]) && !ref($theta)) {
            $theta=Math::BigFloat->new($theta) if !ref($theta);
        }
    }

    for (my $i=0;$i<scalar(@{$self->{myidangerranges}});$i++) {
        if (($theta > $self->{myidangerranges}->[$i]->[0] || $theta eq $self->{myidangerranges}->[$i]->[0]) && ($theta < $self->{myidangerranges}->[$i]->[1] || $theta eq $self->{myidangerranges}->[$i]->[1]) && !ref($theta)) {
            $theta=Math::BigFloat->new($theta) if !ref($theta);
        }
    }

    #ellipse formulas derived from SVG spec
    # x = rx * cos(theta) * cos(phi) + ry * sin(theta) * -sin(phi) + Cx
    # y = rx * cos(theta) * sin(phi) + ry * sin(theta) *  cos(phi) + Cy

    # We had a BigFloat version of the math here, but it was commented out
    # and looked messy and unfinished, so deleted it.
    # You may need to reintroduce that at some point.
    # BigFloat math was done if $theta was a BigFloat.
    # Until that's needed again, downgrade any $theta that's a BigFloat.
    if (ref($theta)) {$theta=0 + sprintf("%.18f",$theta->bstr);}
    
    my $ct=cos($arc_theta);
    my $st=sin($arc_theta);
    return [$self->{rx} * $ct * cos($self->{phi_radians}) + $self->{ry} * $st * -sin($self->{phi_radians}) + $self->{cx},
            $self->{rx} * $ct * sin($self->{phi_radians}) + $self->{ry} * $st *  cos($self->{phi_radians}) + $self->{cy}];
}

sub point_offset() {
    my ($self, $t, $distance) = @_;
    my ($x, $y) = @{$self->point($t)};
    my $a = $self->angleNormal_byTheta($t);
    $x += cos($a) * $distance;
    $y += sin($a) * $distance;
    return [$x,$y];
}

sub X_offset {
    my ($self, $t, $distance) = @_;
    my $x = $self->evalXofTheta($t);
    my $a = $self->angleNormal_byTheta($t);
    $x += cos($a) * $distance;
    return $x;
}
sub Y_offset {
    my ($self, $t, $distance) = @_;
    my $y = $self->evalYofTheta($t);
    my $a = $self->angleNormal_byTheta($t);
    $y += sin($a) * $distance;
    return $y;
}

sub t_from_xoff {
    my ($self, $xoff, $distance, $t_bounds, $tprimeofxfunc) = @_;

    # solve for the arc theta where X_offset($distance) == $xoff
    # then to t_of_theta(theta) of course
    # you should pass in the subseg LUT entry stuff if you have it; otherwise, could cycle through all in LUT looking for solutions.
    # this is likely analytic. should work out on paper.
    # Tried to work it out. Seemed close. But couldn't finish.
    # So root finding approach for now.

    my $rfsub = sub {
        my $t = $_[0];
        my $x = $self->evalXofTheta($t);
        my $offset_x = $x + cos($self->angleNormal_byTheta($t)) * $distance;
        return $offset_x - $xoff;
    };
    $t_bounds = [$t_bounds->[1], $t_bounds->[0]] if $t_bounds->[0] > $t_bounds->[1];
    my ($t, $msg) = BrentsMethod($rfsub,$t_bounds,0.000001,undef,'finding t_from_xoff()');
    die "t_from_xoff() root find fail: $msg\n" if $msg;
    return $t;
}

sub t_from_xoff_prime {
    my ($self, $xoff, $distance, $t_bounds, $tprimeofxfunc) = @_;
}


sub evalYofTheta { # this "Theta" means "t" - fix
    my ($self, $t) = @_;
    my $theta = $self->theta_of_t($t);
    return if ! $self->isWithinThetaRange($theta);
    return $self->{rx} * (0 + sprintf("%.14f",cos($theta))) * sin($self->{phi_radians}) + $self->{ry} * (0 + sprintf("%.14f",sin($theta))) *  cos($self->{phi_radians}) + $self->{cy};
}

sub evalYofArcTheta {
    my ($self, $theta) = @_;
    return if !$self->isWithinThetaRange($theta);
    return $self->{rx} * (0 + sprintf("%.14f",cos($theta))) * sin($self->{phi_radians}) + $self->{ry} * (0 + sprintf("%.14f",sin($theta))) *  cos($self->{phi_radians}) + $self->{cy};
}

sub evalXofTheta { # this "Theta" means "t" - fix
    my ($self, $t) = @_;
    my $theta = $self->theta_of_t($t);
    return if !$self->isWithinThetaRange($theta);
    return $self->{rx} * (0 + sprintf("%.14f",cos($theta))) * cos($self->{phi_radians}) + $self->{ry} * (0 + sprintf("%.14f",sin($theta))) * -sin($self->{phi_radians}) + $self->{cx};
}

sub evalXofArcTheta {
    my ($self, $theta) = @_;
    return if !$self->isWithinThetaRange($theta);
    return $self->{rx} * (0 + sprintf("%.14f",cos($theta))) * cos($self->{phi_radians}) + $self->{ry} * (0 + sprintf("%.14f",sin($theta))) * -sin($self->{phi_radians}) + $self->{cx};
}

sub evalYPrimeofTheta { # this "Theta" means "t" - fix
    my ($self, $t) = @_;
    my $theta = $self->theta_of_t($t);
    return if !$self->isWithinThetaRange($theta);
    # Y'(arc_theta) * arc_theta'(t)
    # but arc_theta'(t) is just the constant delta_theta of the arc
    return (  $self->{rx} * (0 + sprintf("%.14f",-sin($theta))) * sin($self->{phi_radians})
            + $self->{ry} * (0 + sprintf("%.14f", cos($theta))) * cos($self->{phi_radians})
           )
           * $self->{delta_theta}
    ;
}
sub evalYPrimeofArcTheta {
    my ($self, $theta) = @_;
    return if !$self->isWithinThetaRange($theta);
    return (  $self->{rx} * (0 + sprintf("%.14f",-sin($theta))) * sin($self->{phi_radians})
            + $self->{ry} * (0 + sprintf("%.14f",cos($theta))) *  cos($self->{phi_radians})
           ) * ($self->{sweep_flag} ? 1:-1)
    ;
}
sub evalXPrimeofTheta { # this "Theta" means "t" - fix
    my ($self, $t) = @_;
    my $theta = $self->theta_of_t($t);
    return if !$self->isWithinThetaRange($theta);
    return (  $self->{rx} * (0 + sprintf("%.14f",-sin($theta))) *  cos($self->{phi_radians})
            + $self->{ry} * (0 + sprintf("%.14f", cos($theta))) * -sin($self->{phi_radians})
           )
           * $self->{delta_theta}
    ;
}
sub evalXPrimeofArcTheta {
    my ($self, $theta) = @_;
    return if !$self->isWithinThetaRange($theta);
    return (  $self->{rx} * (0 + sprintf("%.14f",-sin($theta))) *  cos($self->{phi_radians})
            + $self->{ry} * (0 + sprintf("%.14f", cos($theta))) * -sin($self->{phi_radians})
           ) * ($self->{sweep_flag} ? 1:-1)
    ;
}
sub evalYDoublePrimeofTheta { # this "Theta" means "t" - fix
    my ($self, $t) = @_;
    my $theta = $self->theta_of_t($t);
    return if !$self->isWithinThetaRange($theta);
    return (  $self->{rx} * (0 + sprintf("%.14f",-cos($theta))) * sin($self->{phi_radians})
            + $self->{ry} * (0 + sprintf("%.14f",-sin($theta))) * cos($self->{phi_radians})
           )
           * $self->{delta_theta}
    ;
}
sub evalXDoublePrimeofTheta { # this "Theta" means "t" - fix
    my ($self, $t) = @_;
    my $theta = $self->theta_of_t($t);
    return if !$self->isWithinThetaRange($theta);
    return (  $self->{rx} * (0 + sprintf("%.14f",-cos($theta))) *  cos($self->{phi_radians})
            + $self->{ry} * (0 + sprintf("%.14f",-sin($theta))) * -sin($self->{phi_radians})
           )
           * $self->{delta_theta}
    ;
}

# 6/16/2018 - WE HAVE NEW VERSIONS OF THESE ABOVE - 
# t(x) and t(y) math finally worked out, 
# along with breaking arc into monotonic sections -
# so these need to be obsolete,
# or they could go through the LUT and collect all ts for the Y (or X)
sub solveYforTheta { # this "Theta" means "t", but inside this function "$t" is used for arc_theta! - rename stuff to clarify
    my $self = shift;
    my $y = shift;
    my @xs = $self->F($y);
    my @ts;
    for (my $i=0;$i<@xs;$i++) {

        my $unrot = [$xs[$i] - $self->{cx}, $y - $self->{cy}];
        $unrot->[0] = ($unrot->[0]*cos(-$self->{phi_radians}) - $unrot->[1]*sin(-$self->{phi_radians}));
        $unrot->[1] = ($unrot->[0]*sin(-$self->{phi_radians}) + $unrot->[1]*cos(-$self->{phi_radians}));

        #now corresponding ellipse formulas are
        #x=$self->{rx} * cos(theta);
        #y=$self->{ry} * sin(theta);
        #so that's solvable for theta, right?
        #y=$self->{ry} * sin(theta);
        #theta=asin(y/ry);

        my $t=asin($unrot->[1]/$self->{ry});

        $t*=-1 if $self->{sweep_flag};

        my $other_t=$t + (($t>0)?-1:1) * $pi;
        my $other_t2=$t + (($t>0)?-1:1) * 2*$pi;
        my $other_t3=$other_t + (($other_t>0)?-1:1) * 2*$pi;
        push(@ts,$t,$other_t,$other_t2,$other_t3);
    }
    return map {$self->t_of_theta($_)} grep {$self->isWithinThetaRange($_) && abs($self->evalYofTheta($_) - $y)<0.0000001} @ts;
}
sub solveXforTheta { # this "Theta" means "t", but inside this function "$t" is used for arc_theta! - rename stuff to clarify
    my $self = shift;
    my $x = shift;
    my @ys = $self->f($x);
    my @ts;
    for (my $i=0;$i<@ys;$i++) {

        my $unrot = [$x - $self->{cx}, $ys[$i] - $self->{cy}];
        $unrot->[0] = ($unrot->[0]*cos(-$self->{phi_radians}) - $unrot->[1]*sin(-$self->{phi_radians}));
        $unrot->[1] = ($unrot->[0]*sin(-$self->{phi_radians}) + $unrot->[1]*cos(-$self->{phi_radians}));

        #now corresponding ellipse formulas are
        #x=$self->{rx} * cos(theta);
        #y=$self->{ry} * sin(theta);
        #so that's solvable for theta, right?
        #x=$self->{rx} * cos(theta);
        #theta=acos(x/rx);
        #theta=pi/2 - asin(x/rx)

        my $t=($pi/2) - asin($unrot->[0]/$self->{rx});

        $t*=-1 if $self->{sweep_flag};

# TODO
# MAKE SURE ALL THIS IS SIMILAR IN solveYforTheta ABOVE AND IN JAVASCRIPT VERSION

        my $other_t  =       $t + (($t      >0) ?-1:1) *   $pi;
        my $other_t2 =       $t + (($t      >0) ?-1:1) * 2*$pi;
        my $other_t3 = $other_t + (($other_t>0) ?-1:1) * 2*$pi;

        push(@ts,$t,$other_t,$other_t2,$other_t3);
    }

    return map  {$self->t_of_theta($_)}
           grep {   $self->isWithinThetaRange($_)
                 && abs($self->evalXofArcTheta($_) - $x)<0.0000001
           } @ts;
}
sub solveYPrimeforTheta { # this "Theta" means "t", but inside this function "$t" is used for arc_theta! - rename stuff to clarify
    my $self = shift;
    my $yp = shift;
    #first of all, that slope would be this slope, if ellipse were "unrotated":
    #yp_unrot = tan(atan(yp) - self->phi)
    my $yp_unrot = sin(atan2($yp,1) - $self->{phi_radians})/cos(atan2($yp,1) - $self->{phi_radians});
    #the unrotated, centered ellipse formulas are
    #x=rx * cos(theta)
    #y=ry * sin(theta)
    #derivatives of those, with respect to theta:
    #x'=rx * -sin(theta)
    #y'=ry * cos(theta)
    #and dy/dx is
    #(dy/dt) / (dx/dt)
    #so you want to solve this for theta:
    #yp_unrot = (ry * cos(theta)) / (rx * -sin(theta))
    #which I hope comes out this simply:
    #theta = atan((ry/-rx) * 1/yp_unrot)
    my $t = atan2($self->{ry},-$self->{rx} * $yp_unrot);
    my $other_t=$t + (($t>0 || $t eq 0)?-1:1) * $pi;
    my $other_t2=$t + (($t>0 || $t eq 0)?-1:1) * 2*$pi;
    my $other_t3=$other_t + (($other_t>0 || $other_t eq 0)?-1:1) * 2*$pi;
    return map {$self->t_of_theta($_)} grep {$self->isWithinThetaRange($_)} ($t,$other_t,$other_t2,$other_t3);
}
sub solveXPrimeforTheta { # this "Theta" means "t", but inside this function "$t" is used for arc_theta! - rename stuff to clarify
    my $self = shift;
    my $xp = shift;
    #first of all, that slope would be this slope, if ellipse were "unrotated":
    #xp_unrot = tan(atan(xp) - self->phi)
    my $xp_unrot = sin(atan2($xp,1) - $self->{phi_radians})/cos(atan2($xp,1) - $self->{phi_radians});
    #the unrotated, centered ellipse formulas are
    #x=rx * cos(theta)
    #y=ry * sin(theta)
    #derivatives of those, with respect to theta:
    #x'=rx * -sin(theta)
    #y'=ry * cos(theta)
    #and dx/dy is
    #(dx/dt) / (dy/dt)
    #so you want to solve this for theta:
    #xp_unrot = (rx * -sin(theta)) / (ry * cos(theta))
    #which I hope comes out this simply:
    #theta = atan((ry*xp_unrot)/-rx/)
    my $t = atan2($self->{ry}*$xp_unrot,-$self->{rx});
    my $other_t=$t + (($t>0 || $t eq 0)?-1:1) * $pi;
    my $other_t2=$t + (($t>0 || $t eq 0)?-1:1) * 2*$pi;
    my $other_t3=$other_t + (($other_t>0 || $other_t eq 0)?-1:1) * 2*$pi;
    return map {$self->t_of_theta($_)} grep {$self->isWithinThetaRange($_)} ($t,$other_t,$other_t2,$other_t3);
}

sub getLength {
    my ($self, $res, $start_t, $end_t) = @_;

    if (!defined($res)) {$res=1000;}
    if (!defined($start_t)) {$start_t=0;}
    if (!defined($end_t)) {$end_t=1;}



    # if the two radii are equal, it's a circular arc
    # (which is usually how I use this ellipse stuff)
    # and we have 10th grade math for that
    if ($self->{rx} eq $self->{ry}) {
        return abs($self->{rx} * $self->{delta_theta} * ($end_t - $start_t));
    }


    my $sum = 0;
    my $point1 = $self->point($start_t);
    my $point2;
    my $t_span = ($end_t - $start_t);
    my $t_inc = $t_span/$res;
    for (my $i=0;$i<$t_span;$i+=$t_inc) {
        $point2 = $point1;
        $point1 = $self->point($start_t + $i);
        $sum += sqrt(($point2->[0]-$point1->[0])**2 + ($point2->[1]-$point1->[1])**2);
    }
    return $sum;
    }

sub getFeet {
    my $self=shift;
    my $x=shift;
    my $y=shift;
    # here I used "i" to mean "t" - the 0 to 1 parameter
    #for each interval between critical points - critical due to features of x(i) and y(i). - hueristic and then root find to get any 90 degree intersections
    my @feet=();
    my %dupsieve_i_set;
    my @extreme_i_set= grep {!$dupsieve_i_set{$_}++} sort {$a<=>$b} (@{$self->{extremeys_is}},@{$self->{extremexs_is}});
    for (my $i=@extreme_i_set - 1;$i>0;$i--) {
        splice(@extreme_i_set,$i,0,(($extreme_i_set[$i]+$extreme_i_set[$i - 1])/2))
    }

    for (my $i=1;$i<@extreme_i_set;$i++) {
        next if ($extreme_i_set[$i - 1] eq $extreme_i_set[$i]);
        my $boundl=$extreme_i_set[$i - 1];
        my $boundr=$extreme_i_set[$i];
        my $find90 = sub {
            #dotproduct equals zero for two perpendicular vectors
            my $rayvec=[$self->evalXofTheta($_[0]) - $x,$self->evalYofTheta($_[0]) - $y];
            my $tanvec = [1,0 + $self->slopeTangent_byTheta($_[0])];
            my $ret = $rayvec->[0] * $tanvec->[0] + $rayvec->[1] * $tanvec->[1];
            return $ret;
        };
        my $find90_2 = sub {
            #dotproduct equals zero for two perpendicular vectors
            my $rayvec=[$self->evalYofTheta($_[0]) - $y,-1 * ($self->evalXofTheta($_[0]) - $x)];
            my $tanvec = [1,0 + $self->slopeNormal_byTheta($_[0])];
            my $ret = $rayvec->[0] * $tanvec->[0] + $rayvec->[1] * $tanvec->[1];
            return $ret;
        };
        my $subtouse=$find90;
        if (abs(cos($self->angleTangent_byTheta($boundl))) < 0.001 || abs(cos($self->angleTangent_byTheta($boundr))) < 0.001) {
            #print "using find90_2\n";
            $subtouse=$find90_2;
        }
        else {
            #print "using find90\n";
        }
        my ($foot_i,$msg);
        ($foot_i,$msg)=BrentsMethod($subtouse,[$boundl,$boundr],$self->{precision},undef,'trying to find feet on Elliptical Arc in MPath with precision:'.($self->{precision}).'');
        if (!defined($msg)) {
            my $retpoint = $self->point($foot_i);
            push(@{$retpoint},$foot_i);
            push(@feet,$retpoint);
        }
    }
    return @feet;
}

sub angleTangent {
    my ($self, $x, $y, $t) = @_;

    my @ts;
    my @ret;

    if (defined($x)) {
        push @ts, map $_->[0]->[0]->($x),
                  grep {   ($x > $_->[1]->[0]  || $x eq $_->[1]->[0])
                        && ($x < $_->[1]->[-1] || $x eq $_->[1]->[-1])
                  } @{$self->{XtoTLUT}};
    }
    elsif (defined($y)) {
        push @ts, map $_->[0]->[3]->($y),
                  grep {   ($y > $_->[4]->[0]  || $y eq $_->[4]->[0])
                        && ($y < $_->[4]->[-1] || $y eq $_->[4]->[-1])
                   } @{$self->{XtoTLUT}};
    }
    elsif (defined($t)) {
        push @ts, $t;
    }

    foreach my $this_t (@ts) {
        my $arctheta = $self->theta_of_t($this_t);
        push @ret, atan2( $self->evalYPrimeofArcTheta($arctheta),
                          $self->evalXPrimeofArcTheta($arctheta)
                        );
    }

    return wantarray ? @ret : $ret[0];
}
sub slopeTangent {
    my @ats = $_[0]->angleTangent($_[1],$_[2],$_[3]);
    my @ret;
    for (my $i=0;$i<@ats;$i++) {
        push(@ret, sin($ats[$i])/cos($ats[$i]));
    } 
    return wantarray ? @ret : $ret[0];
}
sub slopeNormal  {
    my @ats = $_[0]->angleTangent($_[1],$_[2],$_[3]);
    my @ret;
    for (my $i=0;$i<@ats;$i++) {
        push(@ret, -cos($ats[$i])/sin($ats[$i]));
    }
    return wantarray ? @ret : $ret[0];
}
sub angleNormal  {
    my @ret = map {$_ + $pi/2} $_[0]->angleTangent($_[1],$_[2],$_[3]);
    @ret = map {angle_reduce($_)} @ret;
    return wantarray ? @ret : $ret[0];
}
sub angleTangent_byTheta { # this "Theta" means "t" - rename stuff to clarify - this and below, and etc. everywhere
    my ($self,$t) = @_;
    return $self->angleTangent(undef,undef,$t);
    }
sub angleNormal_byTheta {
    my ($self,$t) = @_;
    return $self->angleNormal(undef,undef,$t);
    }
sub slopeTangent_byTheta {
    my ($self,$t) = @_;
    return $self->slopeTangent(undef,undef,$t);
}
sub slopeNormal_byTheta {
    my ($self,$t) = @_;
    return $self->slopeNormal(undef,undef,$t);
}

sub t_t_prime_of_y {
    my ($self, $y, $t)=@_;
    $y //= $self->evalYofTheta($t);
    my @xspans = grep {($y > $_->[4]->[0] || $y eq $_->[4]->[0]) && ($y < $_->[4]->[-1] || $y eq $_->[4]->[-1])} @{$self->{XtoTLUT}};
    if (defined $t) {
        @xspans = grep {
            my ($lt,$ht)=$_->[2]->[0] < $_->[2]->[-1]?($_->[2]->[0] < $_->[2]->[-1]):($_->[2]->[-1] < $_->[2]->[0]);
            $t >= $lt && $t <= $ht;
        } @xspans;
    }
    my @ret = map {[$_->[0]->[3]->($y),$_->[0]->[4]->($y),$_->[3]]} @xspans;
    return wantarray ? @ret : $ret[0];
}

sub t_t_prime_t_2prime_of_y {
    my ($self, $y, $t)=@_;
    $y //= $self->evalYofTheta($t);
    my @xspans = grep {($y > $_->[4]->[0] || $y eq $_->[4]->[0]) && ($y < $_->[4]->[-1] || $y eq $_->[4]->[-1])} @{$self->{XtoTLUT}};
    if (defined $t) {
        @xspans = grep {
            my ($lt,$ht)=$_->[2]->[0] < $_->[2]->[-1]?($_->[2]->[0] < $_->[2]->[-1]):($_->[2]->[-1] < $_->[2]->[0]);
            $t >= $lt && $t <= $ht;
        } @xspans;
    }
    my @ret = map {[$_->[0]->[3]->($y),$_->[0]->[4]->($y),$_->[0]->[5]->($y),$_->[3]]} @xspans;
    return wantarray ? @ret : $ret[0];
}

sub F_prime {
    my ($self,$y,$t) = @_;

    my @t_t_prime_of_y = $self->t_t_prime_of_y($y,$t);
    $y //= $self->evalYofTheta($t);
    my @ret;
    foreach my $t_t_prime_of_y (@t_t_prime_of_y) {
        my ($t_of_y,$t_prime_of_y) = @$t_t_prime_of_y;
        # F'(y) = X'(t) * t'(y)
        my $x_prime_of_y = $self->evalXPrimeofTheta($t_of_y) * $t_prime_of_y;
        push @ret, $x_prime_of_y;
    }
    if (@ret>0) {
        return wantarray ? @ret : $ret[0];
    }
    else {
        return;
    }
}

sub F_2prime {
    my ($self,$y,$t) = @_;

    my @t_t_prime_t_2prime_of_y = $self->t_t_prime_t_2prime_of_y($y);
    $y //= $self->evalYofTheta($t);
    my @ret;
    foreach my $t_t_prime_t_2prime_of_y (@t_t_prime_t_2prime_of_y) {
        my ($t,$t_prime,$t_2prime) = @$t_t_prime_t_2prime_of_y;
        # F''(y) = [X'(t) * t'(y)]'
        #        = ( ( X''(t) * t'(y) ) * t'(y) + X'(t) * t''(y))
        #        = ( X''(t) * t'(y)^2 + X'(t) * t''(y))
        my $x_2prime_of_y = $self->evalXDoublePrimeofTheta($t) * $t_prime**2
                          + $self->evalXPrimeofTheta($t) * $t_2prime;
        push @ret, $x_2prime_of_y;
    }
    if (@ret>0) {
        return wantarray ? @ret : $ret[0];
    }
    else {
        return;
    }

}


sub secondDerivative { #think we might call it something else. follow what you did/do in cubic Bezier stuff
    my $self = shift;
    die "Need to work out second derivative with respect to x for elliptical arc when rx != ry";
}

sub isWithinSweep {
    my $self = shift;
    my $p = shift;
    my $leg1 = shift;
    my $leg2 = shift;
    my $leftness_1 = _howleft([[$leg1->[0]->[0],$leg1->[0]->[1]],[$leg1->[1]->[0],$leg1->[1]->[1]]],[$p->[0],$p->[1]]);
    my $leftness_2 = _howleft([[$leg2->[0]->[0],$leg2->[0]->[1]],[$leg2->[1]->[0],$leg2->[1]->[1]]],[$p->[0],$p->[1]]);
    if (($leftness_1 < 0 && $leftness_2 > 0) || $leftness_1==0 || $leftness_2==0) {return 1;}
    else {return 0;}

}
sub _howleft { #just copied from CAD::Calc
    my ($line, $pt) = @_;
    my $isleft = ($line->[1]->[0] - $line->[0]->[0]) *
                        ($pt->[1] - $line->[0]->[1]) -
                 ($line->[1]->[1] - $line->[0]->[1]) *
                         ($pt->[0] - $line->[0]->[0]);
    return($isleft);
}
sub asin {
    #Based on Wolfram MathWorld
    #http://mathworld.wolfram.com/InverseSine.html
    #but it's supposed to use complex numbers, hmmm
    if ($_[0] eq -1)    {return $pi/2 * -1} #eq 2
    elsif ($_[0] eq 1)  {return $pi/2     } #eq 4
    elsif ($_[0] eq 0)  {return 0;             } #eq 3
    else {                                       #eq 15
        if (abs($_[0])>1) {warn("giving something bigger than +/-1 to your asin:",$_[0],"\n");}
        return atan2($_[0],sqrt(1-$_[0]**2));
    }
}
sub angle_reduce { # Copied from CAD::Calc
    my $ang = shift;
    my $angbef = $ang;
    while($ang > $pi) {
        $ang -= 2*$pi;
    }
    while($ang <= -$pi) {
        $ang += 2*$pi;
    }
    return($ang);
}
sub dimensionalStepFromTheta {

    #same code as for bezier
    #and maybe even line segment if you've been to lazy to do a simpler version for the simple line case
    #it's been fast so far
    #didn't test arc case after pasting this, but line segment case worked with no mods, so this
    #has a good chance of working too.

    my ($self,$dim,$theta,$direction) = @_;
    if ( ! defined ($direction) ) { $direction = 1; }

    my $pt_last = $self->point($theta);

    my $findnexttheta = sub {
        my $ret;
        if (!ref($_[0])) {
            my $pt_new  = $self->point($_[0]);
            $ret = $dim - CORE::sqrt(($pt_new->[0] - $pt_last->[0])**2 + ($pt_new->[1] - $pt_last->[1])**2);
            #print "$ret = $dim - CORE::sqrt(($pt_new->[0] - $pt_last->[0])**2 + ($pt_new->[1] - $pt_last->[1])**2) = $ret\n";
        }
        else {
            #warn "I don't think you want to be here - not sure if this mess is debugged.\n";
            my $pt_new  = $self->point($_[0]);
            #print "using BigFloat - \n";
            my $dx=(ref($pt_new->[0]))?$pt_new->[0]->copy()->bsub($pt_last->[0]):$pt_new->[0] - $pt_last->[0];
            my $dxsqrd=(ref($dx))?$dx->bpow(2):$dx**2;
            my $dy=(ref($pt_new->[1]))?$pt_new->[1]->copy()->bsub($pt_last->[1]):$pt_new->[1] - $pt_last->[1];
            my $dysqrd=(ref($dy))?$dy->bpow(2):$dy**2;
            my $distsqrd=(ref($dysqrd))?$dysqrd->copy()->badd($dxsqrd):$dysqrd + $dxsqrd;
            my $dist = (ref($distsqrd))?bigsqrt($distsqrd):sqrt($distsqrd);
            $ret = (ref($dim) ? $dim:Math::BigFloat->new(''.$dim)) - $dist;
        }
        return $ret;
    };

    my $newtheta;
    my $er;
    ($newtheta,$er) = FalsePosition($findnexttheta,($direction ? [$theta,1]:[0,$theta]),$self->{precision},($direction ? ($theta + (1-$theta)/2):($theta/2)),'dimensionalStepFromTheta for Elliptical Arc segment');
    if (defined($er)) {
        #probably just reached the end
        if (abs(&{$findnexttheta}(($direction ? 1:0))) < $dim) {
            $newtheta=($direction ? 1:0);
        }
        #otherwise the error might be real
        else {
            warn "dimstep er: $er";
        }
    }
    return ($self->point($newtheta),$newtheta);
}

} # end package
1;
