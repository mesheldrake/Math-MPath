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
    $self->{phi} = ((1000000*$self->{phi}) % (360000000))/1000000; # angle starts as degrees, mod 360
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

    #Make all these theta mods rational, in accordance with what's supposed to happen, overcoming any angle flattenning that comes from atan sqrt usage above

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
    $self->{f1} = _rotate2d([0,0],$self->{f1},$self->{phi_radians});
    $self->{f1}->[0]+=$self->{cx};
    $self->{f1}->[1]+=$self->{cy};
    $self->{f2} = _rotate2d([0,0],$self->{f2},$self->{phi_radians});
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

    $self->{length}=getLength($self,1000,0,1);

    return $self;
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
sub arcThetaToNormalizedTheta {
    my $self=shift;
    my $arcTheta=shift;
    my $num = ($arcTheta - $self->{theta1});
    return 0 if $num == 0; # avoids divide by zero error when below comes out as 0/0
    return $num / $self->{delta_theta};
}
sub normalizedThetaToArcTheta {
    my $self = shift;
    my $normalizedTheta = shift;
    return $self->{theta1} + $self->{delta_theta} * $normalizedTheta;
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
    my $self = shift;
    my $x = shift;
    my @intersections;
    $x -= $self->{cx};
    if ($x < -$self->{rx} || $x > $self->{rx}) {return;}

    for (my $i=0;$i<scalar(@{$self->{fxdangerranges}});$i++) {
        if (   ($x > $self->{'fxdangerranges'}->[$i]->[0] || $x eq $self->{'fxdangerranges'}->[$i]->[0])
            && ($x < $self->{'fxdangerranges'}->[$i]->[1] || $x eq $self->{'fxdangerranges'}->[$i]->[1])
            && !ref($x)
           ) {
            #print "DANGER ZONE FOR f(x) (arc)! x:$x\n";
            $x=Math::BigFloat->new($x) if !ref($x);
        }
    }

    my $rot_line_slope=sin($pi/2 - $self->{phi_radians})/cos($pi/2 - $self->{phi_radians});

    if (abs($rot_line_slope) > 1.0 * 10**6 || $rot_line_slope eq 'inf' || $rot_line_slope eq '-inf') {
        my $y=sqrt($self->{ry}**2 * (1 - ($x**2)/($self->{rx}**2)));#vertical line. use ellipse formula to get the +/- y vals
        push(@intersections,[$x,$y],[$x,-$y]);
    }
    else {
        my $rot_line = _rotate2d([0,0],[$x,0],-$self->{phi_radians}); # point on a vertical x=C line getting tilted into ellipse frame, where the line will have a slope equal to -tan(phi)
        
        my $a = (($rot_line_slope)**2/$self->{ry}**2) + 1/$self->{rx}**2;
        my $b = ( 2 * ($rot_line_slope) * ($rot_line->[1] - ($rot_line_slope)*$rot_line->[0]))/$self->{ry}**2;
        my $c =(($rot_line->[1] - ($rot_line_slope)*$rot_line->[0])**2 / $self->{ry}**2 ) - 1;
        my @xs = &quadraticformula($a,$b,$c,1);
        for (my $i=0;$i<@xs;$i++) {
            my $y=$rot_line_slope * $xs[$i] + ($rot_line->[1] - $rot_line_slope * $rot_line->[0]); #line formula
            push @intersections, [ $xs[$i], $y ];
        }
    }

    for (my $i=0;$i<@intersections;$i++) {
        my $h=sqrt($intersections[$i]->[0]**2 + $intersections[$i]->[1]**2);
        $intersections[$i] = _rotate2d([0,0],$intersections[$i],$self->{phi_radians});
        $intersections[$i]->[0]+=$self->{cx};
        $intersections[$i]->[1]+=$self->{cy};
    }

    #Now check to see of those intersections are within bounds - within sweep

    my $leg1;
    my $leg2;
    if ($self->{large_arc_flag}==0) {
        if ($self->{sweep_flag} == 0) {
            $leg1=[[$self->{cx},$self->{cy}],[$self->{p1}->[0],$self->{p1}->[1]]];
            $leg2=[[$self->{cx},$self->{cy}],[$self->{p2}->[0],$self->{p2}->[1]]];
        }
        else {
            $leg1=[[$self->{cx},$self->{cy}],[$self->{p2}->[0],$self->{p2}->[1]]];
            $leg2=[[$self->{cx},$self->{cy}],[$self->{p1}->[0],$self->{p1}->[1]]];
        }
    }
    else {
        if ($self->{sweep_flag} == 0) {
            $leg1=[[$self->{cx},$self->{cy}],[$self->{p2}->[0],$self->{p2}->[1]]];
            $leg2=[[$self->{cx},$self->{cy}],[$self->{p1}->[0],$self->{p1}->[1]]];
        }
        else {
            $leg1=[[$self->{cx},$self->{cy}],[$self->{p1}->[0],$self->{p1}->[1]]];
            $leg2=[[$self->{cx},$self->{cy}],[$self->{p2}->[0],$self->{p2}->[1]]];
        }
    }
    @intersections = grep { ($self->{large_arc_flag} && !$self->isWithinSweep($_,$leg1,$leg2)) || (!$self->{large_arc_flag} && $self->isWithinSweep($_,$leg1,$leg2)) } @intersections;
    return wantarray ? (map {$_->[1]} @intersections) : $intersections[0]->[1];
}

sub F {
    my $self = shift;
    my $y = shift;
    my @intersections;
    $y -= $self->{cy};
    if ($y < -$self->{ry} || $y > $self->{ry}) {return;}

    for (my $i=0;$i<scalar(@{$self->{Fydangerranges}});$i++) {
        if (   ($y > $self->{Fydangerranges}->[$i]->[0] || $y eq $self->{Fydangerranges}->[$i]->[0]) 
            && ($y < $self->{Fydangerranges}->[$i]->[1] || $y eq $self->{Fydangerranges}->[$i]->[1])
            && !ref($y)
           ) {
            #print "DANGER ZONE FOR F(y) (arc) ! y:$y ";
            if ($Math::MPath::enableCarefulFofy) {
                $y=Math::BigFloat->new(''.$y) if !ref($y);
            }
            else {
                #print " BUT ignoring because \$Math::MPath::enableCarefulFofy = $Math::MPath::enableCarefulFofy";
            }
            #print "\n";
        }
    }

    my $rot_line_slope = sin(-$self->{phi_radians})/cos(-$self->{phi_radians});

    if (abs($rot_line_slope) > 1.0 * 10**6 || $rot_line_slope eq 'inf' || $rot_line_slope eq '-inf' || abs($rot_line_slope) < 1.0 * 10**-10) {
        my $x = sqrt($self->{rx}**2 * (1 - ($y**2)/($self->{ry}**2)));#vertical line. use ellipse formula to get the +/- y vals
        push(@intersections,[$x,$y],[-$x,$y]);
    }
    elsif (abs($rot_line_slope) < 1.0 * 10**-10) {
        my $x = sqrt($self->{rx}**2 * (1 - ($y**2)/($self->{ry}**2)));#vertical line. use ellipse formula to get the +/- y vals
        push(@intersections,[$x,$y],[-$x,$y]);
    }
    else {
        my $rot_line = _rotate2d([0,0],[0,$y],-$self->{phi_radians}); # point on a vertical x=C line getting tilted into ellipse frame, where the line will have a slope equal to -tan(phi)
        my $a = (1/$self->{ry}**2) + 1/($self->{rx}**2 * $rot_line_slope**2);
        my $b = (2*($rot_line->[0] - ($rot_line->[1]/$rot_line_slope)))/($self->{rx}**2 * $rot_line_slope);
        my $c = (($rot_line->[0] - ($rot_line->[1]/$rot_line_slope))**2 / $self->{rx}**2) - 1;
        my @ys = &quadraticformula($a,$b,$c,1);
        for (my $i=0;$i<@ys;$i++) {
            my $x=(($ys[$i] - $rot_line->[1])/$rot_line_slope) + $rot_line->[0]; #line formula
            push @intersections, [ $x, $ys[$i] ];
        }
    }

    for (my $i=0;$i<@intersections;$i++) {
        my $h=sqrt($intersections[$i]->[0]**2 + $intersections[$i]->[1]**2);
        $intersections[$i] = _rotate2d([0,0],$intersections[$i],$self->{phi_radians});
        $intersections[$i]->[0]+=$self->{cx};
        $intersections[$i]->[1]+=$self->{cy};
    }

    #Now check to see of those intersections are within bounds - within sweep

    my $leg1;
    my $leg2;
    if ($self->{large_arc_flag}==0) {
        if ($self->{sweep_flag} == 0) {
            $leg1=[[$self->{cx},$self->{cy}],[$self->{p1}->[0],$self->{p1}->[1]]];
            $leg2=[[$self->{cx},$self->{cy}],[$self->{p2}->[0],$self->{p2}->[1]]];
        }
        else {
            $leg1=[[$self->{cx},$self->{cy}],[$self->{p2}->[0],$self->{p2}->[1]]];
            $leg2=[[$self->{cx},$self->{cy}],[$self->{p1}->[0],$self->{p1}->[1]]];
        }
    }
    else {
        if ($self->{sweep_flag} == 0) {
            $leg1=[[$self->{cx},$self->{cy}],[$self->{p2}->[0],$self->{p2}->[1]]];
            $leg2=[[$self->{cx},$self->{cy}],[$self->{p1}->[0],$self->{p1}->[1]]];
        }
        else {
            $leg1=[[$self->{cx},$self->{cy}],[$self->{p1}->[0],$self->{p1}->[1]]];
            $leg2=[[$self->{cx},$self->{cy}],[$self->{p2}->[0],$self->{p2}->[1]]];
        }
    }
    @intersections = grep { ($self->{large_arc_flag} && !$self->isWithinSweep($_,$leg1,$leg2)) || (!$self->{large_arc_flag} && $self->isWithinSweep($_,$leg1,$leg2)) } @intersections;
    return wantarray ? (map {$_->[0]} @intersections) : $intersections[0]->[0];
}
sub point {
    my $self = shift;
    my $theta = shift;

    my $arc_theta = $self->normalizedThetaToArcTheta($theta);
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

sub evalYofTheta {
    my $self = shift;
    #print "evalYofTheta     theta: $_[0]\n";
    my $theta = $self->normalizedThetaToArcTheta(shift);
    #print "evalYofTheta arc theta: $theta\n";
    #if (!$self->isWithinThetaRange($theta)) {print "evalYofTheta (OUT OF RANGE)\n";}
    return if !$self->isWithinThetaRange($theta);
    #my $ret=$self->{rx} * eval(sprintf("%.14f",cos($theta))) * sin($self->{phi_radians}) + $self->{ry} * eval(sprintf("%.14f",sin($theta))) *  cos($self->{phi_radians}) + $self->{cy};
    #print "evalYofTheta $ret=$self->{rx} * cos($theta) * sin($self->{phi_radians}) + $self->{ry} * sin($theta) *  cos($self->{phi_radians}) + $self->{cy};\n";
    return $self->{rx} * (0 + sprintf("%.14f",cos($theta))) * sin($self->{phi_radians}) + $self->{ry} * (0 + sprintf("%.14f",sin($theta))) *  cos($self->{phi_radians}) + $self->{cy};
}
# same as above, but provided theta is alread the arc theta and not the normalized theta
# so should probably have above jjust do the theat conversion then call this
sub evalYofArcTheta {
    my $self = shift;
    my $theta = shift;
    return if !$self->isWithinThetaRange($theta);
    return $self->{rx} * (0 + sprintf("%.14f",cos($theta))) * sin($self->{phi_radians}) + $self->{ry} * (0 + sprintf("%.14f",sin($theta))) *  cos($self->{phi_radians}) + $self->{cy};
}
sub evalXofTheta {
    my $self = shift;
    my $theta = $self->normalizedThetaToArcTheta(shift);
    return if !$self->isWithinThetaRange($theta);
    return $self->{rx} * (0 + sprintf("%.14f",cos($theta))) * cos($self->{phi_radians}) + $self->{ry} * (0 + sprintf("%.14f",sin($theta))) * -sin($self->{phi_radians}) + $self->{cx};
}
# same as above, but provided theta is already the arc theta and not the normalized theta
# so should probably have above just do the theata conversion then call this
sub evalXofArcTheta {
    my $self = shift;
    my $theta = shift;
    return if !$self->isWithinThetaRange($theta);
    return $self->{rx} * (0 + sprintf("%.14f",cos($theta))) * cos($self->{phi_radians}) + $self->{ry} * (0 + sprintf("%.14f",sin($theta))) * -sin($self->{phi_radians}) + $self->{cx};
}
sub evalYPrimeofTheta {
    my $self = shift;
    my $theta = $self->normalizedThetaToArcTheta(shift);
    return if !$self->isWithinThetaRange($theta);
    return $self->{rx} * (0 + sprintf("%.14f",-sin($theta))) * sin($self->{phi_radians}) + $self->{ry} * (0 + sprintf("%.14f",cos($theta))) *  cos($self->{phi_radians});
}
sub evalXPrimeofTheta {
    my $self = shift;
    my $theta = $self->normalizedThetaToArcTheta(shift);
    return if !$self->isWithinThetaRange($theta);
    return $self->{rx} * (0 + sprintf("%.14f",-sin($theta))) * cos($self->{phi_radians}) + $self->{ry} * (0 + sprintf("%.14f",cos($theta))) * -sin($self->{phi_radians});
}
sub evalYDoublePrimeofTheta {
    my $self = shift;
    my $theta = $self->normalizedThetaToArcTheta(shift);
    return if !$self->isWithinThetaRange($theta);
    return $self->{rx} * (0 + sprintf("%.14f",-cos($theta))) * sin($self->{phi_radians}) + $self->{ry} * (0 + sprintf("%.14f",-sin($theta))) *  cos($self->{phi_radians});
}
sub evalXDoublePrimeofTheta {
    my $self = shift;
    my $theta = $self->normalizedThetaToArcTheta(shift);
    return if !$self->isWithinThetaRange($theta);
    return $self->{rx} * (0 + sprintf("%.14f",-cos($theta))) * cos($self->{phi_radians}) + $self->{ry} * (0 + sprintf("%.14f",-sin($theta))) * -sin($self->{phi_radians});
}
sub solveYforTheta {
    my $self = shift;
    my $y = shift;
    my @xs = $self->F($y);
    my @ts;
    for (my $i=0;$i<@xs;$i++) {
        my $unrot=_rotate2d([$self->{cx},$self->{cy}],[$xs[$i],$y],-$self->{phi_radians});
        $unrot->[0] -= $self->{cx};
        $unrot->[1] -= $self->{cy};
        #now corresponding ellipse formulas are
        #x=$self->{rx} * cos(theta);
        #y=$self->{ry} * sin(theta);
        #so that's solvable for theta, right?
        #y=$self->{ry} * sin(theta);
        #theta=asin(y/ry);

        my $t=asin($unrot->[1]/$self->{ry});

# SIMILAR TO WHAT WE DEBUGGED FIXED FOR solveXforTheta BELOW, BUT UNTESTED/VERIFIED HERE SO FAR
# 5/10/2017
        $t*=-1 if $self->{sweep_flag};


        my $other_t=$t + (($t>0)?-1:1) * $pi;
        my $other_t2=$t + (($t>0)?-1:1) * 2*$pi;
        my $other_t3=$other_t + (($other_t>0)?-1:1) * 2*$pi;
        push(@ts,$t,$other_t,$other_t2,$other_t3);
    }
    return map {$self->arcThetaToNormalizedTheta($_)} grep {$self->isWithinThetaRange($_) && abs($self->evalYofTheta($_) - $y)<0.0000001} @ts;
}
sub solveXforTheta {
    my $self = shift;
    my $x = shift;
    my @ys = $self->f($x);
    my @ts;
    for (my $i=0;$i<@ys;$i++) {
        my $unrot=_rotate2d([$self->{cx},$self->{cy}],[$x,$ys[$i]],-$self->{phi_radians});
        #warn "unrot same?:\n$unrot->[0],$unrot->[1] sameas\n$x,$ys[$i]\n";
        $unrot->[0] -= $self->{cx};
        $unrot->[1] -= $self->{cy};
        #warn "aftershift?:\n$unrot->[0],$unrot->[1]\n [center was $self->{cx},$self->{cy}]\n";

        #now corresponding ellipse formulas are
        #x=$self->{rx} * cos(theta);
        #y=$self->{ry} * sin(theta);
        #so that's solvable for theta, right?
        #x=$self->{rx} * cos(theta);
        #theta=acos(x/rx);
        #theta=pi/2 - asin(x/rx)

# DO WE NEED TO DO A DIFF CALC IF SWEEP FLAG IS ONE WAY OR THE OTHER?
# 5/10/2017

        my $t=($pi/2) - asin($unrot->[0]/$self->{rx});

# NOT SURE ABOUT THIS, BUT MAYBE - looking good! for test problem case
# could be I only ever used this before with sweep_flag == 0 type arcs?
# 5/10/2017
        $t*=-1 if $self->{sweep_flag};

# TODO
# MAKE SURE ALL THIS IS SIMILAR IN solveYforTheta ABOVE AND IN JAVASCRIPT VERSION

        my $other_t  =       $t + (($t      >0) ?-1:1) *   $pi;
        my $other_t2 =       $t + (($t      >0) ?-1:1) * 2*$pi;
        my $other_t3 = $other_t + (($other_t>0) ?-1:1) * 2*$pi;

        push(@ts,$t,$other_t,$other_t2,$other_t3);
    }

    #warn "theta candidates: \n",
    #    join("\n",map {
    #        $_.' nrm:'.$self->arcThetaToNormalizedTheta($_).' '
    #        .($self->isWithinThetaRange($_)?'win':'NOTwin').' ['.$self->{theta1}.', '.$self->{theta2}.'] '
    #        .abs($self->evalXofArcTheta($_) - $x).'<?0.0000001 '
    #        ."\n".$self->evalXofArcTheta($_).','.$self->evalYofArcTheta($_)."\n"
    #        } @ts),"\n";

    return map  {$self->arcThetaToNormalizedTheta($_)}
           grep {   $self->isWithinThetaRange($_)
                 && abs($self->evalXofArcTheta($_) - $x)<0.0000001
           } @ts;
}
sub solveYPrimeforTheta {
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
    return map {$self->arcThetaToNormalizedTheta($_)} grep {$self->isWithinThetaRange($_)} ($t,$other_t,$other_t2,$other_t3);
}
sub solveXPrimeforTheta {
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
    return map {$self->arcThetaToNormalizedTheta($_)} grep {$self->isWithinThetaRange($_)} ($t,$other_t,$other_t2,$other_t3);
}

sub getLength {
    my $self = shift;
    my $res=shift;
    my $start_theta = shift;
    my $end_theta = shift;
    #also need to take theta range options for arc lengths
    if (!defined($res)) {$res=1000;}
#    if (!defined($start_theta)) {$start_theta=$self->{theta1};}
#    if (!defined($end_theta)) {$end_theta=$self->{theta2};}
    if (!defined($start_theta)) {$start_theta=0;}
    if (!defined($end_theta)) {$end_theta=1;}



    # if the two radii are equal, it's a circular arc
    # (which is usually how I use this ellipse stuff)
    # and we have 10th grade math for that
    if ($self->{rx} eq $self->{ry}) {
        #warn "arc length: ",($self->{rx} * $self->{delta_theta} * ($end_theta - $start_theta))," = $self->{rx} * $self->{delta_theta} * ($end_theta - $start_theta)";
        return abs($self->{rx} * $self->{delta_theta} * ($end_theta - $start_theta));
    }

# this isn't set up yet to take arbitrary thetas for sub-lengths of ellipse arc

    my $sum=0;
    my $point1=$self->point($start_theta);
    my $point2;
# here's the problem: for anything but the full arc from theta1 to theta2, using $self->{delta_theta} is wrong
# - you need to calc a different delta - and if you use full ellipse angles, you don't know right off if
# you should use the big delta or the small one - and you would need to screen the angles passed in,
# and what if one of them was out of range? etc.
#    my $thetainc=$self->{delta_theta}/$res;
    my $thetainc=1/$res;
    for (my $i=0;$i<$res;$i++) {
        $point2=$point1;
        #$point1=$self->point($self->{theta1}+$i*$thetainc);
        $point1=$self->point($start_theta+$i*$thetainc);
        $sum += sqrt(($point2->[0]-$point1->[0])**2 + ($point2->[1]-$point1->[1])**2);
    }
    return $sum;
#Those comments might not be right. I think this might work now with the alt parameterization I did for arcs. Still suspect until verified working, but glancing at it now, seems like it might be right and working.
    }

sub getFeet {
    my $self=shift;
    my $x=shift;
    my $y=shift;
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
sub angleTangent_byTheta {
    my $self = shift;
    my $theta = shift;
    return if !$self->isWithinThetaRange($self->normalizedThetaToArcTheta($theta));
    my @ret=$self->angleTangent(undef,undef,$theta);
    return wantarray ? @ret : $ret[0];
    }
sub angleNormal_byTheta {
    my $self = shift;
    my $theta = shift;
    return if !$self->isWithinThetaRange($self->normalizedThetaToArcTheta($theta));
    my @ret=$self->angleNormal(undef,undef,$theta);
    return wantarray ? @ret : $ret[0];
    }
sub slopeTangent_byTheta {
    my $self = shift;
    my $theta = shift;
    return if !$self->isWithinThetaRange($self->normalizedThetaToArcTheta($theta));
    my $yp=$self->evalYPrimeofTheta($theta);
    my $xp=$self->evalXPrimeofTheta($theta);
    if ($xp eq '0') {return (($yp < 0)?'-':'+').'inf';}
    if ($yp eq '0') {return (($xp < 0)?'-':'+').'0';}
    else {
        return $yp/$xp;
    }
}

sub slopeNormal_byTheta {
    my $self = shift;
    my $theta = shift;
    return if !$self->isWithinThetaRange($self->normalizedThetaToArcTheta($theta));
    my $slopeTangent = $self->slopeTangent_byTheta($theta);
    if ($slopeTangent =~ /([\-\+]?)inf$/i) {
        my $sign = '';
        if (length($1)) {if ($1 eq '-') {$sign='+';} else {$sign='-';}}
        return $sign.'0';
    }
    elsif ($slopeTangent =~ /^([\-\+]?)0$/) {
        my $sign = '';
        if (length($1)) {if ($1 eq '-') {$sign='+';} else {$sign='-';}}
        return $sign.'inf';
    }
    else {
        return -1/$slopeTangent;
    }
}
sub angleTangent {
    my $self = shift;
    my $x    = shift;
    my $y    = @_?shift:undef;
    my $theta= @_?shift:undef;
    #get intersect points
    my @intersects;
    if (defined($x)) {my @ys;@ys=$self->f($x);foreach (sort {$a<=>$b} @ys) {push(@intersects,[$x,$_])}}
    elsif (defined($y)) {my @xs;@xs=$self->F($y);foreach (sort {$a<=>$b} @xs) {push(@intersects,[$_,$y])}}
    elsif (defined($theta)) {push(@intersects,$self->point($theta));}
    #then use the foci calculated in the ellipse setup
    # and the info here: http://mysite.du.edu/~jcalvert/math/ellipse.htm
    # to make lines and figure angles to get tangent angle...

    #The tangent line at point P on the ellipse is perpendicular to the line bisecting
    #the angle between the two lines extending from point P to the two foci of the ellipse. (So the bisector is the normal.)
    #That angle is given by the arccosine of the dot product over the product of the magnitude of the vectors (lines) between the two lines:
    # arccos( (line1 dot line2) / (|line1|*|line2|) )
    #arccos(x) is eqivalent to pi/2 - arcsin(x)

    #really you're calculating an inward pointing normal angle and adding 90 deg to get the tangent

    # ... much later : but you should take sweep_flag into account
    #    The elliptical arcs here have direction - start point and end point - and a tangent should go in start-to-end direction of elliptical path
    #    and a normal should point off to the right (left! right in +y points down coord sys, but left in a +y goes up coord sys)
    #    so added one sweep flag controlled * 1 or -1 in this stuff, and took a negative sign off the cos/sin in slopeNormal function
    #    and that looks right on my test page

    my @ret;
    for (my $i=0;$i<@intersects;$i++) {
        my $line1=[$intersects[$i],$self->{f1}];
        my $line2=[$intersects[$i],$self->{f2}];
        my $a1 = atan2($line1->[1]->[1] - $line1->[0]->[1],$line1->[1]->[0] - $line1->[0]->[0]);
        my $a2 = atan2($line2->[1]->[1] - $line2->[0]->[1],$line2->[1]->[0] - $line2->[0]->[0]);
        push(@ret,
            $pi/2 * ($self->{sweep_flag}?-1:1) + #add +/- 90 deg from the normal angle you calculate below to get the tangent
            $a1 + #angle of the line/vector from point P on ellipse to focus 1
            (($a2 - $a1)>0 || ($a2 - $a1) eq 0?1:-1) * # hmm..., whether we add or subtract the half angle below from the line1 angle. Is this okay?
            0.5 * # this and below calculated half the angle between the two lines
            ($pi/2 -
                asin(
                        (
                            ($line1->[1]->[0] - $line1->[0]->[0]) * ($line2->[1]->[0] - $line2->[0]->[0])  +
                            ($line1->[1]->[1] - $line1->[0]->[1]) * ($line2->[1]->[1] - $line2->[0]->[1])
                        ) /
                        (
                            sqrt(($line1->[1]->[1] - $line1->[0]->[1])**2 + ($line1->[1]->[0] - $line1->[0]->[0])**2) *
                            sqrt(($line2->[1]->[1] - $line2->[0]->[1])**2 + ($line2->[1]->[0] - $line2->[0]->[0])**2)
                        )
                )
            )
        );
    }
    @ret = map {angle_reduce($_)} @ret;
    return wantarray ? @ret : $ret[0];
}
sub secondDerivative {
    my $self = shift;
    die "Need to work out second derivative with respect to x for elliptical arc when rx != ry";
}
sub slopeTangent {my @ats;@ats=$_[0]->angleTangent($_[1],$_[2],$_[3]);my @ret;for (my $i=0;$i<@ats;$i++) {push(@ret, sin($ats[$i])/cos($ats[$i]))} return wantarray ? @ret : $ret[0];}
sub slopeNormal  {my @ats;@ats=$_[0]->angleTangent($_[1],$_[2],$_[3]);my @ret;for (my $i=0;$i<@ats;$i++) {push(@ret, cos($ats[$i])/sin($ats[$i]))} return wantarray ? @ret : $ret[0];}
sub angleNormal  {my @ret = (map {$_ + $pi/2} $_[0]->angleTangent($_[1],$_[2],$_[3]));@ret = map {angle_reduce($_)} @ret;return wantarray ? @ret : $ret[0];}

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
sub _rotate2d {
    my $origin = shift;
    my $point = shift;
    my $angle = shift;
    my $dx=($point->[0]-$origin->[0]);
    my $dy=($point->[1]-$origin->[1]);
    #{a c-b d, a d+b c}
    return [$origin->[0] + ($dx*cos($angle) - $dy*sin($angle)),$origin->[1] + ($dx*sin($angle) + $dy*cos($angle))];
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

    my $self=shift;

    my $dim=shift;
    my $theta=shift;
    my $direction=scalar(@_)?shift:1; # 0 or 1

    my $pt_last = $self->point($theta); # shouldn't this be outside the function def?

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