####################################################################################
###      Math::MPath::BezierCubicSegment   ######################################
package Math::MPath::BezierCubicSegment;
{
use POSIX qw(nexttoward); # lets you step floats by smallest possible increment
use Math::MPath::CubicFormula qw(cubicformula);
use Math::MPath::QuadraticFormula qw(quadraticformula);
use Math::MPath::Function::Root qw(BrentsMethod Bisection FalsePosition);
our $pi = 4 * atan2(1,1);
our $pi_big = Math::BigFloat->new('3.141592653589793238462643383279502884197169399375105820974944592308'); # for use possibly with new built-in cubic solver, unless we end up factoring it out
our $piovertwo = 2 * atan2(1,1);
our $piovertwo_big = $pi_big->copy()->bdiv(2); # for use possibly with new built-in cubic solver, unless we end up factoring it out
our $twopi = 8 * atan2(1,1);
our $twopi_big = $pi_big->copy()->bmul(2); # for use possibly with new built-in cubic solver, unless we end up factoring it out
our $fourpi = 16 * atan2(1,1);
our $fourpi_big = $pi_big->copy()->bmul(4); # for use possibly with new built-in cubic solver, unless we end up factoring it out

sub new {
    my $class = shift;
    my $self={};
    bless $self,$class;
    $self->{p1} = shift;
    $self->{cp1}= shift;
    $self->{cp2}= shift;
    $self->{p2} = shift;
    $self->{precision} = shift;
    $self->{isLite} = @_ ? shift:0;
    $self->{maxdiaglength} = (reverse sort {$a<=>$b} map {sqrt(($_->[3] - $_->[2])**2 + ($_->[1] - $_->[0])**2)} ([@{$self->{p1}},@{$self->{cp2}}],[@{$self->{cp1}},@{$self->{p2}}]))[0];
    $self->{A}  = $self->{p2}->[0] - 3.0 * $self->{cp2}->[0] + 3.0 * $self->{cp1}->[0] -       $self->{p1}->[0];
    $self->{B}  =                    3.0 * $self->{cp2}->[0] - 6.0 * $self->{cp1}->[0] + 3.0 * $self->{p1}->[0];
    $self->{C}  =                                              3.0 * $self->{cp1}->[0] - 3.0 * $self->{p1}->[0];
    $self->{D}  =                                                                              $self->{p1}->[0];
    $self->{E}  = $self->{p2}->[1] - 3.0 * $self->{cp2}->[1] + 3.0 * $self->{cp1}->[1] -       $self->{p1}->[1];
    $self->{F}  =                    3.0 * $self->{cp2}->[1] - 6.0 * $self->{cp1}->[1] + 3.0 * $self->{p1}->[1];
    $self->{G}  =                                              3.0 * $self->{cp1}->[1] - 3.0 * $self->{p1}->[1];
    $self->{H}  =                                                                              $self->{p1}->[1];
    my $bA;
    if (abs($self->{A}) < 0.0001) {$bA=Math::BigFloat->new(''.$self->{A})} #just discovered this is a rare problem that pushes the limits of arctan at +/-1 in CubicFormula() stuff
    if ($self->{B} eq 0) {$self->{BdA}=0;} else {$self->{BdA}=(!defined($bA))?$self->{B}/$self->{A} : (Math::BigFloat->new(''.$self->{B}))->bdiv($bA) ;}
    if ($self->{C} eq 0) {$self->{CdA}=0;} else {$self->{CdA}=(!defined($bA))?$self->{C}/$self->{A} : (Math::BigFloat->new(''.$self->{C}))->bdiv($bA) ;}
    my $bE;
    if (abs($self->{E}) < 0.0001) {$bE=Math::BigFloat->new(''.$self->{E})} # same as for A above ....  you know, if E or A somehow came out to exactly zero... what kind of bezier does that?... does it degrade to quadratic bezier? or to a circle/arc formula? worth playing with at some point
    if ($self->{F} eq 0) {$self->{FdE}=0;} else {$self->{FdE}=(!defined($bE))?$self->{F}/$self->{E} : (Math::BigFloat->new(''.$self->{F}))->bdiv($bE) ;}
    if ($self->{G} eq 0) {$self->{GdE}=0;} else {$self->{GdE}=(!defined($bE))?$self->{G}/$self->{E} : (Math::BigFloat->new(''.$self->{G}))->bdiv($bE) ;}

    $self->{Am3}=$self->{A} * 3.0;
    $self->{Bm2}=$self->{B} * 2.0;
    $self->{Em3}=$self->{E} * 3.0;
    $self->{Fm2}=$self->{F} * 2.0;
    $self->{Am6}=$self->{A} * 6.0;
    $self->{Em6}=$self->{E} * 6.0;

    # TODO might be nice to try to replace this with lazy init per big needed?
    if (!$self->{isLite}) {
        $self->initBigs();
        }

    my @extremexs_is=(0,((!$self->{isLite})?$self->solveXPrimeforThetaBig(Math::BigFloat->bzero()):$self->solveXPrimeforTheta(0)),1);
    my @extremexs = map {(ref($_) && !$self->{isLite}) ? $self->bezierEvalXofTBig($_):$self->bezierEvalXofT($_)} @extremexs_is;
    my @extremexs_sorted = sort {$a<=>$b} @extremexs;
    my @extremeys_is=(0,((!$self->{isLite})?$self->solveYPrimeforThetaBig(Math::BigFloat->bzero()):$self->solveYPrimeforTheta(0)),1);
    my @extremeys = map {(ref($_) && !$self->{isLite}) ? $self->bezierEvalYofTBig($_):$self->bezierEvalYofT($_)} @extremeys_is;
    my @extremeys_sorted = sort {$a<=>$b} @extremeys;

    $self->{extremexs_is}=[(@extremexs_is)];
    $self->{extremeys_is}=[(@extremeys_is)];
    $self->{extremexs}=[(@extremexs)];
    $self->{extremeys}=[(@extremeys)];
    $self->{minx} = $extremexs_sorted[0];
    $self->{maxx} = $extremexs_sorted[$#extremexs];
    $self->{miny} = $extremeys_sorted[0];
    $self->{maxy} = $extremeys_sorted[$#extremeys];

    #FOLLOWING (and above and other repeated in initDangerRanges function)
    #profiled slice making task with DProfLB and found that
    #20% to 60% of time is taking doing inRange checks for beziers
    #probably because those checks are using these min/max x/y vals
    #that are BigFloats. So let's make regular PERL scalar number versions
    #and use those in the inRange tests. Will probably reintroduce
    #some boundary problems, but it might also give a 50%+ speedup.

    # numify if refs (BigFloats)
    $self->{minx_notbig} = ref($self->{minx}) ? 0 + sprintf("%.20f",$self->{minx}->bstr()):$self->{minx};
    $self->{maxx_notbig} = ref($self->{maxx}) ? 0 + sprintf("%.20f",$self->{maxx}->bstr()):$self->{maxx};
    $self->{miny_notbig} = ref($self->{miny}) ? 0 + sprintf("%.20f",$self->{miny}->bstr()):$self->{miny};
    $self->{maxy_notbig} = ref($self->{maxy}) ? 0 + sprintf("%.20f",$self->{maxy}->bstr()):$self->{maxy};


    # An old comment about trying to deal with precision problems in special cases
    # (reveals not quite fully formed understanding of how to deal with that stuff):

    # There are two cases you need to guard against
    #1. f'(x) or F'(y) are too small to produce unique answers using specified and available precision
    #2. m_ix or m_iy are too small to give unique answers for current-precision-sized increments of x or y
    #I think the best way to ensure good results while maintaining most of the speed that results from ignoring this issue
    #is to determine ahead of time what level of precision is needed over different intervals of the xy curve. And then
    #ensure that that level of precision is in effect when each point lookup is made.
    #There will be a generally safe intervals where the required precision will usually be present. Here you will simply check for presence of required precision, and then upgrade if necessary.
    #In the more dangerous intervals, we will immediately go into calculating how much extra precision is needed.
    #Default for PERL appears to be to use DOUBLEs, whose constant accuracy you can roughly translate into current precision
    #for a given number by subtracting the number of digits left of decimal from total avail (15 I think).
    #Another thing to consider in making safe and unsafe intervals is, then, magnitude of x or y on the left side of decimal.
    #Interval[0,99] has 14 digit precision. Interval [100,999] has 13.
    #It would be nice to be able to pass PERL native floats into calculations, but still be able to demand calculation to certain precision, which may exceed PERL float native precision.
    #This would allow delay of upgrade to some higher precision number format. If the calculation subroutine implements its own method of doing exact arithmatic with the PERL native floats, those passed-in floats could stay PERL native (and fast, right?).
    #Only return vals would have to be made arb. precision.

    if (!$self->{isLite}) {
        $self->initDangerRanges();
        }


    # The new decompostion of the cubic Bezier into piecewise functions of x

    $self->{cubQ} = (3.0*$self->{CdA} - $self->{BdA}**2)/9.0;
    $self->{cubsqrtofnegQ} = sqrt(($self->{cubQ} < 0 ? -1:1) * $self->{cubQ}); # only used in cubD < 0 case, where cubQ is always < 0

    $self->{cubQ_ofy} = (3.0*$self->{GdE} - $self->{FdE}**2)/9.0;
    $self->{cubsqrtofnegQ_ofy} = sqrt(($self->{cubQ_ofy} < 0 ? -1:1) * $self->{cubQ_ofy}); # only used in cubD < 0 case, where cubQ is always < 0


    # These subs seem to be working
    my $theta = sub { # acos( R / sqrt( -(Q^3) ) )
                     my $toAcos =
                      # R
                      (
                        (   9.0 * $self->{BdA}*$self->{CdA}
                          - 27.0 *
                                   #C
                                   (($self->{D} - $_[0]) / $self->{A})

                          - 2.0  * $self->{BdA}**3
                        )
                        / 54.0
                      )

                      /

                      # sqrt( -(Q^3) )
                      sqrt(-(((3.0*$self->{CdA} - $self->{BdA}**2)/9.0)**3));

                     #return acos($toAcos);

                     $toAcos = 1.0 if $toAcos > 1;
                     $toAcos = -1.0 if $toAcos < -1;

                     #      this is arccosine
                     return atan2(sqrt(1.0 - $toAcos * $toAcos),$toAcos);

    };

    my $theta_ofy = sub { # acos( R / sqrt( -(Q^3) ) )
                     my $toAcos =
                      # R
                      (
                        (   9.0 * $self->{FdE}*$self->{GdE}
                          - 27.0 *
                                   #C
                                   (($self->{H} - $_[0]) / $self->{E})

                          - 2.0  * $self->{FdE}**3
                        )
                        / 54.0
                      )

                      /

                      # sqrt( -(Q^3) )
                      sqrt(-(((3.0*$self->{GdE} - $self->{FdE}**2)/9.0)**3));

                     #return acos($toAcos);

                     $toAcos = 1.0 if $toAcos > 1;
                     $toAcos = -1.0 if $toAcos < -1;

                     #      this is arccosine
                     return atan2(sqrt(1.0 - $toAcos * $toAcos),$toAcos);

    };


    my $D = sub { return
                       ( (3*$self->{CdA} - $self->{BdA}**2)/9.0 )**3 #   Q^3
                     + (                                               # + R^2
                        # R
                        (   9.0  * $self->{BdA}*$self->{CdA}
                          - 27.0 *

                                 #C
                                 (($self->{D} - $_[0]) / $self->{A})

                          - 2.0  * $self->{BdA}**3
                        )
                        / 54.0
                       )**2;
    };

    my $D_ofy = sub { return
                       ( (3*$self->{GdE} - $self->{FdE}**2)/9.0 )**3 #   Q^3
                     + (                                               # + R^2
                        # R
                        (   9.0  * $self->{FdE}*$self->{GdE}
                          - 27.0 *

                                 #C
                                 # TODO : does returning Inf here do what you expect in every case? Should you consider -Inf too?
                                 ( $self->{E} == 0 ? Inf : (($self->{H} - $_[0]) / $self->{E}) )

                          - 2.0  * $self->{FdE}**3
                        )
                        / 54.0
                       )**2;
    };

    my $sqrtD = sub {
                     my $tosqrt =
                       ( (3*$self->{CdA} - $self->{BdA}**2)/9.0 )**3 #   Q^3
                     + (                                               # + R^2
                        # R
                        (   9.0  * $self->{BdA}*$self->{CdA}
                          - 27.0 *

                                 #C
                                 (($self->{D} - $_[0]) / $self->{A})

                          - 2.0  * $self->{BdA}**3
                        )
                        / 54.0
                       )**2;
                     #$tosqrt = 0 if abs($tosqrt) < 1e-14;

if ($tosqrt<0) {
my $x=$_[0];
#warn "1-4 fixup\n";
            my $d1 = $D->($x=nexttoward($x,$self->{maxx}));
            my $dir = ($d1>$tosqrt) ? $self->{maxx} : $self->{minx};
            while ($tosqrt<0) {$tosqrt = $D->($x=nexttoward($x,$dir)); }
}

                     return (ref($tosqrt)?bigsqrt($tosqrt):sqrt($tosqrt));
                    };

    my $sqrtD_ofy = sub {
                     my $tosqrt =
                       ( (3*$self->{GdE} - $self->{FdE}**2)/9.0 )**3 #   Q^3
                     + (                                               # + R^2
                        # R
                        (   9.0  * $self->{FdA}*$self->{GdE}
                          - 27.0 *

                                 #C
                                 (($self->{H} - $_[0]) / $self->{E})

                          - 2.0  * $self->{FdE}**3
                        )
                        / 54.0
                       )**2;
                     #$tosqrt = 0 if abs($tosqrt) < 1e-14;

if ($tosqrt<0) {
my $y=$_[0];
#warn "1-4 fixup\n";
            my $d1 = $D_ofy->($y=nexttoward($y,$self->{maxy}));
            my $dir = ($d1>$tosqrt) ? $self->{maxy} : $self->{miny};
            while ($tosqrt<0) {$tosqrt = $D_ofy->($y=nexttoward($y,$dir)); }
}

                     return (ref($tosqrt)?bigsqrt($tosqrt):sqrt($tosqrt));
                    };


    my $R = sub {
                 return
                 (   9.0  * $self->{BdA}*$self->{CdA}
                    - 27.0 *

                             #C
                             (($self->{D} - $_[0]) / $self->{A})

                    - 2.0  * $self->{BdA}**3
                  )
                  / 54.0
    };
    my $R_ofy = sub {
                 return
                 (   9.0  * $self->{FdE}*$self->{GdE}
                    - 27.0 *

                             #C
                             # TODO : does returning Inf here do what you expect in every case? Should you consider -Inf too?
                             ( $self->{E} == 0 ? Inf : (($self->{H} - $_[0]) / $self->{E}))

                    - 2.0  * $self->{FdE}**3
                  )
                  / 54.0
    };

    my $preS = sub {return $R->($_[0])+$sqrtD->($_[0])};
    my $preT = sub {return $R->($_[0])-$sqrtD->($_[0])};
    my $preS_ofy = sub {return $R_ofy->($_[0])+$sqrtD_ofy->($_[0])};
    my $preT_ofy = sub {return $R_ofy->($_[0])-$sqrtD_ofy->($_[0])};

    my $eqn_1 = sub {return -($self->{BdA}/3) + (  (( ($R->($_[0])+$sqrtD->($_[0])))**(1/3)) + ( ( ($R->($_[0])-$sqrtD->($_[0])))**(1/3)) );};
    my $eqn_2 = sub {return -($self->{BdA}/3) + (  (( ($R->($_[0])+$sqrtD->($_[0])))**(1/3)) + (-(-($R->($_[0])-$sqrtD->($_[0])))**(1/3)) );}; # yes 2 and 4 are the same
    my $eqn_3 = sub {return -($self->{BdA}/3) + ( -((-($R->($_[0])+$sqrtD->($_[0])))**(1/3)) + (-(-($R->($_[0])-$sqrtD->($_[0])))**(1/3)) );};
    my $eqn_4 = sub {return -($self->{BdA}/3) + (  (( ($R->($_[0])+$sqrtD->($_[0])))**(1/3)) + (-(-($R->($_[0])-$sqrtD->($_[0])))**(1/3)) );};

    my $eqn_1_ofy = sub {return -($self->{FdA}/3) + (  (( ($R_ofy->($_[0])+$sqrtD_ofy->($_[0])))**(1/3)) + ( ( ($R_ofy->($_[0])-$sqrtD_ofy->($_[0])))**(1/3)) );};
    my $eqn_2_ofy = sub {return -($self->{FdA}/3) + (  (( ($R_ofy->($_[0])+$sqrtD_ofy->($_[0])))**(1/3)) + (-(-($R_ofy->($_[0])-$sqrtD_ofy->($_[0])))**(1/3)) );}; # yes 2 and 4 are the same
    my $eqn_3_ofy = sub {return -($self->{FdA}/3) + ( -((-($R_ofy->($_[0])+$sqrtD_ofy->($_[0])))**(1/3)) + (-(-($R_ofy->($_[0])-$sqrtD_ofy->($_[0])))**(1/3)) );};
    my $eqn_4_ofy = sub {return -($self->{FdA}/3) + (  (( ($R_ofy->($_[0])+$sqrtD_ofy->($_[0])))**(1/3)) + (-(-($R_ofy->($_[0])-$sqrtD_ofy->($_[0])))**(1/3)) );};

    #my $eqn_8 = sub {return one of above + ($self->{BdA}/3) then /2, then negated, then minus - ($self->{BdA}/3) } # the D eq 0 duplicate real root

    my $eqn_1to4_prime  = sub {
        my $x=$_[0];
        my $preS = $preS->($x);
        my $preT = $preT->($x);
        my $thissqrtD = $sqrtD->($x);
        return Inf if $thissqrtD == 0;
        return
        (1/(6.0  * $thissqrtD*$self->{A}))
        * (  ($preS/(($preS**2)**(1/3)))
           - ($preT/(($preT**2)**(1/3)))
          )
        ;
    };

    my $eqn_1to4_2prime = sub {
        my $x=$_[0];
        my $preS = $preS->($x);
        my $preT = $preT->($x);
        my $sqrtD = $sqrtD->($x);
        my $R = $R->($x);
        return
          (1/(12.0 * $D->($x)    *$self->{A}**2))
        * (  ($preS/(($preS**2)**(1/3)))
           * ( (1/3) - $R/($sqrtD) )
           +
             ($preT/(($preT**2)**(1/3)))
           * ( (1/3) + $R/($sqrtD) )
          )
        ;
    };


    my $eqn_1to4_prime_ofy  = sub {
        my $preS = $preS_ofy->($_[0]);
        my $preT = $preT_ofy->($_[0]);
        my $thissqrtD = $sqrtD_ofy->($_[0]);
        return Inf if $thissqrtD == 0;
        return
        (1/(6.0  * $thissqrtD*$self->{E}))
        * (  ($preS/(($preS**2)**(1/3)))
           - ($preT/(($preT**2)**(1/3)))
          )
        ;
    };

    my $eqn_1to4_2prime_ofy = sub {
        my $preS = $preS_ofy->($_[0]);
        my $preT = $preT_ofy->($_[0]);
        my $sqrtD = $sqrtD_ofy->($_[0]);
        my $R = $R_ofy->($_[0]);
        return
          (1/(12.0 * $D_ofy->($_[0]) * $self->{E}**2))
        * (  ($preS/(($preS**2)**(1/3)))
           * ( (1/3) - $R/($sqrtD) )
           +
             ($preT/(($preT**2)**(1/3)))
           * ( (1/3) + $R/($sqrtD) )
          )
        ;
    };



    my $eqn_5 = sub {return 2.0*sqrt(-((3.0*$self->{CdA} - $self->{BdA}**2)/9.0)) * cos( ($theta->($_[0])          ) / 3.0) - $self->{BdA}/3.0};
    my $eqn_6 = sub {return 2.0*sqrt(-((3.0*$self->{CdA} - $self->{BdA}**2)/9.0)) * cos( ($theta->($_[0]) + $twopi ) / 3.0) - $self->{BdA}/3.0};
    my $eqn_7 = sub {return 2.0*sqrt(-((3.0*$self->{CdA} - $self->{BdA}**2)/9.0)) * cos( ($theta->($_[0]) + $fourpi) / 3.0) - $self->{BdA}/3.0};

    my $eqn_5y = sub {return 2.0*sqrt(-((3.0*$self->{GdA} - $self->{FdA}**2)/9.0)) * cos( ($theta_ofy->($_[0])          ) / 3.0) - $self->{FdA}/3.0};
    my $eqn_6y = sub {return 2.0*sqrt(-((3.0*$self->{GdA} - $self->{FdA}**2)/9.0)) * cos( ($theta_ofy->($_[0]) + $twopi ) / 3.0) - $self->{FdA}/3.0};
    my $eqn_7y = sub {return 2.0*sqrt(-((3.0*$self->{GdA} - $self->{FdA}**2)/9.0)) * cos( ($theta_ofy->($_[0]) + $fourpi) / 3.0) - $self->{FdA}/3.0};

    my $eqn_5_prime  = sub {
        my $x=$_[0] - $self->{minx}*0;
        my $thisD = $D->($x);

        if ($thisD>0) {
#warn "inf5 ";
#return ($self->{cubsqrtofnegQ}/(3.0*$self->{A})) * sin($theta->($x)/3.0) > 0 ? Inf:-Inf;
            my $d1 = $D->($x=nexttoward($x,$self->{maxx}));
            while ($d1 == $thisD) {$d1 = $D->($x=nexttoward($x,$self->{maxx}));}
            my $dir = ($d1<$thisD) ? $self->{maxx} : $self->{minx};
            while ($thisD>0) {$thisD = $D->($x=nexttoward($x,$dir)); }
            }

        return ($self->{cubsqrtofnegQ}/(3.0*$self->{A})) * (sin($theta->($x)/3.0) / (sqrt(-$thisD)));
    };
    my $eqn_5_2prime = sub {
        my $x=$_[0] - $self->{minx}*0;
        my $D = $D->($x);
        my $theta = $theta->($x);
        return
        -1 * # because I didn't keep good enough track of the square root of -1 while deriving this
        ($self->{cubsqrtofnegQ}/(6.0 * $D * $self->{A}**2))
        *
        ( (sin($theta/3.0) / sqrt(-$D))  * $R->($x) - cos($theta/3.0)/3.0 )
        ;
    };

    my $eqn_6_prime  = sub {
        my $x=$_[0] - $self->{minx}*0;
        my $thisD = $D->($x);

        if ($thisD>0) {
#warn "inf6 ";
#return ($self->{cubsqrtofnegQ}/(3.0*$self->{A})) * sin(($theta->($x) + $twopi)/3.0) > 0 ? Inf:-Inf;
            my $d1 = $D->($x=nexttoward($x,$self->{maxx}));
            while ($d1 == $thisD) {$d1 = $D->($x=nexttoward($x,$self->{maxx}));}
            my $dir = ($d1<$thisD) ? $self->{maxx} : $self->{minx};
            while ($thisD>0) {$thisD = $D->($x=nexttoward($x,$dir)); }

            }

        return ($self->{cubsqrtofnegQ}/(3.0*$self->{A})) * (sin(($theta->($x) + $twopi)/3.0) / (sqrt(-$thisD)));
    };
    my $eqn_6_2prime = sub {
        my $x=$_[0] - $self->{minx}*0;
        my $D = $D->($x);
        my $theta = $theta->($x);
        return
        -1 * # because I didn't keep good enough track of the square root of -1 while deriving this
        ($self->{cubsqrtofnegQ}/(6.0 * $D * $self->{A}**2))
        *
        ( (sin(($theta + $twopi)/3.0) / sqrt(-$D))  * $R->($x) - cos(($theta + $twopi)/3.0)/3.0 )
        ;
    };

    my $eqn_7_prime  = sub {
        my $x=$_[0] - $self->{minx}*0;
        my $thisD = $D->($x);

        if ($thisD>0) {
#warn "inf7 ";
#return ($self->{cubsqrtofnegQ}/(3.0*$self->{A})) * sin(($theta->($x) + $fourpi)/3.0) > 0 ? Inf:-Inf;
            my $d1 = $D->($x=nexttoward($x,$self->{maxx}));
            while ($d1 == $thisD) {$d1 = $D->($x=nexttoward($x,$self->{maxx}));}
            my $dir = ($d1<$thisD) ? $self->{maxx} : $self->{minx};
            while ($thisD>0) {$thisD = $D->($x=nexttoward($x,$dir)); }
            }

        return ($self->{cubsqrtofnegQ}/(3.0*$self->{A})) * (sin(($theta->($x) + $fourpi)/3.0) / (sqrt(-$thisD)));
    };
    my $eqn_7_2prime = sub {
        my $x=$_[0] - $self->{minx}*0;
        my $D = $D->($x);
        my $theta = $theta->($x);
        return
        -1 * # because I didn't keep good enough track of the square root of -1 while deriving this
        ($self->{cubsqrtofnegQ}/(6.0 * $D * $self->{A}**2))
        *
        ( (sin(($theta + $fourpi)/3.0) / sqrt(-$D))  * $R->($x) - cos(($theta + $fourpi)/3.0)/3.0 )
        ;
    };


    ### start _ofy versions

    my $eqn_5_prime_ofy  = sub {
        my $y=$_[0];
        my $thisD = $D_ofy->($y);

        if ($thisD>0) {
            my $d1 = $D_ofy->($y=nexttoward($y,$self->{maxy}));
            my $dir = ($d1<$thisD) ? $self->{maxy} : $self->{miny};
            while ($thisD>0) {$thisD = $D_ofy->($y=nexttoward($y,$dir)); }
            }

        return ($self->{cubsqrtofnegQ_ofy}/(3.0*$self->{E})) * (sin($theta_ofy->($y)/3.0) / (sqrt(-$thisD)));
    };
    my $eqn_5_2prime_ofy = sub {
        my $y=$_[0];
        my $D = $D_ofy->($y);
        my $theta = $theta_ofy->($y);
        return
        -1 * # because I didn't keep good enough track of the square root of -1 while deriving this
        ($self->{cubsqrtofnegQ_ofy}/(6.0 * $D * $self->{E}**2))
        *
        ( (sin($theta/3.0) / sqrt(-$D))  * $R_ofy->($y) - cos($theta/3.0)/3.0 )
        ;
    };

    my $eqn_6_prime_ofy  = sub {
        my $y=$_[0];
        my $thisD = $D_ofy->($y);

        if ($thisD>0) {
            my $d1 = $D_ofy->($y=nexttoward($y,$self->{maxy}));
            my $dir = ($d1<$thisD) ? $self->{maxy} : $self->{miny};
            while ($thisD>0) {$thisD = $D_ofy->($y=nexttoward($y,$dir)); }
            }

        return ($self->{cubsqrtofnegQ_ofy}/(3.0*$self->{E})) * (sin(($theta_ofy->($y) + $twopi)/3.0) / (sqrt(-$thisD)));
    };
    my $eqn_6_2prime_ofy = sub {
        my $y=$_[0];
        my $D = $D_ofy->($y);
        my $theta = $theta_ofy->($y);
        return
        -1 * # because I didn't keep good enough track of the square root of -1 while deriving this
        ($self->{cubsqrtofnegQ_ofy}/(6.0 * $D * $self->{E}**2))
        *
        ( (sin(($theta + $twopi)/3.0) / sqrt(-$D))  * $R_ofy->($y) - cos(($theta + $twopi)/3.0)/3.0 )
        ;
    };

    my $eqn_7_prime_ofy  = sub {
        my $y=$_[0];
        my $thisD = $D_ofy->($y);

        if ($thisD>0) {
            my $d1 = $D_ofy->($y=nexttoward($y,$self->{maxy}));
            my $dir = ($d1<$thisD) ? $self->{maxy} : $self->{miny};
            while ($thisD>0) {$thisD = $D_ofy->($y=nexttoward($y,$dir)); }
            }

        return ($self->{cubsqrtofnegQ_ofy}/(3.0*$self->{E})) * (sin(($theta_ofy->($y) + $fourpi)/3.0) / (sqrt(-$thisD)));
    };
    my $eqn_7_2prime_ofy = sub {
        my $y=$_[0];
        my $D = $D_ofy->($y);
        my $theta = $theta_ofy->($y);
        return
        -1 * # because I didn't keep good enough track of the square root of -1 while deriving this
        ($self->{cubsqrtofnegQ_ofy}/(6.0 * $D * $self->{E}**2))
        *
        ( (sin(($theta + $fourpi)/3.0) / sqrt(-$D))  * $R_ofy->($y) - cos(($theta + $fourpi)/3.0)/3.0 )
        ;
    };

    ### done _ofy versions

    my %eqn_names = (
        $eqn_1 => 'eqn_1',
        $eqn_2 => 'eqn_2',
        $eqn_3 => 'eqn_3',
        $eqn_4 => 'eqn_4',
        $eqn_5 => 'eqn_5',
        $eqn_6 => 'eqn_6',
        $eqn_7 => 'eqn_7',
#        $eqn_8 => 'eqn_8', # still havn't done D eq 0 case(s)
    );

    # Find bezier parameters at key points that bound where the above equations apply.

    # there may be more than these ts in the end
    my @all_ts = sort {$a<=>$b} (0 ,
                                 # these are where slope goes infinite vertical
                                 #((!$self->{isLite})?$self->solveXPrimeforThetaBig_noFilter(Math::BigFloat->bzero()):$self->solveXPrimeforTheta_noFilter(0)) ,
                                 # think we went past idea of using out of bounds ts
                                 ((!$self->{isLite})?$self->solveXPrimeforThetaBig(Math::BigFloat->bzero()):$self->solveXPrimeforTheta(0)) ,
                                 1
                                );

    #warn "\nallts\n",join("\n",map {$_.' : '.$self->bezierEvalXofT($_)} @all_ts),"\n";


    # This is ultimately to find ts where D(x(t)) crosses 0 (not just touches, but crosses)
    my $x_d0_1 = ($self->{cubQ}>0)? 'nan' :  2.0*$self->{A}*sqrt(-( (3.0*$self->{CdA} - $self->{BdA}**2)/9.0 )**3) - $self->{A}*($self->{BdA}*$self->{CdA})/(3.0) + $self->{A}*(2.0*$self->{BdA}**3)/(27.0) + $self->{D};
    my $x_d0_2 = ($self->{cubQ}>0)? 'nan' : -2.0*$self->{A}*sqrt(-( (3.0*$self->{CdA} - $self->{BdA}**2)/9.0 )**3) - $self->{A}*($self->{BdA}*$self->{CdA})/(3.0) + $self->{A}*(2.0*$self->{BdA}**3)/(27.0) + $self->{D};

    # same for dy/dt == 0 stuff - but not using yet
    #my $y_d0_1 =  2.0*$self->{E}*sqrt(-( (3.0*$self->{GdE} - $self->{FdE}**2)/9.0 )**3) - $self->{E}*($self->{FdE}*$self->{GdE})/(3.0) + $self->{E}*(2.0*$self->{FdE}**3)/(27.0) + $self->{H};
    #my $y_d0_2 = -2.0*$self->{E}*sqrt(-( (3.0*$self->{GdE} - $self->{FdE}**2)/9.0 )**3) - $self->{E}*($self->{FdE}*$self->{GdE})/(3.0) + $self->{E}*(2.0*$self->{FdE}**3)/(27.0) + $self->{H};

    #warn "x_d0s: $x_d0_1, $x_d0_2\n";

    my @more_t;
    my @more_x;

    my $whilewatch = 0;

    if ($x_d0_1 && $x_d0_1 !~ /nan/i && $x_d0_1 > $self->{minx} && $x_d0_1 < $self->{maxx}) {
        # then figure and add the right xs and ts to capture that D(x)==0 crossing
        # The "zeroing" xs you have, if in bez xrange, probably don't exactly zero the D(x) function
        # If D(x) > 0, you figure which of eqns 1-4 to use to get the theta we need to add to the list.
        # If D(x) == 0 for real, maybe get that eqn8 working to get the theta.
        # If D(x) < 0, there are three thetas in play - two almost duplicate, and you should already have the one that's close
        # to those, that they are heading towards, from you solve XPrime for theta deal.
        #  You don't want to work with that side of the zero crossing though, since you don't know which of the eqns 5-7 to use
        #  to get the one new theta you want. Actually, It's probably eqn 5 or 6. But it's easier to work on the other side where you know there's just one eqn.

        # so, adapt the loop below
        # find the first x makes D(x) land on the positive side of zero
        # then get the theta for that, from the one eqn form 1-4 that applies.
        # then do again like the loop below, but working with thetas, incrementing thetas, generating x(t)s
        # calcing D(x(t)) to get the complete bracket.
        # Remember this is all to make t(x) eqns work - especially t'(x) - to get the D(x) that happens in those
        # to always come out with the correct sign.


        my $d1 = $D->($x_d0_1);
        my $d2x=nexttoward($x_d0_1,0);
        my $d2 = $D->($d2x);
        while ( $d1 == $d2 ) { $d2x = nexttoward($d2x,0); $d2 = $D->($d2x); }

        my ($dirxDneg,$dirxDpos);

        #warn "        if ($d2 < $d1) [",($D->(nexttoward($x_d0_1,0)) < $d1),"]\n";

        if ($d2 < $d1) {
            $dirxDneg = $self->{minx};
            $dirxDpos = $self->{maxx};
        }
        else {
            $dirxDneg = $self->{maxx};
            $dirxDpos = $self->{minx};
        }

        warn "\nwhilewhatch\n" if $whilewatch;

        # step x so D(x) goes to the first zero or negative result
        while ($D->($x_d0_1) >=  0) {$x_d0_1=nexttoward($x_d0_1,$dirxDneg);}
        my $x_d0_1_Dneg = $x_d0_1;
        warn "ww1\n" if $whilewatch;
        # then step x so D(x) goes to the first positive result
        while ($D->($x_d0_1) <= 0) {$x_d0_1=nexttoward($x_d0_1,$dirxDpos);}
        warn "ww2\n" if $whilewatch;

        # $x_d0_1 is now the first x val that makes D(x) come out positive above zero

        # find the t for that x from the one eqn that applies
        my $t1;
        if ($R->($x_d0_1) >= 0) {
            if ($self->{cubQ}<=0) { $t1 = $eqn_1->($x_d0_1); }
            else                  { $t1 = $eqn_2->($x_d0_1); }
        }
        else { #elsif ($R < 0) {
            if ($self->{cubQ}<=0) { $t1 = $eqn_3->($x_d0_1); }
            else                  { $t1 = $eqn_4->($x_d0_1); }
        }

        # that t may not exactly produce the x we derived it from
        # in order to find a good safe t span end, make sure the t we use produces an x that produces a D(x) with correct sign.
        
        $d1 = $D->($self->bezierEvalXofT($t1));
        my ($dirtDneg,$dirtDpos);
        my $t2 = nexttoward($t1,0);
        my $d2 = $D->($self->bezierEvalXofT(nexttoward($t1,0)));
        while ($d2 == $d1) { $d2 = $D->($self->bezierEvalXofT($t2 = nexttoward($t2,0))); }
        if ($d2 < $d1) {
            $dirtDneg = 0;
            $dirtDpos = 1;
        }
        else {
            $dirtDneg = 1;
            $dirtDpos = 0;
        }

        # step t so D(x(t)) goes to the first zero or negative result
        while ($D->($self->bezierEvalXofT($t1)) >=  0) {$t1=nexttoward($t1,$dirtDneg);}
        warn "ww3\n" if $whilewatch;
        my $t1neg = $t1;
        # then step t so D(x(t)) goes to the first positive result
        while ($D->($self->bezierEvalXofT($t1)) <= 0) {$t1=nexttoward($t1,$dirtDpos);}
        warn "ww4\n" if $whilewatch;

        # $t1 is now the first t val that makes D(x(t)) come out positive above zero

        if (0) {
        warn "bracket1 [$x_d0_1_Dneg,$x_d0_1],[$t1neg,$t1]\n";
        warn "[",$D->($self->bezierEvalXofT($t1neg)),",",$D->($self->bezierEvalXofT($t1)),"]\n";    
        warn "[",$D->($x_d0_1_Dneg),",",$D->($x_d0_1),"]\n";
        warn "$x_d0_1_Dneg vs1\n",$self->bezierEvalXofT($t1neg),"\n";
        warn "$x_d0_1 vs2\n",$self->bezierEvalXofT($t1),"\n";
        }

        push @more_t, [$t1neg,$t1];
        push @more_x, [$self->bezierEvalXofT($t1neg),$self->bezierEvalXofT($t1)];
        #push @more_x, [$x_d0_1_Dneg,$x_d0_1];

    }
    else {
        #warn "no D(x)==0 [1]\n";
    }

    ################### COPY PASTE ADAPT REPEAT OF ABOVE CODE
    # do this better, instead of copy paste adapt
    # make it a loop I guess? or a sub{}
    #

    if ($x_d0_2 && $x_d0_2 !~ /nan/i && $x_d0_2 > $self->{minx} && $x_d0_2 < $self->{maxx}) {

        my $d1 = $D->($x_d0_2);
        my ($dirxDneg,$dirxDpos);
        if ($D->(nexttoward($x_d0_2,0)) < $d1) {
            $dirxDneg = $self->{minx};
            $dirxDpos = $self->{maxx};
        }
        else {
            $dirxDneg = $self->{maxx};
            $dirxDpos = $self->{minx};
        }

        # step x so D(x) goes to the first zero or negative result
        while ($D->($x_d0_2) >=  0) {$x_d0_2=nexttoward($x_d0_2,$dirxDneg);}
        my $x_d0_2_Dneg = $x_d0_2;
        warn "ww5\n" if $whilewatch;

        # then step x so D(x) goes to the first positive result
        while ($D->($x_d0_2) <= 0) {$x_d0_2=nexttoward($x_d0_2,$dirxDpos);}
        warn "ww6\n" if $whilewatch;

        # $x_d0_2 is now the first x val that makes D(x) come out positive above zero

        # find the t for that x from the one eqn that applies
        my $t1;
        if ($R->($x_d0_2) >= 0) {
            if ($self->{cubQ}<=0) { $t1 = $eqn_1->($x_d0_2); }
            else                  { $t1 = $eqn_2->($x_d0_2); }
        }
        else { #elsif ($R < 0) {
            if ($self->{cubQ}<=0) { $t1 = $eqn_3->($x_d0_2); }
            else                  { $t1 = $eqn_4->($x_d0_2); }
        }

        # that t may not exactly produce the x we derived it from
        # in order to find a good safe t span end, make sure the t we use produces an x that produces a D(x) with correct sign.
        
        $d1 = $D->($self->bezierEvalXofT($t1));
        my ($dirtDneg,$dirtDpos);
        my $t2 = nexttoward($t1,0);
        my $d2 = $D->($self->bezierEvalXofT(nexttoward($t1,0)));
        while ($d2 == $d1) { $d2 = $D->($self->bezierEvalXofT($t2 = nexttoward($t2,0))); }
        if ($d2 < $d1) {
            $dirtDneg = 0;
            $dirtDpos = 1;
        }
        else {
            $dirtDneg = 1;
            $dirtDpos = 0;
        }

        # step t so D(x(t)) goes to the first zero or negative result
        while ($D->($self->bezierEvalXofT($t1)) >=  0) {$t1=nexttoward($t1,$dirtDneg);}
        my $t1neg = $t1;
        warn "ww7\n" if $whilewatch;

        # then step t so D(x(t)) goes to the first positive result
        while ($D->($self->bezierEvalXofT($t1)) <= 0) {$t1=nexttoward($t1,$dirtDpos);}
        warn "ww8\n" if $whilewatch;

        # $t1 is now the first t val that makes D(x(t)) come out positive above zero

        if (0) {
        warn "bracket2 [$x_d0_2_Dneg,$x_d0_2],[$t1neg,$t1]\n";
        warn "[",$D->($self->bezierEvalXofT($t1neg)),",",$D->($self->bezierEvalXofT($t1)),"]\n";    
        warn "[",$D->($x_d0_2_Dneg),",",$D->($x_d0_2),"]\n";
        warn "$x_d0_2_Dneg vs1\n",$self->bezierEvalXofT($t1neg),"\n";
        warn "$x_d0_2 vs2\n",$self->bezierEvalXofT($t1),"\n";
        }

        push @more_t, [$t1neg,$t1];
        push @more_x, [$self->bezierEvalXofT($t1neg),$self->bezierEvalXofT($t1)];
        #push @more_x, [$x_d0_2_Dneg,$x_d0_2];

    }
    else {
        #warn "no D(x)==0 [2]\n";
    }
    #
    ################### end COPY PASTE ADAPT REPEAT OF ABOVE CODE

    if ($y_d0_1 && $y_d0_1 !~ /nan/i) {
        # then figure and add the right xs and ts to capture that D(y)==0 crossing
    }
    else {
        #warn "no D(y)==0\n";
    }

    # THIS @interval_xs SETUP IS THE OLD STUFF WE'RE MOVING AWAY FROM
    # HERE JUST FOR REFERENCE. DELETE WHEN DONE WITH OTHER REPLACEMENT STUFF
    my @interval_xs;
    foreach my $t (@all_ts) {
        # bez eval has no t filter, so a t that is out of the 0 to 1 range should still work
        my $x = (ref($t) && !$self->{isLite}) ? $self->bezierEvalXofTBig($t):$self->bezierEvalXofT($t);
        if ( (($t>0 || $t eq 0) && ($t<1 || $t eq 1)) || (($x>$self->{minx} || $x eq $self->{minx}) && ($x<$self->{maxx} || $x eq $self->{maxx})) ) {
            push @interval_xs, $x;
        }
    }


    # TODO: if you have BigFloats at this point in @interval_xs
    #       (probably because some @all_ts were BigFloats)
    #       then lots of the following is probably buggy?
    #       Was getting partly wrong results when I had some BigFloats get in there.
    #       Could be that nexttoward() doesn't work on BigFloats?

    @interval_xs = sort {$a<=>$b} @interval_xs;
        


    my @t_intervals = ([0,undef]);
    my @x_intervals = ([$self->bezierEvalXofT(0),undef]);

    for (my $tind=1;$tind < @all_ts - 1;$tind++) {

        my $newt = $all_ts[$tind];
        my $lastnewt = $newt;
        my $m0 = $self->bezierEvalXPrimeofT($newt);
        my $m1 = $self->bezierEvalXPrimeofT(nexttoward($newt,0));
        my $dir = $m0 > $m1 ? 0 : 1;
        my $e0=abs($m0-$m1);
        my $e1=abs($m0-$m1);

        #warn "\nstart $all_ts[$tind]\nnewt,e1,m1: $newt, $e1, $m1\n";

        my $cnt=0; # debug safety limit

        while (
               $cnt++ < 100 &&        
               $m1 != 0 && 
               ($e1 <= $e0 || $e1 != 0) && 
               $thisdir == $dir
              ) {
            $lastnewt = $newt;
            $newt = nexttoward($newt,$dir);
            $m0 = $m1;
            $m1 = $self->bezierEvalXPrimeofT($newt);
            $e0 = $e1;
            $e1=abs($m0-$m1);
            #warn "X[$cnt $dir] newt,e1,m1: $lastnewt to $newt, $e1, $m1\n";
            $thisdir = $m0 > $m1 ? 0 : 1;
        }
        

        $t_intervals[$tind - 1]->[1] = $dir == 1 ? $lastnewt : $newt;
        $t_intervals[$tind]->[0]     = $dir == 1 ? $newt : $lastnewt;

        $x_intervals[$tind - 1]->[1] = $self->bezierEvalXofT($dir == 1 ? $lastnewt : $newt);
        $x_intervals[$tind]->[0] = $self->bezierEvalXofT($dir == 1 ? $newt : $lastnewt);


        #warn "new t: $newt\n";
        $all_ts[$tind] = $newt;
    }

    $t_intervals[-1]->[1] = 1;
    $x_intervals[-1]->[1] = $self->bezierEvalXofT(1);



    my @m0ts = ((!$self->{isLite})?$self->solveYPrimeforThetaBig(Math::BigFloat->bzero()):$self->solveYPrimeforTheta(0));


    my @m0xs;


    for (my $tind=$#m0ts;$tind > -1;$tind--) {
        my $newt = $m0ts[$tind];
        my $lastnewt = $newt;
        my $m0 = $self->bezierEvalYPrimeofT($newt);
        my $m1 = $self->bezierEvalYPrimeofT(nexttoward($newt,0));
        my $dir = $m0 > $m1 ? 0 : 1;
        my $e0=abs(0);
        my $e1=abs($m0-$m1);
        my $thisdir = $dir;

        #warn "\nstart $m0ts[$tind]\ne1,m1: $e1, $m1, $dir\n";

        my $cnt=0; # debug safety limit
        while (
               $cnt++ < 100 &&        
               $m1 != 0 && 
               ($e1 <= $e0 || $e1 != 0) && 
               $thisdir == $dir
              ) {
            $lastnewt = $newt;
            $newt = nexttoward($newt,$dir);
            $m0 = $m1;
            $m1 = $self->bezierEvalYPrimeofT($newt);
            $e0 = $e1;
            $e1=abs($m0-$m1);
            #warn "Y[$cnt] e1,m1: $e1, $m1, $m0\n";
            $thisdir = $m0 > $m1 ? 0 : 1;
        }

        splice @m0ts, $tind, 1, [$lastnewt,$newt];
        $m0xs[$tind] = [$self->bezierEvalXofT($lastnewt),$self->bezierEvalXofT($newt)];

        #warn "new t: $newt\n";
        $all_ts[$tind] = $newt;
    }

    # split intervals that contain dy/dt = 0 ts
    for (my $ind=$#t_intervals;$ind > -1;$ind--) {
        for (my $m0ts_ind = 0; $m0ts_ind < @m0ts; $m0ts_ind++) {
            if ($m0ts[$m0ts_ind]->[0] > $t_intervals[$ind]->[0] && $m0ts[$m0ts_ind]->[0] < $t_intervals[$ind]->[1]) {
                splice @t_intervals, $ind+1, 0, [$m0ts[$m0ts_ind]->[1],$t_intervals[$ind]->[1]];
                splice @x_intervals, $ind+1, 0, [$m0xs[$m0ts_ind]->[1],$x_intervals[$ind]->[1]];
                $t_intervals[$ind]->[1] = $m0ts[$m0ts_ind]->[0];
                $x_intervals[$ind]->[1] = $m0xs[$m0ts_ind]->[0];
            }
        }

    }

    # split intervals that contain D(x) zero crossings
    for (my $ind=$#t_intervals;$ind > -1;$ind--) {
        for (my $more_t_ind = 0; $more_t_ind < @more_t; $more_t_ind++) {
            if ($more_t[$more_t_ind]->[0] > $t_intervals[$ind]->[0] && $more_t[$more_t_ind]->[0] < $t_intervals[$ind]->[1]) {
                splice @t_intervals, $ind+1, 0, [$more_t[$more_t_ind]->[1],$t_intervals[$ind]->[1]];
                splice @x_intervals, $ind+1, 0, [$more_x[$more_t_ind]->[1],$x_intervals[$ind]->[1]];
                $t_intervals[$ind]->[1] = $more_t[$more_t_ind]->[0];
                $x_intervals[$ind]->[1] = $more_x[$more_t_ind]->[0];
            }
        }

    }

    my @reversed;

    for (my $ind=0;$ind < @x_intervals;$ind++) {
        if ($x_intervals[$ind]->[0] <= $x_intervals[$ind]->[1]) {
            $reversed[$ind] = 0;
        }
        else {
            $x_intervals[$ind] = [$x_intervals[$ind]->[1], $x_intervals[$ind]->[0]];
            $t_intervals[$ind] = [$t_intervals[$ind]->[1], $t_intervals[$ind]->[0]];
            $reversed[$ind] = 1;
        }
    }

    warn join("\n",map {'xintrvl['.join(',',@{$x_intervals[$_]})."]\n".
                        'tintrvl['.join(',',@{$t_intervals[$_]})."]"
               } (0 .. $#x_intervals)),"\n"
    if (0);


    my %equation_map;

    for (my $i=0; $i < @x_intervals; $i++) {
        my $xa = $x_intervals[$i]->[0];
        my $xb = $x_intervals[$i]->[1];

        my $xmid = ($xa + $xb) / 2.0;
        my $tmid = ($t_intervals[$i]->[0] + $x_intervals[$i]->[1]) / 2.0;
        my $ymid = $self->bezierEvalYofT($tmid);

        my $Rmid = $R->($xmid);
        my $Dmid = $D->($xmid);

        my $Rmid_ofy = $R_ofy->($ymid);
        my $Dmid_ofy = $D_ofy->($ymid);

        # first figure the new _ofy eqs, so we can stuff them into lut entries when it's time, after the eqns for x

        my @yeqns;

        if ($Dmid_ofy > 0 || $Dmid_ofy eq 0) {

            if ($Dmid_ofy eq 0) {
                die "D == 0 not handled yet (for y)";
            }

            if ($Rmid_ofy >= 0) {
                if ($self->{cubQ_ofy}<=0) {    # use equation 1
                    push @yeqns, $eqn_1_ofy,$eqn_1to4_prime_ofy,$eqn_1to4_2prime_ofy;
                }
                else {                     # use equation 2
                    push @yeqns, $eqn_2_ofy,$eqn_1to4_prime_ofy,$eqn_1to4_2prime_ofy;
                }
            }
            elsif ($Rmid_ofy < 0) {
                if ($self->{cubQ_ofy}<=0) {    # use equation 3
                    push @yeqns, $eqn_3_ofy,$eqn_1to4_prime_ofy,$eqn_1to4_2prime_ofy;
                }
                else {                     # use equation 4
                    push @yeqns, $eqn_4_ofy,$eqn_1to4_prime_ofy,$eqn_1to4_2prime_ofy;
                }
            }

        }
        elsif ($Dmid_ofy<0) {
            foreach my $eqns ([$eqn_5_ofy,$eqn_5_prime_ofy,$eqn_5_2prime_ofy],
                              [$eqn_6_ofy,$eqn_6_prime_ofy,$eqn_6_2prime_ofy],
                              [$eqn_7_ofy,$eqn_7_prime_ofy,$eqn_7_2prime_ofy] ) {
                my $eqn = $eqns->[0];
                my ($tlow,$thigh) = ($t_intervals[$i]->[0]<$t_intervals[$i]->[1]) ? ($t_intervals[$i]->[0],$t_intervals[$i]->[1]) : ($t_intervals[$i]->[1],$t_intervals[$i]->[0]);
                if ($tmid > $tlow && $tmid < $thigh) {
                    push @yeqns, @$eqns;
                }
            }
        }
        else {die "Dmid_ofy not < 0 and Dmid_ofy not >= 0 ???"}

        # now set up the lut entries that were originally more x oriented
        # and stuff those y eqs in ther now
        if ($Dmid > 0 || $Dmid eq 0) {
            if ($Dmid eq 0) {
                die "D == 0 not handled yet";
            }

            if ($Rmid >= 0) {
                if ($self->{cubQ}<=0) {    # use equation 1
                    push @{$equation_map{$eqn_1}}, [[$eqn_1,$eqn_1to4_prime,$eqn_1to4_2prime,@yeqns],$x_intervals[$i],$t_intervals[$i],$reversed[$i]];
                }
                else {                     # use equation 2
                    push @{$equation_map{$eqn_2}}, [[$eqn_2,$eqn_1to4_prime,$eqn_1to4_2prime,@yeqns],$x_intervals[$i],$t_intervals[$i],$reversed[$i]];
                }
            }
            elsif ($Rmid < 0) {
                if ($self->{cubQ}<=0) {    # use equation 3
                    push @{$equation_map{$eqn_3}}, [[$eqn_3,$eqn_1to4_prime,$eqn_1to4_2prime,@yeqns],$x_intervals[$i],$t_intervals[$i],$reversed[$i]];
                }
                else {                     # use equation 4
                    push @{$equation_map{$eqn_4}}, [[$eqn_4,$eqn_1to4_prime,$eqn_1to4_2prime,@yeqns],$x_intervals[$i],$t_intervals[$i],$reversed[$i]];
                }
            }

        }
        elsif ($Dmid<0) {
            foreach my $eqns ([$eqn_5,$eqn_5_prime,$eqn_5_2prime,@yeqns],
                              [$eqn_6,$eqn_6_prime,$eqn_6_2prime,@yeqns],
                              [$eqn_7,$eqn_7_prime,$eqn_7_2prime,@yeqns] ) {
                my $eqn = $eqns->[0];
                
                my $tmid = ($eqn->($xa) + $eqn->($xb))/2;

                my ($tlow,$thigh) = ($t_intervals[$i]->[0]<$t_intervals[$i]->[1]) ? ($t_intervals[$i]->[0],$t_intervals[$i]->[1]) : ($t_intervals[$i]->[1],$t_intervals[$i]->[0]);

                if ($tmid > $tlow && $tmid < $thigh) {
                    if (!exists $equation_map{$eqn}) { $equation_map{$eqn} = []; }
                    push @{$equation_map{($eqn)}}, [$eqns,$x_intervals[$i],$t_intervals[$i],$reversed[$i]];
                    @{$equation_map{($eqn)}} = sort {$a->[1]->[0]<=>$b->[1]->[0]} @{$equation_map{($eqn)}};
                }
            }
        }
        else {die "Dmid not < 0 and Dmid not >= 0 ???"}

    }

    # make eqn LUT
    my @XtoTLUT;
    foreach my $xspan0 (map {@$_} values %equation_map) {

        push @XtoTLUT, $xspan0;

    }

    # sort in ascending t order
    @XtoTLUT = sort {$a->[2]->[$a->[3]?0:1] <=> $b->[2]->[$b->[3]?0:1]} @XtoTLUT;

    $self->{XtoTLUT} = \@XtoTLUT;

    #warn(join(',',(map {'['.join(',',map {ref($_)?$eqn_names{$_}:$_} @$_).']'} @{$_}[0,1,2]),$_->[3]),"\n") for @XtoTLUT;

    return $self;

    }


sub get_xspans_from_x_range {
    my ($self,$xa,$xb) = @_;
    die "2nd xrange arg s/b >= to first [ $xb < $xa ]" if ($xb<$xa);
    my @xspans_a = grep {($xa > $_->[1]->[0] || $xa eq $_->[1]->[0]) && ($xa < $_->[1]->[-1]                        )} @{$self->{XtoTLUT}};
    my @xspans_b = grep {($xb > $_->[1]->[0]                       ) && ($xb < $_->[1]->[-1] || $xb eq $_->[1]->[-1])} @{$self->{XtoTLUT}};
    return (@xspans_a, @xspans_b);
}

# these are all t(x)
# ideally you should also have t(y) available, but you didn't think through that for when setting this up
#
# As I ttryto start doing offset points, I see you always need the flag for if t is increasing or decreasing over the corresponding x ranges
# That flag is the 4th item in XtoTLUT entries. Need it to know if "left" offset is a positive or negative 90 degree rotation from the tangent.
# That's why we're returning $_->[3] in most of these t*() functions.
# And that's why t_prime() is returning array refs, instead of just the t_prime vals. Which is irksome.
# There's probably a next version of all this needed soon, to do t(x) AND t(y) and be aware of that flag, etc.
#
# You might do a special BezierCubicSegmentMonotonic or something like that,
# where you've broken up a normal Bez into new Beziers that only have a single t(x) in effect over their x or y ranges.
# Or where there's only one each of t(x) and t(y) over the range (that's even more broken up? isn't it? or same? not sure at moment.)

sub t {
    #warn "t(x)\n";
    my $self = shift;
    my $x = shift;
    my @xspans = grep {($x > $_->[1]->[0] || $x eq $_->[1]->[0]) && ($x < $_->[1]->[-1] || $x eq $_->[1]->[-1])} @{$self->{XtoTLUT}};
    my @ret = map {$_->[0]->[0]->($x)} @xspans;
    return wantarray ? @ret : $ret[0];
}
sub t_prime {
    my $self = shift;
    my $x = shift;
    my @xspans = grep {($x > $_->[1]->[0] || $x eq $_->[1]->[0]) && ($x < $_->[1]->[-1] || $x eq $_->[1]->[-1])} @{$self->{XtoTLUT}};
    my @ret = map {[$_->[0]->[1]->($x),$_->[3]]} @xspans;
    return wantarray ? @ret : $ret[0];
}
sub t_2prime {
    my $self = shift;
    my $x = shift;
    my @xspans = grep {($x > $_->[1]->[0] || $x eq $_->[1]->[0]) && ($x < $_->[1]->[-1] || $x eq $_->[1]->[-1])} @{$self->{XtoTLUT}};
    my @ret = map {[$_->[0]->[2]->($x),$_->[3]]} @xspans;
    return wantarray ? @ret : $ret[0];
}
sub t_t_prime {
    my $self = shift;
    my $x = shift;
    my @xspans = grep {($x > $_->[1]->[0] || $x eq $_->[1]->[0]) && ($x < $_->[1]->[-1] || $x eq $_->[1]->[-1])} @{$self->{XtoTLUT}};
    my @ret = map {[$_->[0]->[0]->($x),$_->[0]->[1]->($x),$_->[3]]} @xspans;
    return wantarray ? @ret : $ret[0];
}
sub t_t_prime_of_y {
    my ($self, $y, $t)=@_;
    $y //= $self->bezierEvalYofT($t);
    my @xspans = grep {($x > $_->[1]->[0] || $x eq $_->[1]->[0]) && ($x < $_->[1]->[-1] || $x eq $_->[1]->[-1])} @{$self->{XtoTLUT}};
    if (defined $t) {
        @xspans = grep {
            my ($lt,$ht)=$_->[2]->[0] < $_->[2]->[-1]?($_->[2]->[0] < $_->[2]->[-1]):($_->[2]->[-1] < $_->[2]->[0]);
            $t >= $lt && $t <= $ht;
        } @xspans;
    }
    my @ret = map {[$_->[0]->[3]->($y),$_->[0]->[4]->($y),$_->[3]]} @xspans;
    return wantarray ? @ret : $ret[0];
}
sub t_t_prime_t_2prime {
    my $self = shift;
    my $x = shift;
    my @xspans = grep {($x > $_->[1]->[0] || $x eq $_->[1]->[0]) && ($x > $_->[1]->[-1] || $x eq $_->[1]->[-1])} @{$self->{XtoTLUT}};
    my @ret = map {[$_->[0]->[0]->($x),$_->[0]->[1]->($x),$_->[0]->[2]->($x),$_->[3]]} @xspans;
    return wantarray ? @ret : $ret[0];
}
sub t_t_prime_t_2prime_of_y {
    my ($self, $y, $t)=@_;
    $y //= $self->bezierEvalYofT($t);
    my @xspans = grep {($x > $_->[1]->[0] || $x eq $_->[1]->[0]) && ($x < $_->[1]->[-1] || $x eq $_->[1]->[-1])} @{$self->{XtoTLUT}};
    if (defined $t) {
        @xspans = grep {
            my ($lt,$ht)=$_->[2]->[0] < $_->[2]->[-1]?($_->[2]->[0] < $_->[2]->[-1]):($_->[2]->[-1] < $_->[2]->[0]);
            $t >= $lt && $t <= $ht;
        } @xspans;
    }
    my @ret = map {[$_->[0]->[3]->($y),$_->[0]->[4]->($y),$_->[0]->[5]->($y),$_->[3]]} @xspans;
    return wantarray ? @ret : $ret[0];
}

sub f {
    my $self = shift;
    my $x = shift;
    my @t = $self->t($x);
    my @ret;
    foreach my $t (@t) {
        my $y = $self->bezierEvalYofT($t);
        push @ret, $y;
    }
    if (@ret > 0) {
        return wantarray ? @ret : $ret[0];
    }
    else {
        return;
    }
}

sub f_prime {
    my $self = shift;
    my $x = shift;
    my @t_t_prime = $self->t_t_prime($x);
    my @ret;
    foreach my $t_t_prime (@t_t_prime) {
        my ($t,$t_prime) = @$t_t_prime;
        my $y_prime_of_x = $self->{Em3} * $t**2 * $t_prime  +  $self->{Fm2} * $t * $t_prime + $self->{G} * $t_prime;
        push @ret, $y_prime_of_x;
    }
    if (@ret>0) {
        return wantarray ? @ret : $ret[0];
    }
    else {
        return;
    }
}

sub F_prime {
    my ($self,$y,$t) = @_;
    my @t_t_prime_of_y = $self->t_t_prime_of_y($y,$t);
    $y //= $self->bezierEvalYofT($t);
    my @ret;
    foreach my $t_t_prime_of_y (@t_t_prime_of_y) {
        my ($t_of_y,$t_prime_of_y) = @$t_t_prime_of_y;
        my $x_prime_of_y = $self->{Am3} * $t_of_y**2 * $t_prime_of_y  +  $self->{Bm2} * $t_of_y * $t_prime_of_y + $self->{C} * $t_prime_of_y;
        push @ret, $x_prime_of_y;
    }
    if (@ret>0) {
        return wantarray ? @ret : $ret[0];
    }
    else {
        return;
    }
}

sub f_2prime {
    my $self = shift;
    my $x = shift;
    my @t_t_prime_t_2prime = $self->t_t_prime_t_2prime($x);
    my @ret;
    foreach my $t_t_prime_t_2prime (@t_t_prime_t_2prime) {
        my ($t,$t_prime,$t_2prime) = @$t_t_prime_t_2prime;
        my $y_2prime_of_x = $self->{Em3} * ($t**2 * $t_2prime + 2*$t*$t_prime**2)
                          + $self->{Fm2} * ($t    * $t_2prime +      $t_prime**2)
                          + $self->{G}   * $t_2prime;
        push @ret, $y_2prime_of_x;
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

    # need to figure single return value version when $t is provided
    if (!defined $y && defined $t) {$y=$self->bezierEvalYofT($t);}

    my @t_t_prime_t_2prime_of_y = $self->t_t_prime_t_2prime_of_y($y);
    my @ret;
    foreach my $t_t_prime_t_2prime_of_y (@t_t_prime_t_2prime_of_y) {
        my ($t,$t_prime,$t_2prime) = @$t_t_prime_t_2prime_of_y;
        my $x_2prime_of_y = $self->{Am3} * ($t**2 * $t_2prime + 2*$t*$t_prime**2)
                          + $self->{Bm2} * ($t    * $t_2prime +      $t_prime**2)
                          + $self->{C}   * $t_2prime;
        push @ret, $x_2prime_of_y;
    }
    if (@ret>0) {
        return wantarray ? @ret : $ret[0];
    }
    else {
        return;
    }
}


# combo of above two
sub f_prime_f_2prime {
    my $self = shift;
    my $x = shift;
    my @t_t_prime_t_2prime = $self->t_t_prime_t_2prime($x);
    my @ret;
    foreach my $t_t_prime_t_2prime (@t_t_prime_t_2prime) {
        my ($t,$t_prime,$t_2prime) = @$t_t_prime_t_2prime;
        my $y_prime_of_x  = $self->{Em3} * $t**2 * $t_prime  +  $self->{Fm2} * $t * $t_prime + $self->{G} * $t_prime;
        my $y_2prime_of_x = $self->{Em3} * ($t**2 * $t_2prime + 2*$t*$t_prime**2)
                          + $self->{Fm2} * ($t    * $t_2prime +      $t_prime**2)
                          + $self->{G}   * $t_2prime;
        push @ret, [$y_prime_of_x, $y_2prime_of_x];
    }
    if (@ret>0) {
        return wantarray ? @ret : $ret[0];
    }
    else {
        return;
    }
}

sub curvature {
    my $self = shift;
    my $x = shift;
    my @f_prime_f_2prime = $self->f_prime_f_2prime($x);
    my @ret;
    foreach my $f_prime_f_2prime (@f_prime_f_2prime) {
        my ($yp,$ypp) = @$f_prime_f_2prime;
        my $k = $ypp / (sqrt(1 + $yp**2)**3);
        push @ret, $k;
    }
    if (@ret>0) {
        return wantarray ? @ret : $ret[0];
    }
    else {
        return;
    }
}

sub radius_of_curvature {
    my $self = shift;
    my $x = shift;
    my @ret = map {1/$_} $self->curvature($x);
    if (@ret>0) {
        return wantarray ? @ret : $ret[0];
    }
    else {
        return;
    }
}

# hmm, does this make sense?
sub f_offset {

    my $self = shift;
    my $x = shift;
    my $distance = shift;
    my @t_t_prime = $self->t_t_prime($x);
    my @ret;
    foreach my $t_t_prime (@t_t_prime) {
      my ($t,$t_prime,$reversed) = @$t_t_prime;
      my $y = $self->bezierEvalYofT($t);
      my $y_prime_of_x = $self->{Em3} * $t**2 * $t_prime  +  $self->{Fm2} * $t * $t_prime + $self->{G} * $t_prime;
      my $a = atan2($y_prime_of_x,1) - $pi/2 * ($reversed?-1:1);
      push @ret, [$x + $distance * cos($a), $y + $distance * sin($a)];
    }
    if (@ret>0) {
      return wantarray ? @ret : $ret[0];
    }
    else {
      return;
    }

}

sub get_tspan {
    my ($self, $t) = @_;
    my @tspans = grep {
        my ($low_t,$high_t) = ($_->[2]->[0] < $_->[2]->[-1] ? ($_->[2]->[0],$_->[2]->[-1]) : ($_->[2]->[-1],$_->[2]->[0]) );
        #warn "t situation: ($t > $low_t || $t eq $low_t) && ($t < $high_t || $t eq $high_t)\n";
        #warn " 1: ",($t > $low_t),"\n";
        #warn " 2: ",($t eq $low_t),"\n";
        #warn " 3: ",($t < $high_t),"\n";
        #warn " 4: ",($t eq $high_t),"\n";
        ($t > $low_t || $t eq $low_t) && ($t < $high_t || $t eq $high_t)
    } @{$self->{XtoTLUT}} ;

    return $tspans[0];
}

sub X_offset{
    my ($self,$t,$offset,$tprimefunc,$reversed) = @_;
    if (!defined($tprimefunc) || !defined($reversed)) {
        my $tspan = $self->get_tspan($t);
        $tprimefunc //= $tspan->[0]->[1];
        $reversed //= $tspan->[3];
    }
    my $x = $self->bezierEvalXofT($t);
    my $tp = $tprimefunc->($x);
    my $y_prime_of_x = $tp == Inf ? Inf : $self->{Em3} * $t**2 * $tp  +  $self->{Fm2} * $t * $tp + $self->{G} * $tp;
    my $a = atan2($y_prime_of_x,1) - $pi/2 * ($reversed?-1:1);
    my $retx = $x + $offset * cos($a);
    return $retx;
}
sub Y_offset{
    my ($self,$t,$offset,$tprimefunc,$reversed) = @_;
    if (!defined($tprimefunc) || !defined($reversed)) {
        my $tspan = $self->get_tspan($t);
        $tprimefunc //= $tspan->[0]->[1];
        $reversed //= $tspan->[3];
    }
    my $x = $self->bezierEvalXofT($t);
    my $y = $self->bezierEvalYofT($t);
    my $tp = $tprimefunc->($x);
    my $y_prime_of_x = $tp == Inf ? Inf : $self->{Em3} * $t**2 * $tp  +  $self->{Fm2} * $t * $tp + $self->{G} * $tp;
    my $a = atan2($y_prime_of_x,1) - $pi/2 * ($reversed?-1:1);
    my $rety = $y + $offset * sin($a);
    return $rety;
}

sub point_offset {
    my ($self, $t, $distance) = @_;
    my $x = $self->bezierEvalXofT($t);
    my $y = $self->bezierEvalYofT($t);
    # only want the first lut entry we get matching on t span. (Usually there's only one, but there could be two if $t is right on a t span boundary.)
    #warn "grep ts\n";
    my @tspans = grep {
        my ($low_t,$high_t) = ($_->[2]->[0] < $_->[2]->[-1] ? ($_->[2]->[0],$_->[2]->[-1]) : ($_->[2]->[-1],$_->[2]->[0]) );
        #warn "t situation: ($t > $low_t || $t eq $low_t) && ($t < $high_t || $t eq $high_t)\n";
        #warn " 1: ",($t > $low_t),"\n";
        #warn " 2: ",($t eq $low_t),"\n";
        #warn " 3: ",($t < $high_t),"\n";
        #warn " 4: ",($t eq $high_t),"\n";
        ($t > $low_t || $t eq $low_t) && ($t < $high_t || $t eq $high_t)
        } @{$self->{XtoTLUT}} ;
    my $tspan = $tspans[0];
    my $t_prime = $tspan->[0]->[1]->($x);
    my $reversed = $tspan->[3];

    my $y_prime_of_x = $self->{Em3} * $t**2 * $t_prime  +  $self->{Fm2} * $t * $t_prime + $self->{G} * $t_prime;
    my $a = atan2($y_prime_of_x,1) + $pi/2 * ($reversed?-1:1);

    my $ret=[$x + $distance * cos($a), $y + $distance * sin($a)];

    return $ret;

}

# belongs with offset related functions
# used in Intersections.pm
sub t_from_xoff {
    my ($self, $xoff, $offset, $tbounds, $tprimefunc, $reversed) = @_;

    my $find_matching_offset_x = sub {
        #warn "findloop\n";
        my $t0 = $_[0];
        my $x0 = $self->bezierEvalXofT($t0);
        my $t_prime_0 = $tprimefunc->($x0);
        my $y_prime_of_x_0 = $t_prime_0 == Inf ? Inf : $self->{Em3} * $t0**2 * $t_prime_0  +  $self->{Fm2} * $t0 * $t_prime_0 + $self->{G} * $t_prime_0;
        my $a0 = atan2($y_prime_of_x_0,1) - $pi/2 * ($reversed?-1:1);

        #warn "findloop: ", ($xoff - ($x0 + $offset * cos($a0))) , " = $xoff - ($x0 + $offset * cos($a0))\n";

        return $xoff - ($x0 + $offset * cos($a0));
    };

    $bounds = [$tbounds->[0] < $tbounds->[-1] ? ($tbounds->[0], $tbounds->[-1]) : ($tbounds->[-1], $tbounds->[0])];
#    $bounds->[0] = 0 if $bounds->[0] < 0;
#    $bounds->[1] = 1 if $bounds->[1] > 1;

#warn("[$_] ",$find_matching_offset_x->(($bounds->[0]+0.9) + (($bounds->[1]-($bounds->[0]+0.9))*($_/300))),"\n") for (0..300);

    my ($t, $msg) = FalsePosition($find_matching_offset_x,$bounds,0.000001,($tbounds->[0]+$tbounds->[-1])/2,'t_from_xoff');

    if ($msg) {warn "false position message: $msg";}    

    $t = 0 if $t < 0;
    $t = 1 if $t > 1;

    return $t;
}

# belongs with offset related functions
sub x_from_xoff {
    my ($self, $xoff, $offset, $tbounds) = @_;
    my $t = $self->t_from_xoff($xoff, $offset, $tbounds);
    my $x = $self->bezierEvalXofT($t);
    return wantarray ? ($x, $t) : $x;
}

# lazy creation of some BigFloats
# keep these around to reuse maybe above when you start to hit precision problems up there
# (but 'till then, these not used, because we're going to comment out or delete the old code below that used them)
sub cubtwotimessqrtofnegQ_big {
    return $self->{cubtwotimessqrtofnegQ_big} if defined($self->{cubtwotimessqrtofnegQ_big});
    $self->{cubtwotimessqrtofnegQ_big} = ($self->cubQ_big())->copy()->bmul(-1.0);
    $self->{cubtwotimessqrtofnegQ_big} = bigsqrt($self->{cubtwotimessqrtofnegQ_big});
    $self->{cubtwotimessqrtofnegQ_big}->bmul(2.0);
    return $self->{cubtwotimessqrtofnegQ_big};
    }
sub cubQ_big {
    return $self->{cubQ_big} if defined($self->{cubQ_big});
    $self->{cubQ_big} = Math::BigFloat->new($self->{cubQ});
    return $self->{cubQ_big};
    }
sub cubsqrtofnegQcubed_big {
    return $self->{cubsqrtofnegQcubed_big} if defined($self->{cubsqrtofnegQcubed_big});
    my $presqrtofnegQcubed = ($self->cubQ_big())->copy()->bpow(3)->bmul(-1.0);
    $self->{cubsqrtofnegQcubed_big} = bigsqrt($presqrtofnegQcubed);
    return $self->{cubsqrtofnegQcubed_big};
    }

sub initDangerRanges {
    my $self=shift;


#redo these in all big - this might muck stuff up, so might delete this section later
    my @extremexs_is=(0,($self->solveXPrimeforThetaBig(Math::BigFloat->bzero())),1);
    #my @extremexs_is=(0,($self->solveXPrimeforThetaBig(Math::BigFloat->bzero())),Math::BigFloat->bone());
    my @extremexs = map {ref($_)?$self->bezierEvalXofTBig($_):$self->bezierEvalXofT($_)} @extremexs_is;
    my @extremexs_sorted = sort {$a<=>$b} @extremexs;

    #my @extremeys_is=(0,($self->solveYPrimeforThetaBig(Math::BigFloat->bzero())),1);
    my @extremeys_is=(Math::BigFloat->bzero(),($self->solveYPrimeforThetaBig(Math::BigFloat->bzero())),Math::BigFloat->bone());
    #my @extremeys_is=(0,($self->solveYPrimeforThetaBig(Math::BigFloat->bzero())),Math::BigFloat->bone());
    my @extremeys = map {ref($_)?$self->bezierEvalYofTBig($_):$self->bezierEvalYofT($_)} @extremeys_is;
    my @extremeys_sorted = sort {$a<=>$b} @extremeys;

    #print "extremeys:\n  ",join("\n  ",@extremeys),"\n";

    $self->{extremexs_is}=[(@extremexs_is)];
    $self->{extremeys_is}=[(@extremeys_is)];
    $self->{extremexs}=[(@extremexs)];
    $self->{extremeys}=[(@extremeys)];

    # Making these for getFeet(), since that gets called so much and shouldn't
    # have to do the expensive bstr() every time.
    my %eisnb_dups;
    $self->{extreme_is_sorted_notbig} = [(sort {$a<=>$b} grep {!$eisnb_dups{$_}++} ((map {ref($_)?$_->bstr():$_} @{$self->{extremexs_is}}),(map {ref($_)?$_->bstr():$_} @{$self->{extremeys_is}})))];
    $self->{tangents_at_sorted_extremes} = [(map {ref($_)?$_->bstr():$_} map {$self->slopeTangent_byTheta($_)} ( @{$self->{extreme_is_sorted_notbig}} ))];

    $self->{minx} = $extremexs_sorted[0];
    $self->{maxx} = $extremexs_sorted[$#extremexs];
    $self->{miny} = $extremeys_sorted[0];
    $self->{maxy} = $extremeys_sorted[$#extremeys];

    #profiled slice making task with DProfLB and found that
    #20% to 60% of time is taking doing inRange checks for beziers
    #probably because those checks are using these min/max x/y vals
    #that are BigFloats. So let's make regular PERL scalar number versions
    #and use those in the inRange tests. Will probably reintroduce
    #some boundary problems, but it might also give a 50%+ speedup.

# SPEED ATTEMPTS
#    $self->{minx_notbig} = ref($self->{minx}) ? eval(sprintf("%.20f",$self->{minx}->bstr())):$self->{minx};
#    $self->{maxx_notbig} = ref($self->{maxx}) ? eval(sprintf("%.20f",$self->{maxx}->bstr())):$self->{maxx};
#    $self->{miny_notbig} = ref($self->{miny}) ? eval(sprintf("%.20f",$self->{miny}->bstr())):$self->{miny};
#    $self->{maxy_notbig} = ref($self->{maxy}) ? eval(sprintf("%.20f",$self->{maxy}->bstr())):$self->{maxy};
    $self->{minx_notbig} = ref($self->{minx}) ? 0 + sprintf("%.20f",$self->{minx}->bstr()):$self->{minx};
    $self->{maxx_notbig} = ref($self->{maxx}) ? 0 + sprintf("%.20f",$self->{maxx}->bstr()):$self->{maxx};
    $self->{miny_notbig} = ref($self->{miny}) ? 0 + sprintf("%.20f",$self->{miny}->bstr()):$self->{miny};
    $self->{maxy_notbig} = ref($self->{maxy}) ? 0 + sprintf("%.20f",$self->{maxy}->bstr()):$self->{maxy};

    #print "minx_notbig: $self->{minx_notbig}\n";
    #print "maxx_notbig: $self->{maxx_notbig}\n";
     #print "miny_notbig: $self->{miny_notbig}\n";
    #print "maxy_notbig: $self->{maxy_notbig}\n";
#end redo

    #my $perlprecision = Math::BigFloat->new('0.00000000000001');
    my $perlprecision =                       0.00000000000001;
    my $maxdim=(sort {$b<=>$a} map {abs($_)} ($self->{maxx},$self->{maxy},$self->{minx},$self->{miny}))[0];
    $maxdim=~/^([0-9]+)\./;
    #$self->{curveprecision}=$perlprecision * Math::BigFloat->new(''.(10**(length($1)))) ;
    $self->{curveprecision}=$perlprecision * (10**(length($1))) ;
    #print "curve precision: ",$self->{curveprecision},"\n";
    #was doing ...PrimeforThetaBig() for all these, but that takes about 1.5 seconds rather than something like 0.0003

#print "back to precision stuff: infinite slope range is greater than : self->{curveprecision}/perlprecision = $self->{curveprecision}/$perlprecision = ",($self->{curveprecision}/$perlprecision),"\n";
#print "in getting xdangeris:\n";
#print "        ",join(",",$self->solveXPrimeforTheta($perlprecision/$self->{curveprecision}))  ," # near zero slope, positive\n";
#print "        ",join(",",$self->solveXPrimeforTheta(0)),"                  # zero slope\n";
#print "        ",join(",",$self->solveXPrimeforTheta(-$perlprecision/$self->{curveprecision})) ," # near zero slope, negative\n";
#print "        ",join(",",$self->solveXPrimeforTheta($self->{curveprecision}/$perlprecision)) ,"  # toward infinite slope, positive\n";
#print "        ",join(",",$self->solveXPrimeforTheta(-$self->{curveprecision}/$perlprecision)) ," # toward infinite slope, negative\n";

    my @xdangeris = sort {$a<=>$b} (
        ($self->solveXPrimeforTheta($perlprecision/$self->{curveprecision}))  , # near zero slope, positive
       #($self->solveXPrimeforTheta(Math::BigFloat->bzero())),                  # zero slope
        ($self->solveXPrimeforTheta(0)),                  # zero slope
        ($self->solveXPrimeforTheta(-$perlprecision/$self->{curveprecision})) , # near zero slope, negative
        ($self->solveXPrimeforTheta($self->{curveprecision}/$perlprecision)) ,  # toward infinite slope, positive
        ($self->solveXPrimeforTheta(-$self->{curveprecision}/$perlprecision)) , # toward infinite slope, negative
        );
    my @ydangeris = sort {$a<=>$b} (
        ($self->solveYPrimeforTheta($perlprecision/$self->{curveprecision}))  ,
       #($self->solveYPrimeforTheta(Math::BigFloat->bzero())) ,
        ($self->solveYPrimeforTheta(0)) ,
        ($self->solveYPrimeforTheta(-$perlprecision/$self->{curveprecision})) ,
        ($self->solveYPrimeforTheta($self->{curveprecision}/$perlprecision)) ,
        ($self->solveYPrimeforTheta(-$self->{curveprecision}/$perlprecision)) ,
        );

    my @mxidangerranges;
    for (my $i=0;$i<@xdangeris;$i++) {
        my $p=$self->bezierEvalXPrimeofT($xdangeris[$i]);
        my $pp=$self->bezierEvalXDoublePrimeofT($xdangeris[$i]);
        if (($p < 0 && $pp < 0) || (($p > 0 || $p eq 0) && ($pp > 0 || $pp eq 0))) {
            #is end of range
            if ($i eq 0)          {push(@mxidangerranges,[0,(sort {$a<=>$b} (1,($xdangeris[$i] + $self->{curveprecision})))[0]]);}
            else                  {push(@mxidangerranges,[(sort {$b<=>$a} (0,($xdangeris[$i-1] - $self->{curveprecision})))[0],(sort {$a<=>$b} (1,($xdangeris[$i] + $self->{curveprecision})))[0]]);}
            }
        elsif ($i eq $#xdangeris) {push(@mxidangerranges,[(sort {$b<=>$a} (0,($xdangeris[$i] - $self->{curveprecision})))[0],1]);}
        }
    $self->{mxidangerranges} = \@mxidangerranges;

    $self->{xofidangerranges} = [];
    push @{$self->{xofidangerranges}}, [$self->{p1}->[0],$self->{p1}->[0]]; #zero-length "range" to
    push @{$self->{xofidangerranges}}, [$self->{p2}->[0],$self->{p2}->[0]]; #try to be more exact at endpoints

    foreach (@mxidangerranges) {
        #$self->{xofidangerranges}->[scalar(@{$self->{xofidangerranges}})]=[eval($self->bezierEvalXofTBig($_->[0])->bstr()),eval($self->bezierEvalXofTBig($_->[1])->bstr())];
        push @{$self->{xofidangerranges}}, [$self->bezierEvalXofT($_->[0]),$self->bezierEvalXofT($_->[1])];
        }
    #print "danger xi ranges: ", join('  ',map {'['.$_->[0].','.$_->[1].']'} @mxidangerranges),' from xis: (',join(',',@xdangeris),')' ,"\n";
    #print "danger xofi ranges: ", join('  ',map {'['.$_->[0].','.$_->[1].']'} @{$self->{xofidangerranges}}) , "\n";

    my @myidangerranges;
    $self->{yofidangerranges} = [];
    $self->{yofidangerranges}->[scalar(@{$self->{yofidangerranges}})]=[$self->{p1}->[1],$self->{p1}->[1]];
    $self->{yofidangerranges}->[scalar(@{$self->{yofidangerranges}})]=[$self->{p2}->[1],$self->{p2}->[1]]; #try to be more exact at endpoints
    for (my $i=0;$i<@ydangeris;$i++) {
        my $p=$self->bezierEvalYPrimeofT($ydangeris[$i]);
        my $pp=$self->bezierEvalYDoublePrimeofT($ydangeris[$i]);
        if (($p < 0 && $pp < 0) || (($p > 0 || $p eq 0) && ($pp > 0 || $pp eq 0))) {
            #is end of range
            if ($i == 0)          {push(@myidangerranges,[0,(sort {$a<=>$b} (1,$ydangeris[$i] + $self->{curveprecision}))[0]]);}
            else                  {push(@myidangerranges,[(sort {$b<=>$a} (0,$ydangeris[$i-1] - $self->{curveprecision}))[0],(sort {$a<=>$b} (1,$ydangeris[$i] + $self->{curveprecision}))[0]]);}
            }
        elsif ($i == $#ydangeris) {push(@myidangerranges,[(sort {$b<=>$a} (0,$ydangeris[$i] - $self->{curveprecision}))[0],1]);}
        }
    $self->{myidangerranges} = \@myidangerranges;
    foreach (@myidangerranges) {
        #$self->{yofidangerranges}->[scalar(@{$self->{yofidangerranges}})]=[eval($self->bezierEvalYofTBig($_->[0])->bstr()),eval($self->bezierEvalYofTBig($_->[1])->bstr())];
        # for f(x) you're using non-Big versions in this case. Want to do that here?
        $self->{yofidangerranges}->[scalar(@{$self->{yofidangerranges}})]=[$self->bezierEvalYofTBig($_->[0]),$self->bezierEvalYofTBig($_->[1])];
        }
    #print "danger yi ranges: ", join('  ',map {'['.$_->[0].','.$_->[1].']'} @myidangerranges),' from yis: (',join(',',@ydangeris),')' ,"\n";
    #print "danger yofi ranges: ", join('  ',map {'['.$_->[0].','.$_->[1].']'} @{$self->{yofidangerranges}}) , "\n";
    }
sub initBigs_old {
    my $self=shift;
    $self->{A_Big}  = new Math::BigFloat '' . ($self->{p2}->[0] - 3 * $self->{cp2}->[0] + 3 * $self->{cp1}->[0] -     $self->{p1}->[0]);
    $self->{B_Big}  = new Math::BigFloat '' . (                   3 * $self->{cp2}->[0] - 6 * $self->{cp1}->[0] + 3 * $self->{p1}->[0]);
    $self->{C_Big}  = new Math::BigFloat '' . (                                           3 * $self->{cp1}->[0] - 3 * $self->{p1}->[0]);
    $self->{D_Big}  = new Math::BigFloat '' . (                                                                       $self->{p1}->[0]);
    $self->{E_Big}  = new Math::BigFloat '' . ($self->{p2}->[1] - 3 * $self->{cp2}->[1] + 3 * $self->{cp1}->[1] -     $self->{p1}->[1]);
    $self->{F_Big}  = new Math::BigFloat '' . (                   3 * $self->{cp2}->[1] - 6 * $self->{cp1}->[1] + 3 * $self->{p1}->[1]);
    $self->{G_Big}  = new Math::BigFloat '' . (                                           3 * $self->{cp1}->[1] - 3 * $self->{p1}->[1]);
    $self->{H_Big}  = new Math::BigFloat '' . (                                                                       $self->{p1}->[1]);
    if (eval($self->{B_Big}->bstr) eq 0) {$self->{BdA_Big}=Math::BigFloat->bzero;} else {$self->{BdA_Big}=new Math::BigFloat '' . ($self->{B_Big}/$self->{A_Big});}
    if (eval($self->{C_Big}->bstr) eq 0) {$self->{CdA_Big}=Math::BigFloat->bzero;} else {$self->{CdA_Big}=new Math::BigFloat '' . ($self->{C_Big}/$self->{A_Big});}
    if (eval($self->{F_Big}->bstr) eq 0) {$self->{FdE_Big}=Math::BigFloat->bzero;} else {$self->{FdE_Big}=new Math::BigFloat '' . ($self->{F_Big}/$self->{E_Big});}
    if (eval($self->{G_Big}->bstr) eq 0) {$self->{GdE_Big}=Math::BigFloat->bzero;} else {$self->{GdE_Big}=new Math::BigFloat '' . ($self->{G_Big}/$self->{E_Big});}
    $self->{Am3_Big}=$self->{A_Big} * 3;
    $self->{Bm2_Big}=$self->{B_Big} * 2;
    $self->{Em3_Big}=$self->{E_Big} * 3;
    $self->{Fm2_Big}=$self->{F_Big} * 2;
    $self->{Am6_Big}=$self->{A_Big} * 6;
    $self->{Em6_Big}=$self->{E_Big} * 6;
    $self->{bigInitted}=1;
    }
sub initBigs {
    my $self=shift;
    my $bp1x =new Math::BigFloat '' .$self->{p1}->[0];
    my $bp1y =new Math::BigFloat '' .$self->{p1}->[1];
    my $bcp1x=new Math::BigFloat '' .$self->{cp1}->[0];
    my $bcp1y=new Math::BigFloat '' .$self->{cp1}->[1];
    my $bcp2x=new Math::BigFloat '' .$self->{cp2}->[0];
    my $bcp2y=new Math::BigFloat '' .$self->{cp2}->[1];
    my $bp2x =new Math::BigFloat '' .$self->{p2}->[0];
    my $bp2y =new Math::BigFloat '' .$self->{p2}->[1];
    $self->{A_Big}  =  ($bp2x - $bcp2x * 3 + $bcp1x * 3 -     $bp1x);
    $self->{B_Big}  =  (                   $bcp2x * 3 - $bcp1x * 6 + $bp1x * 3);
    $self->{C_Big}  =  (                                           $bcp1x * 3 - $bp1x * 3);
    $self->{D_Big}  =  (                                                                       $bp1x);
    $self->{E_Big}  =  ($bp2y - $bcp2y * 3 + $bcp1y * 3 -     $bp1y);
    $self->{F_Big}  =  (                   $bcp2y * 3 - $bcp1y * 6 + $bp1y * 3);
    $self->{G_Big}  =  (                                           $bcp1y * 3 - $bp1y * 3);
    $self->{H_Big}  =  (                                                                       $bp1y);
# SPEED ATTEMPTS
#    if (eval($self->{B_Big}->bstr) eq 0) {$self->{BdA_Big}=Math::BigFloat->bzero;} else {$self->{BdA_Big}=($self->{B_Big}/$self->{A_Big});}
#    if (eval($self->{C_Big}->bstr) eq 0) {$self->{CdA_Big}=Math::BigFloat->bzero;} else {$self->{CdA_Big}=($self->{C_Big}/$self->{A_Big});}
#    if (eval($self->{F_Big}->bstr) eq 0) {$self->{FdE_Big}=Math::BigFloat->bzero;} else {$self->{FdE_Big}=($self->{F_Big}/$self->{E_Big});}
#    if (eval($self->{G_Big}->bstr) eq 0) {$self->{GdE_Big}=Math::BigFloat->bzero;} else {$self->{GdE_Big}=($self->{G_Big}/$self->{E_Big});}
# this is iffy, because any non-number string plus zero is also zero, I think
    if (0 + $self->{B_Big}->bstr eq 0) {$self->{BdA_Big}=Math::BigFloat->bzero;} else {$self->{BdA_Big}=($self->{B_Big}/$self->{A_Big});}
    if (0 + $self->{C_Big}->bstr eq 0) {$self->{CdA_Big}=Math::BigFloat->bzero;} else {$self->{CdA_Big}=($self->{C_Big}/$self->{A_Big});}
    if (0 + $self->{F_Big}->bstr eq 0) {$self->{FdE_Big}=Math::BigFloat->bzero;} else {$self->{FdE_Big}=($self->{F_Big}/$self->{E_Big});}
    if (0 + $self->{G_Big}->bstr eq 0) {$self->{GdE_Big}=Math::BigFloat->bzero;} else {$self->{GdE_Big}=($self->{G_Big}/$self->{E_Big});}
    $self->{Am3_Big}=$self->{A_Big} * 3;
    $self->{Bm2_Big}=$self->{B_Big} * 2;
    $self->{Em3_Big}=$self->{E_Big} * 3;
    $self->{Fm2_Big}=$self->{F_Big} * 2;
    $self->{Am6_Big}=$self->{A_Big} * 6;
    $self->{Em6_Big}=$self->{E_Big} * 6;
    $self->{bigInitted}=1;
    }


sub getInfiniteSlopeThetas {
    my $self = shift;
    my @zts = $self->solveXPrimeforTheta(0);
    return wantarray ? @zts:$zts[0];
    }

sub getLength { #approx length
    my $self = shift;
    my $res = shift;
    my $start_theta = shift;
    my $end_theta = shift;
    if (!defined($res)) {$res=1000;}
    if (!defined($start_theta)) {$start_theta=0;}
    if (!defined($end_theta)) {$end_theta=1;}
    if ($end_theta<$start_theta) {my $tmp=$start_theta;$start_theta=$end_theta;$end_theta=$tmp;}
    if ($end_theta == $start_theta) {return 0;}
    my $length=0;
    my $inc = ($end_theta - $start_theta)/$res;
    #print "     $start_theta + $inc , l; $end_theta , inc;$inc\n";
    for (my $t=$start_theta + $inc;$t<=$end_theta;$t+=$inc) {
        $length+=sqrt(($self->bezierEvalXofT($t) - $self->bezierEvalXofT($t - $inc))**2 + ($self->bezierEvalYofT($t) - $self->bezierEvalYofT($t - $inc))**2);
        }
    return $length;
    }
sub precision {
    my $self=shift;
    if (defined($_[0])) {
        $self->{precision}=$_[0];
        #$self->{thetaprecision} = $self->{precision}/$self->{maxcoordval};
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
    #I've spent too much time here trying to handle those boundary cases.
    #I thought I had it handled, but now I find I've been using BigFloats in the min/max x/y comparisons
    # and that this code is taking maybe half of my execution time.
    #So, lets use normal numbers to speed stuff up, unless were really close to a boundary, and then we'll let the BigFloats in.
#print "looking for in range $coords->[0] , $coords->[1] \n";
    if (defined($coords->[0])) {
        if ( #really close to a boundary - allow any BigFloats to stay big
            #1 ||
            abs($coords->[0] - $self->{minx_notbig}) < 0.0001 ||
            abs($coords->[0] - $self->{maxx_notbig}) < 0.0001
            ) {
# SPEED ATTEMPTS - leave these evals here for now since hopefully this is rare?
# not so rare
            #warn "inrange using bigs 1\n";
            if ((eval($self->{minx})       < eval($coords->[0]) || eval($self->{minx})       eq eval($coords->[0])) && (eval($self->{maxx})       > eval($coords->[0]) || eval($self->{maxx})       eq eval($coords->[0]))) {
                $xok=1;
                }
            }
        else { #use definitely PERL size normal numbers for min/max x/y
            if ((    $self->{minx_notbig} < $coords->[0]        ||      $self->{minx_notbig} eq      $coords->[0])  && (     $self->{maxx_notbig} >      $coords->[0]  ||      $self->{maxx_notbig} eq      $coords->[0] )) {
                $xok=1;
                }
            }
        }
    if (defined($coords->[1])) {
        if (
            #1 ||
            abs($coords->[1] - $self->{miny_notbig}) < 0.0001 ||
            abs($coords->[1] - $self->{maxy_notbig}) < 0.0001
            ) {
            #print " uuuuuuubug\n";
            #print "    test1: eval($self->{miny})       < eval($coords->[1]) ?: ",((eval($self->{miny})       < eval($coords->[1]))?'yes':'no'),"\n";
            #print "  ||test2: eval(sprintf(\"%.14f\",$self->{miny}))        eq eval($coords->[1])) ?: ",((eval(sprintf("%.14f",$self->{miny}))        eq eval($coords->[1]))?'yes':'no'),"\n";
            #print "  &&\n";
            #print "    test3:    eval($self->{maxy})       > eval($coords->[1]) ?: ",((eval($self->{maxy})       > eval($coords->[1]))?'yes':'no'),"\n";
            #print "  ||test4:    eval($self->{maxy})       eq eval($coords->[1]) ?: ",((eval($self->{maxy})       eq eval($coords->[1]))?'yes':'no'),"\n";
            #sprintf in here is a recent hack that might get to stay if it doesn't cause other problems - "%.14f" worked for awhile, then had to drop it to "%.13f" to nail a later case of similar problem of a ~0 eq 0 situation where ~0 really should have been understood or calculated as zero. When will I rewrite to whole system with exact math and "robust predicates". After the Mac Aurthur Fellowship award, right? Is that spelled correctly?
# SPEED ATTEMPTS - leave these evals here for now since hopefully this is rare?
#warn "inrange using bigs 2\n";
            if ( (eval($self->{miny})       < eval($coords->[1]) || eval(sprintf("%.13f",$self->{miny}))        eq eval($coords->[1])) && (eval($self->{maxy})       > eval($coords->[1]) || eval($self->{maxy})       eq eval($coords->[1]))) {
                $yok=1;
                }
            }
        else { #use definitely PERL size normal numbers for min/max x/y
            if ( (     $self->{miny_notbig} <      $coords->[1]  ||                      $self->{miny_notbig}   eq      $coords->[1] ) && (     $self->{maxy_notbig} >      $coords->[1]  ||      $self->{maxy_notbig} eq      $coords->[1] )) {
                $yok=1;
                }
            }
        }
    return $xok,$yok;
    }
sub onSegment {
    my $self = shift;
    my $point= shift;
    return scalar(grep {$point->[1] - $_[0]->bezierEvalYofT($_) < $self->{precision}} $self->solveXforTheta($point->[0]));
    }


# dont like that we're making so many copies of two subroutines in getFeet
# try to make those once - but then need a way to get x and y params in them
# so going to put those on self - or could be package globals

sub find90 ($) { #dotproduct equals zero for two perpendicular vectors
#        return  ($Math::MikePath::BezierCubicSegment::activeseg->bezierEvalYofT($_[0]) - $Math::MikePath::BezierCubicSegment::activey)
#              *  $Math::MikePath::BezierCubicSegment::activeseg->slopeTangent_byTheta($_[0])
#              + ($Math::MikePath::BezierCubicSegment::activeseg->bezierEvalXofT($_[0]) - $Math::MikePath::BezierCubicSegment::activex); # * $tanvec->[1] ( which is == 1 )

        # unrolled bezier evals to do dot product math
        return  (  (((($Math::MikePath::BezierCubicSegment::activeseg->{E} * $_[0]) + $Math::MikePath::BezierCubicSegment::activeseg->{F}) * $_[0] + $Math::MikePath::BezierCubicSegment::activeseg->{G}) * $_[0] + $Math::MikePath::BezierCubicSegment::activeseg->{H})
                 - $Math::MikePath::BezierCubicSegment::activey )
              *  (  (($Math::MikePath::BezierCubicSegment::activeseg->{Em3} * $_[0]  +  $Math::MikePath::BezierCubicSegment::activeseg->{Fm2}) * $_[0] + $Math::MikePath::BezierCubicSegment::activeseg->{G})
                  / (($Math::MikePath::BezierCubicSegment::activeseg->{Am3} * $_[0]  +  $Math::MikePath::BezierCubicSegment::activeseg->{Bm2}) * $_[0] + $Math::MikePath::BezierCubicSegment::activeseg->{C})
                 )
              + ( (((($Math::MikePath::BezierCubicSegment::activeseg->{A} * $_[0]) + $Math::MikePath::BezierCubicSegment::activeseg->{B}) * $_[0] + $Math::MikePath::BezierCubicSegment::activeseg->{C}) * $_[0] + $Math::MikePath::BezierCubicSegment::activeseg->{D})
                 - $Math::MikePath::BezierCubicSegment::activex); # * $tanvec->[1] ( which is == 1 )


        }
sub find90_other ($) { #dotproduct equals zero for two perpendicular vectors
        # don't bother to unroll this one - it rarely gets used

        return  ($Math::MikePath::BezierCubicSegment::activeseg->bezierEvalXofT($_[0]) - $Math::MikePath::BezierCubicSegment::activex)
              *  $Math::MikePath::BezierCubicSegment::activeseg->slopeNormal_byTheta($_[0])
              + -1 * ( $Math::MikePath::BezierCubicSegment::activeseg->bezierEvalYofT($_[0]) - $Math::MikePath::BezierCubicSegment::activey ); # * $tanvec->[1] ( which is == 1 )

        }



sub getFeet {
    my $self=shift;
    my $x=shift;
    my $y=shift;
    #for each interval between critical points - critical due to features of x(i) and y(i). - hueristic and then root find to get any 90 degree intersections
    my @feet=();
    my %ds={};

# clean up the pre-speed-optimized mess after you've seen to real
# geometry output
# cut point lookup time roughly in half so far, mostly with changes here

#    my @dangers=sort {$a<=>$b} grep {!$ds{$_}++} map {ref($_)?0 + sprintf("%.20f",$_->bstr):$_} (@{$self->{extremeys_is}},@{$self->{extremexs_is}});
    # let's make non-BigFloat versions of extremeis during new(), so we don't
    # ever have to call bstr in here.
    # Oh nice - cuts point lookup time by 15%.
    #my @dangers=sort {$a<=>$b} grep {!$ds{$_}++} (@{$self->{extreme_is_sorted_notbig}});

    # We feel funny setting temp package variables
    # but this lets us define the two subroutines that use these
    # just once. Previously we created anonymous subroutines right here, on
    # every call to getFeet(). That didn't seem right either, making thousands
    # of subroutines.
    # (Originally this was about looking for speed ups, but this way seems
    # to take the same amount of exe time.)

    $Math::MikePath::BezierCubicSegment::activeseg=$self;
    $Math::MikePath::BezierCubicSegment::activex=$x;
    $Math::MikePath::BezierCubicSegment::activey=$y;

    for (my $i=1;$i<scalar(@{$self->{extreme_is_sorted_notbig}});$i++) {
        my $boundl=$self->{extreme_is_sorted_notbig}->[$i - 1];
        my $boundr=$self->{extreme_is_sorted_notbig}->[$i];

        #my $st1=$self->slopeTangent_byTheta($boundl);
        #my $st2=$self->slopeTangent_byTheta($boundr);
        # maybe shaves a little time off
        my $st1=$self->{tangents_at_sorted_extremes}->[$i - 1];
        my $st2=$self->{tangents_at_sorted_extremes}->[$i];

        my $functouse = (abs($st1) > 1000 || abs($st2) > 1000
                         #|| $st1 =~ /inf/ || $st2 =~ /inf/
                         # is this faster than pattern matches?
                         || $st1 eq 'inf' || $st1 eq '-inf'
                         || $st2 eq 'inf' || $st2 eq '-inf'
                        )
                      ? \&find90_other
                      : \&find90;

        my ($foot_i,$msg);
        # Precision here is pretty high. Slight speedup if you turn it down.
        # But how much is enough? (Too much is enough.)
        ($foot_i,$msg)=BrentsMethod($functouse,[$boundl,$boundr],$self->{precision}/1000,undef,'trying to find feet on Bezier in MikePath');
        if (!defined($msg)) {
            if (ref($foot_i)) { #then downgrade
# SPEED ATTEMPTS
#                $foot_i=eval(substr($foot_i->bstr,0,25));
                $foot_i=0 + sprintf("%.20f",$foot_i->bstr);
                }
            push(@feet,[$self->bezierEvalXofT($foot_i),$self->bezierEvalYofT($foot_i),$foot_i]);
            }
        else {
            # no foot found - possible, okay.
            }
        }
    return @feet;
    }
sub f_old { # f(x) = y
    my $self=shift;
    my $x=shift;
    my $origx=$x;#in case we downgrade x from a BigFloat, but then decide we still want to try using the BigFloat version
    my $needstobebig=0;
    for (my $i=0;$i< @{$self->{'xofidangerranges'}};$i++) {
        next if ($x <  $self->{'xofidangerranges'}->[$i]->[0]);
        next if ($x >  $self->{'xofidangerranges'}->[$i]->[1]);
        if ( (   $x >  $self->{'xofidangerranges'}->[$i]->[0]
              && $x <  $self->{'xofidangerranges'}->[$i]->[1]   )
           ||    $x eq $self->{'xofidangerranges'}->[$i]->[0]
           ||    $x eq $self->{'xofidangerranges'}->[$i]->[1]
           ) {
            $needstobebig=1;
            last;
            }
        }
    if ($needstobebig && !ref($x) && $Math::MikePath::enableCarefulfofx) {
        $x=Math::BigFloat->new(''.$x) if !ref($x);
        }
    elsif (!$needstobebig && ref($x)) {
        #print "NOT danger zone for f(x)! so downgrading x:$x\n";
        #$x=eval("$x");
        }
    my %dupsieve;
    my @ret;
    if ((ref($x) && !$self->{isLite})) {
        @ret = map {$self->bezierEvalYofTBig($_)} grep {!$dupsieve{ref($_)?$_->bstr:$_}++} $self->solveXforThetaBig($x);
        }
    else {
        @ret = map {$self->bezierEvalYofT($_)}    grep {!$dupsieve{$_}++}       $self->solveXforTheta($x);
        }
    if (!scalar(@ret)) {
        #print "nothing for f() - is x w/in range?: (($x>$self->{minx} || $x eq $self->{minx}) && ($x<$self->{maxx} || $x eq $self->{maxx})) : ",(($x>$self->{minx} || $x eq $self->{minx}) && ($x<$self->{maxx} || $x eq $self->{maxx})?'yes':'no'),"\n";
        if (($origx>$self->{minx} || $origx eq $self->{minx}) && ($origx<$self->{maxx} || $origx eq $self->{maxx})) {
            #@ret should have had something - try again somehow
            #print "###### ret should have had something for $origx : try again somehow - \n";
            #print "  (paths minx: $self->{minx})\n";
            #print "  (paths maxx: $self->{maxx})\n";
            #print "  (paths miny: $self->{miny})\n";
            #print "  (paths maxy: $self->{maxy})\n";
            #print "bez spec: \n[ $self->{p1}->[0], $self->{p1}->[1] ],\n[ $self->{cp1}->[0], $self->{cp1}->[1] ],\n[ $self->{cp2}->[0], $self->{cp2}->[1] ],\n[ $self->{p2}->[0],  $self->{p2}->[1] ] \n";
            #print "      or: C$self->{p1}->[0],$self->{p1}->[1],$self->{cp1}->[0],$self->{cp1}->[1],$self->{cp2}->[0], $self->{cp2}->[1],$self->{p2}->[0],$self->{p2}->[1]\n";
            #maybe you hit exactly on an endpoint?
            foreach my $cp ($self->{p1},$self->{p2}) {if ($cp->[0] eq $x) {push(@ret,$cp->[1]);print "used exact control point y\n";}}
            #maybe you hit exactly on minx or maxx?
            if ($x eq $self->{maxx} || $x eq $self->{minx}) { # already know we're not at end points, so should be at extreme or peak of curve
    #warn "temp dis new hit extreme in f(x)";
            #warn "snap to extreme not at endpoint in f(x)";
                my @xPrimeZeros=$self->solveXPrimeforTheta(0);
                if (scalar(@xPrimeZeros)>=1) {
                    #if multiple thetas, sort by distance from $x. Theta we want should give shortest dist, ideally zero
                    @xPrimeZeros = sort {abs($x - $self->bezierEvalXofT( $a )) <=> abs($x - $self->bezierEvalXofT( $b ))} @xPrimeZeros;
                    push(@ret,$self->bezierEvalYofT( $xPrimeZeros[0] ));
                    }
                }
            #then, if that didn't work, try some more
            if (!scalar(@ret)) {
                print "###### ret should have had something for $origx : try again somehow - \n";
                if (!$self->{bigInitted}) {
                    print "was lite, so initted bigs and tried big evals\n";
                    $self->initBigs();
                    $self->initDangerRanges();
                    $self->{isLite}=0;
                    }
                #now try again, using BigFloats all around, with original x (which might have been a BigFloat that we downgraded, since we're not always right about downgrading)
                my @ret_is = grep {!$dupsieve{ref($_)?$_->bstr:$_}++} $self->solveXforThetaBig(ref($origx)?$origx:Math::BigFloat->new(''.$origx));
                for (my $ti=0;$ti<@ret_is;$ti++) {
                    my $needstobebig=0;
                    for (my $i=0;$i<scalar(@{$self->{mxidangerranges}});$i++) {
                        print "x(i) test:  ($ret_is[$ti] > $self->{'mxidangerranges'}->[$i]->[0] || $ret_is[$ti] eq $self->{'mxidangerranges'}->[$i]->[0]) && ($ret_is[$ti] < $self->{'mxidangerranges'}->[$i]->[1] || $ret_is[$ti] eq $self->{'mxidangerranges'}->[$i]->[1])\n";
                        if (($ret_is[$ti] > $self->{'mxidangerranges'}->[$i]->[0] || $ret_is[$ti] eq $self->{'mxidangerranges'}->[$i]->[0]) && ($ret_is[$ti] < $self->{'mxidangerranges'}->[$i]->[1] || $ret_is[$ti] eq $self->{'mxidangerranges'}->[$i]->[1])) {
                            #print "DANGER ZONE FOR x(i)! t:$ret_is[$ti]\n";
                            $needstobebig=1;
                            last;
                            }
                        }
                    for (my $i=0;$i<scalar(@{$self->{myidangerranges}});$i++) {
                        print "y(i) test:  ($ret_is[$ti] > $self->{'myidangerranges'}->[$i]->[0] || $ret_is[$ti] eq $self->{'myidangerranges'}->[$i]->[0]) && ($ret_is[$ti] < $self->{'myidangerranges'}->[$i]->[1] || $ret_is[$ti] eq $self->{'myidangerranges'}->[$i]->[1])\n";
                        if (($ret_is[$ti] > $self->{'myidangerranges'}->[$i]->[0] || $ret_is[$ti] eq $self->{'myidangerranges'}->[$i]->[0]) && ($ret_is[$ti] < $self->{'myidangerranges'}->[$i]->[1] || $ret_is[$ti] eq $self->{'myidangerranges'}->[$i]->[1])) {
                            #print "DANGER ZONE FOR y(i)! t:$ret_is[$ti]\n";
                            $needstobebig=1;
                            last;
                            }
                        }
                    if (!$needstobebig) {
                        #print 'DOWNGRADED A THETA',"\n";
                        #$ret_is[$ti]=eval(substr($ret_is[$ti]->bstr,0,25));
                        }

                    }
                #@ret= map {ref($_)?$self->bezierEvalYofTBig($_):$self->bezierEvalYofT($_)} @ret_is;
                @ret= map {$self->bezierEvalYofTBig($_)} @ret_is;
                }

            if (!scalar(@ret)) {
                my $brute = sub {
                    return $origx - $self->bezierEvalXofTBig($_[0]);
                    };
                for (my $i=1;$i<@{$self->{extremexs_is}};$i++) {
                    my $bl=$self->{extremexs_is}->[$i - 1];
                    my $br=$self->{extremexs_is}->[$i];
                    my $er;
                    my $onetheta;
                    ($onetheta,$er) = Bisection($brute,[$bl,$br],0.0000001,($bl+$br)/2,'Brute Bisection to find at least one theta');
                    if (!defined($er)) {
                        push(@ret,$self->bezierEvalYofTBig($onetheta));
                        print "\nBrute found $onetheta\n\n";
                        }
                    else {
                        #print "Brute couldn't find anything either\n";
                        }
                    }
                }

            if (!scalar(@ret)) {
                warn "you should have found something for f(x) right?";
                }
            }
        }
    return wantarray?@ret:$ret[0];
    #return @ret;
    }
sub F { # F(y) = x
    my $self=shift;
    my $y=shift;
    my $origy=ref($y)?$y->copy():$y;
    my $needstobebig=0;
    for (my $i=0;$i<scalar(@{$self->{'yofidangerranges'}});$i++) {
        if (($y > $self->{'yofidangerranges'}->[$i]->[0] || $y eq $self->{'yofidangerranges'}->[$i]->[0]) && ($y < $self->{'yofidangerranges'}->[$i]->[1] || $y eq $self->{'yofidangerranges'}->[$i]->[1])) {
            #print "DANGER ZONE FOR F(y)! y:$y\n";
            $needstobebig=1;
            last;
            }
        }
    if ($needstobebig && !ref($y) && $Math::MikePath::enableCarefulFofy) {
        $y=Math::BigFloat->new(''.$y) if !ref($y);
        }
    else {
        #downgrade?
        }
    my %dupsieve;
    my @ret;
    if (ref($y) && !$self->{isLite}) {
        @ret = map {$self->bezierEvalXofTBig($_)} grep {!$dupsieve{ref($_)?$_->bstr:$_}++} $self->solveYforThetaBig($y);
        }
    else {
        @ret=map {$self->bezierEvalXofT($_)} grep {!$dupsieve{$_}++} $self->solveYforTheta($y);
        }
    if (!scalar(@ret)) {
        #print "nothing for F() - is y w/in range?: (($y>$self->{miny} || $y eq $self->{miny}) && ($y<$self->{maxy} || $y eq $self->{maxy})) : ",(($y>$self->{miny} || $y eq $self->{miny}) && ($y<$self->{maxy} || $y eq $self->{maxy})?'yes':'no'),"\n";
        if (($origy>$self->{miny} || eval("$origy") eq eval("$self->{miny}")) && ($origy<$self->{maxy} || eval("$origy") eq eval("$self->{maxy}"))) {
            #@ret should have had something - try again somehow
            #print "###### ret should have had something for $origy : try again somehow - \n";
            #print "  (paths minx: $self->{minx})\n";
            #print "  (paths maxx: $self->{maxx})\n";
            #print "  (paths miny: $self->{miny})\n";
            #print "  (paths maxy: $self->{maxy})\n";

            #maybe you hit exactly on an endpoint?
            foreach my $cp (($self->{p1},$self->{p2})) {if ($cp->[1] eq $y) {push(@ret,$cp->[0]);print "used exact control point x\n";}}
            #maybe you hit exactly on miny or maxy?
            if ($y eq $self->{maxy} || $y eq $self->{miny}) {
                my @yPrimeZeros=$self->solveYPrimeforTheta(0);
                if (scalar(@yPrimeZeros)>=1) {
                    #see similar section in f(x) function above for comments
                    #this is just adapted from that, swapping y for x
                    @yPrimeZeros = sort {abs($y - $self->bezierEvalYofT( $a )) <=> abs($y - $self->bezierEvalYofT( $b ))} @yPrimeZeros;
                    push(@ret,$self->bezierEvalXofT( $yPrimeZeros[0] ));
                    #push(@ret,$self->bezierEvalXofT( $yPrimeZeros[0] ) . ' '); # hack! the appended space is to let the duplicate numbers pass through a string-compare-based duplicate filter downstream that (hopefully) sees diffrent string values for the identical number values
                    # no, don't want that hack. think I'm handling single result okay, and if not, should be.
                    #print "experimental F(y) shortcut, where y==maxy or miny, on probation. So is it working okay?\n";
                    #print "  by the way - returned two identical points for this, where probably want to return just one - fix board.pm stuff so it will work with one point in this case, and then don't return second duplicate point here\n";
                    }
                }

            #then, if that didn't work, try some more
            print "###### ret should have had something for $origy : try again somehow - \n";
            if (!scalar(@ret)) {
                if (!$self->{bigInitted}) {
                    print "was lite, so initted bigs and tried big evals\n";
                    $self->initBigs();
                    $self->initDangerRanges();
                    $self->{isLite}=0;
                    }
                #now try again, using BigFloats all around, with original y (which might have been a BigFloat that we downgraded, since we're not always right about downgrading)
                my @ret_is = grep {!$dupsieve{ref($_)?$_->bstr:$_}++} $self->solveYforThetaBig(ref($origy)?$origy:Math::BigFloat->new(''.$origy));
                for (my $ti=0;$ti<@ret_is;$ti++) {
                    my $needstobebig=0;
                    for (my $i=0;$i<scalar(@{$self->{mxidangerranges}});$i++) {
                        #print "x(i) test:  ($ret_is[$ti] > $self->{'mxidangerranges'}->[$i]->[0] || $ret_is[$ti] eq $self->{'mxidangerranges'}->[$i]->[0]) && ($ret_is[$ti] < $self->{'mxidangerranges'}->[$i]->[1] || $ret_is[$ti] eq $self->{'mxidangerranges'}->[$i]->[1])\n";
                        if (($ret_is[$ti] > $self->{'mxidangerranges'}->[$i]->[0] || $ret_is[$ti] eq $self->{'mxidangerranges'}->[$i]->[0]) && ($ret_is[$ti] < $self->{'mxidangerranges'}->[$i]->[1] || $ret_is[$ti] eq $self->{'mxidangerranges'}->[$i]->[1])) {
                            #print "DANGER ZONE FOR x(i)! t:$ret_is[$ti]\n";
                            $needstobebig=1;
                            last;
                            }
                        }
                    for (my $i=0;$i<scalar(@{$self->{myidangerranges}});$i++) {
                        #print "y(i) test:  ($ret_is[$ti] > $self->{'myidangerranges'}->[$i]->[0] || $ret_is[$ti] eq $self->{'myidangerranges'}->[$i]->[0]) && ($ret_is[$ti] < $self->{'myidangerranges'}->[$i]->[1] || $ret_is[$ti] eq $self->{'myidangerranges'}->[$i]->[1])\n";
                        if (($ret_is[$ti] > $self->{'myidangerranges'}->[$i]->[0] || $ret_is[$ti] eq $self->{'myidangerranges'}->[$i]->[0]) && ($ret_is[$ti] < $self->{'myidangerranges'}->[$i]->[1] || $ret_is[$ti] eq $self->{'myidangerranges'}->[$i]->[1])) {
                            #print "DANGER ZONE FOR y(i)! t:$ret_is[$ti]\n";
                            $needstobebig=1;
                            last;
                            }
                        }
                    if (!$needstobebig) {
                        #print 'DOWNGRADED A THETA',"\n";
                        #$ret_is[$ti]=eval(substr($ret_is[$ti]->bstr,0,25));
                        }
                    }
                #@ret= map {ref($_)?$self->bezierEvalYofTBig($_):$self->bezierEvalYofT($_)} @ret_is;
                @ret= map {$self->bezierEvalXofTBig($_)} @ret_is;
                }



            if (!scalar(@ret)) {
                my $brute = sub {
                    return $origy - $self->bezierEvalYofTBig($_[0]);
                    };
                for (my $i=1;$i<@{$self->{extremeys_is}};$i++) {
                    my $bl=$self->{extremeys_is}->[$i - 1];
                    my $br=$self->{extremeys_is}->[$i];
                    my $er;
                    my $onetheta;
                    ($onetheta,$er) = Bisection($brute,[$bl,$br],0.0000001,($bl+$br)/2,'Brute Bisection to find at least one theta');
                    if (!defined($er)) {
                        push(@ret,$self->bezierEvalXofTBig($onetheta));
                        print "\nBrute found $onetheta\n\n";
                        }
                    else {
                        #print "Brute couldn't find anything either\n";
                        }
                    }
                }

            if (!scalar(@ret)) {
                warn "you should have found something for f(x) right?";
                }

            }
        }
    return wantarray?@ret:$ret[0];
    #return @ret;
    }
sub secondDerivative { # THIS SHOULD BE WHERE I FIRST WORKED OUT THE RIGHT APPROACH TO 2ND DERIVATIVE IN XY FRAME (dy/dx)
                       # BUT BELIEVE I'M REDOING ESSENTIALLS OF THIS ABOVE IN f_prime() FUNTCTION ABOVE, SO THIS GOES AWAY WHEN YR SURE THAT'S ALL WORKING.
    my $self = shift;
    my $x    = shift;
    my $y    = @_?shift:undef;
    my @thetas;
#    if    (defined($x)) {@thetas = $self->solveXforTheta($x);}
    if    (defined($x)) {@thetas = $self->solveXforTheta_newself($x,2);} # gives theta AND the first and second derivs of theta with respect to x, by hacking a cubic solver
    elsif (defined($y)) {die "not yet, but s/b cut and paste easy";}
    else {return;}
    my @ret;
    foreach my $theta_and_deriv (@thetas) {

        my $theta = $theta_and_deriv->[0];
        my $thetaprime = $theta_and_deriv->[1];
        my $thetaprime2 = $theta_and_deriv->[2];

        my $yp=$self->bezierEvalYPrimeofT($theta);
        my $xp=$self->bezierEvalXPrimeofT($theta);
        my $ypp=$self->bezierEvalYDoublePrimeofT($theta);
        my $xpp=$self->bezierEvalXDoublePrimeofT($theta);

        # This matches yp/xp !
        # thatis yprime of theta * dtheta/dx = yprime of theta/xprime of theta
        my $ypofx = $self->{Em3} * $theta**2 * $thetaprime  +  $self->{Fm2} * $theta * $thetaprime + $self->{G} * $thetaprime;

        # for second derivative, have to use chain and product rules to differentiate the above.

        my $yppofx = $self->{Em3} * ($theta**2 * $thetaprime2 + 2*$theta*$thetaprime**2)
                   + $self->{Fm2} * ($theta    * $thetaprime2 +          $thetaprime**2)
                   + $self->{G}   * $thetaprime2;

        return $yppofx if (!$main::dodebug);

        # JUST PASSING BY LATER... and want to note that I think the above was/is
        # working fine when I was last deving here. And I just didn't bother to
        # clean up what's below, partly because as it says right there, there's maybe useful debug stuff
        # and partly because I just never have the mental spacetime in this environment to burnish the gems...
        # But, yeah, looks like I figured graph-space 2nd deriv for cubic bezier,
        # and I think that means I can calc graph-space curvature (K) or radius of curvature too.
        # (I used this as part of doing that for rail joiner curvature - but that had extra complicating steps.
        # In other words, calcing curvature for simple cubic bez case is probably simple, since I did curvature for something more complex...)
        # An SVG Path library that offers analytic graph-space 2nd derivitive and curvature calcs is
        # something some people will want.
        # This needs test suite and engaging documentation with app examples.
        # Okay, also, just noticing that looks like all the old calcs above for $yp,$xp,$ypp, and $xpp can go.
        # Originally thought I would be working with them, or at least comparing to them, but
        # you can see the math actually used for the return val didn't use them. (Double check! 'cause I'm just quick passing through.)


        # You'll want some of this debug mess (and more probably)
        # as you hit each of the 5 expressions in the solver - especially
        # the two Det >= 0 cases.

        # no, because thetaprime is really an expression containing x, so
        # you also want thetadoubleprimeofx now
        #

        warn "\n";
        warn "theta        : $theta\n";
        warn "1st derrrrr??: $ypofx            vs yp/xp: ",($yp/$xp),"\n";
        warn "2nd derrrrr??: $yppofx\n";
        warn "\n";


        my $hh=0.001; # approx can actually get worse with higher res hh, so careful
        my ($th1,$dr1,$dr1_2) = map @$_, $self->solveXforTheta_newself($x-$hh,2);
        my ($th2,$dr2,$dr2_2) = map @$_, $self->solveXforTheta_newself($x    ,2);
        my ($th3,$dr3,$dr3_2) = map @$_, $self->solveXforTheta_newself($x+$hh,2);
        my $tm1 = ($th2-$th1)/$hh;
        my $tm2 = ($th3-$th2)/$hh;
        my $tdx_approx = ($tm2-$tm1)/$hh;
        warn "approx test thetas : $th1 , $th2 , $th3\n";
        warn "approx test theta's : $dr1 , $dr2 , $dr3\n";
        warn "approx test theta''s : $dr1_2 , $dr2_2 , $dr3_2\n";
        warn "\n";
        warn "tPrimeOfX: 2 approx slopes vs calc tan in middle:\n";
        warn "$tm1\n$dr2 <-- calc s/b in middle\n$tm2\n";
        warn "\n";
        warn "tDoublePrimeOfX: approx vs 3 s/b bracketing calcs\n";
        warn "$tdx_approx = ($tm2 - $tm1)/$hh vs\n$dr1_2\n$dr2_2\n$dr3_2\n";

        my ($th4,$dr4) = map @$_, $self->solveXforTheta_newself($x+1,'doderiv');
        warn "ifonly: ($th2 + $dr2) = ",($th2 + $dr2)," vs $th4 hey actually close! maybe right\n";


# not quite here...
# must be of x not of i

        my $lowdeehiminushideelow = ($xp*$ypp - $yp*$xpp);
        my $lowlow = $xpp**2;

        #if ($xpp eq '0') {push(@ret,(($lowdeehiminushideelow < 0)?'-':'+').'inf');}
        #if ($lowdeehiminushideelow eq '0') {push(@ret,(($lowlow < 0)?'-':'+').'0');}

        # (yp/xp)' = low dee high minus high dee low over low low
        # (xp*ypp - yp*xpp) / xpp^2



warn "\n\nyp/xp eq? yp * dt/dx ?\n";
warn "",($yp/$xp), " eq? ",($yp*$thetaprime)," or ",$ypofx,"\n";

warn "  ypp stuff $yppofx : ",($ypp/$xpp), " eq? ",($yp*$thetaprime*$thetaprime2),"\n";

die "stop w/in 2nd deriv";
return $yppofx;


        #else {
            push(@ret, ($thetaprime*2) * ($lowdeehiminushideelow / $lowlow));
            # might it be just ?? push(@ret, $ypp / $xpp);
# looking into hodograph approach - which I might have already mimicked? in my first derivative approach? but using the cubic formula. hmmm doesn't seem like it should work
            #push @ret, $ypp / $xpp;
        #    }

        }
    #warn "ret:",join(',',@ret),"\n";
    return @ret;
    }
sub slopeTangent {
    my $self = shift;
    my $x    = shift;
    my $y    = @_?shift:undef;
    my @thetas;
    if    (defined($x)) {@thetas = $self->solveXforTheta($x);}
    elsif (defined($y)) {@thetas = $self->solveYforTheta($y);}
    else {return;}
    my @ret;
    foreach my $theta (@thetas) {
        my $yp=$self->bezierEvalYPrimeofT($theta);
        my $xp=$self->bezierEvalXPrimeofT($theta);
        if ($xp eq '0') {push(@ret,(($yp < 0)?'-':'+').'inf');}
        if ($yp eq '0') {push(@ret,(($xp < 0)?'-':'+').'0');}
        else {
            push(@ret,$yp/$xp);
            }
        }
    return @ret;
    }
sub slopeTangent_byTheta {
    my $self = shift;
    my $theta = shift;
    my $yp=$self->bezierEvalYPrimeofT($theta);
    my $xp=$self->bezierEvalXPrimeofT($theta);
    # adding zero here converts possible '-0' to 0 for this comparison
    if (0 + $xp eq '0') {return (($yp < 0)?'-':'+').'inf';}
    if (0 + $yp eq '0') {return (($xp < 0)?'-':'+').'0';}
    else {
        return $yp/$xp;
        }
    }
sub angleTangent {
    my $self = shift;
    my $x    = shift;
    my $y    = @_?shift:undef;
    my @thetas;
    if    (defined($x)) {@thetas = $self->solveXforTheta($x);}
    elsif (defined($y)) {@thetas = $self->solveYforTheta($y);}
    else {return;}
    my @ret;
    foreach my $theta (@thetas) {
        my $yp=$self->bezierEvalYPrimeofT($theta);
        my $xp=$self->bezierEvalXPrimeofT($theta);
        push(@ret,atan2($yp,$xp));
        }
    return @ret;
    }
sub angleTangent_byTheta {
    my $self = shift;
    my $theta = shift;
    my $yp=$self->bezierEvalYPrimeofT($theta);
    my $xp=$self->bezierEvalXPrimeofT($theta);
    return atan2($yp,$xp);
    }
sub fDoublePrime { # THINK THIS IS ALL WRONG AND WHAT MY NEW STUFF ACTUALLY DOES RIGHT. Simple "lo deehi minus hi deelo over lo lo" was probably wrong wrong, and can't see that I actually used this function anywhere.
    my $self = shift;
    my $x    = shift;
    my $y    = @_?shift:undef;
    my @thetas;
    if    (defined($x)) {@thetas = $self->solveXforTheta($x);}
    elsif (defined($y)) {@thetas = $self->solveYforTheta($y);}
    else {return;}
    my @ret;
    foreach my $theta (@thetas) {
        my $yp=$self->bezierEvalYPrimeofT($theta);
        my $ypp=$self->bezierEvalYDoublePrimeofT($theta);
        my $xp=$self->bezierEvalXPrimeofT($theta);
        my $xpp=$self->bezierEvalXDoublePrimeofT($theta);
        if ($xp eq '0') {push(@ret,(($yp < 0)?'-':'+').'inf');}
        if ($yp eq '0') {push(@ret,(($xp < 0)?'-':'+').'0');}
        else {push(@ret,($xp*$ypp)-($yp*$xpp)/($xp**2));} #lo deehi minus hi deelo over lo lo, right?
        }
    return @ret;
    }
sub slopeNormal {
    my $self = shift;
    my $x    = shift;
    my $y    =@_?shift:undef;
    my @ret;
    for my $slopeTangent ($self->slopeTangent($x,$y)) {
        my $negRecip;
        if ($slopeTangent =~ /([\-\+]?)inf$/i) {
            my $sign = '';
            if (length($1)) {if ($1 eq '-') {$sign='+';} else {$sign='-';}}
            $negRecip=$sign.'0';
            }
        elsif ($slopeTangent =~ /^([\-\+]?)0$/) {
            my $sign = '';
            if (length($1)) {if ($1 eq '-') {$sign='+';} else {$sign='-';}}
            $negRecip=$sign.'inf';
            }
        else {
            $negRecip=-1/$slopeTangent;
            }
        push(@ret,$negRecip);
        }
    return @ret;
    }

sub slopeNormal_byTheta {
    my $self = shift;
    my $theta = shift;
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
sub angleNormal {
    my $self = shift;
    my $x    = shift;
    my $y    = @_?shift:undef;
    my @thetas;
    if    (defined($x)) {@thetas = $self->solveXforTheta($x);}
    elsif (defined($y)) {@thetas = $self->solveYforTheta($y);}
    else {return;}
    my @ret;
    foreach my $theta (@thetas) {
        my $yp=$self->bezierEvalYPrimeofT($theta);
        my $xp=$self->bezierEvalXPrimeofT($theta);
        push(@ret,atan2(-$xp,$yp));
        }
    return @ret;
    }
sub angleNormal_byTheta {
    my $self = shift;
    my $theta = shift;
    my $yp=$self->bezierEvalYPrimeofT($theta);
    my $xp=$self->bezierEvalXPrimeofT($theta);
    return atan2(-$xp,$yp);
    }
sub point {
    my $self=shift;
    my $theta=shift;
    for (my $i=0;$i<scalar(@{$self->{'mxidangerranges'}});$i++) {
        if (($theta > $self->{mxidangerranges}->[$i]->[0] || $theta eq $self->{mxidangerranges}->[$i]->[0]) && ($theta < $self->{mxidangerranges}->[$i]->[1] || $theta eq $self->{mxidangerranges}->[$i]->[1]) && !ref($theta)) {
            #print "DANGER ZONE FOR point(theta) for x. $theta:ref?:",ref($theta),"\n";
            #print "DZ  theta: $theta\n";
            $theta=Math::BigFloat->new($theta) if !ref($theta);
            #print "DZB theta: $theta\n";
            }
        }
    for (my $i=0;$i<scalar(@{$self->{'myidangerranges'}});$i++) {
        if (($theta > $self->{myidangerranges}->[$i]->[0] || $theta eq $self->{myidangerranges}->[$i]->[0]) && ($theta < $self->{myidangerranges}->[$i]->[1] || $theta eq $self->{myidangerranges}->[$i]->[1]) && !ref($theta)) {
            #print "DANGER ZONE FOR point(theta) for y. $theta:$theta\n";
            $theta=Math::BigFloat->new($theta) if !ref($theta);
            }
        }
    #if (ref($theta)) {print "  your XY sir: ",$self->bezierEvalXofT($theta),",",$self->bezierEvalYofT($theta),"\n";}
    my $ret;
    if (ref($theta)) {$ret=[$self->bezierEvalXofTBig($theta),$self->bezierEvalYofTBig($theta)];}
    else {$ret=[$self->bezierEvalXofT($theta),$self->bezierEvalYofT($theta)];}
    return $ret;
    }

sub bezierEvalXofT {
    my $self = shift;
    my $t    = shift;
    return ((($self->{A} * $t) + $self->{B}) * $t + $self->{C}) * $t + $self->{D};
    }
sub bezierEvalXofTBig {
    my $self = shift;
    my $t    = shift;
    return ((($self->{A_Big} * $t) + $self->{B_Big}) * $t + $self->{C_Big}) * $t + $self->{D_Big};
    }
sub bezierEvalYofT {
    my $self = shift;
    my $t    = shift;
    return ((($self->{E} * $t) + $self->{F}) * $t + $self->{G}) * $t + $self->{H};
    }
sub bezierEvalYofTBig {
    my $self = shift;
    my $t    = shift;
    return ((($self->{E_Big} * $t) + $self->{F_Big}) * $t + $self->{G_Big}) * $t + $self->{H_Big};
    }
sub bezierEvalXPrimeofT {
    my $self = shift;
    my $t    = shift;
    # x'= 3At^2+  2Bt  +     C
    #   =(3At  +  2B)t +     C
    my $ret = ($self->{Am3} * $t  +  $self->{Bm2}) * $t + $self->{C};
    if ($ret == 0) {$ret = (($self->bezierEvalXDoublePrimeofT($t) < 0)?'-':'').$ret;} #? is that enough? useful?
    return $ret;
    }
sub bezierEvalYPrimeofT {
    my $self = shift;
    my $t    = shift;
    # y'= 3Et^2+  2Ft  +     G
    #   =(3Et  +  2F)t +     G
    my $ret = ($self->{Em3} * $t  +  $self->{Fm2}) * $t + $self->{G};
    if ($ret eq '0') {$ret = (($self->bezierEvalYDoublePrimeofT($t) < 0)?'-':'').$ret;} #? is that enough? useful?
    return $ret;
    }
sub bezierEvalXDoublePrimeofT {
    my $self = shift;
    my $t    = shift;
    # x''= 6At +  2B
    return $self->{Am6} * $t + $self->{Bm2};
    }
sub bezierEvalYDoublePrimeofT {
    my $self = shift;
    my $t    = shift;
    # y''= 6Et +  2F
    return $self->{Em6} * $t + $self->{Fm2};
    }
sub bezierSolve { #obsolete?
    my $self = shift;
    my $x    = shift;
    my $y    = shift;
    return [$self->solveX($x),$self->solveY($y)];
    }
sub solveXforTheta {
    my $self = shift;
    my $x    = shift;
    # x = At^3 +  Bt^2 +     Ct +     D
    # 0 = At^3 +  Bt^2 +     Ct +     D -x
    # 0 =  t^3 + (B/A)t^2 + (C/A)t + (D-x)/A
    #print "for cubic formula: $self->{BdA} , $self->{CdA} , ",(($self->{D} - ((ref($x))?$x:$x))/$self->{A}),"\n";
    #print "pre filter cubic formula returns:",join(" , ",&cubicformula($self->{BdA},$self->{CdA},($self->{D} - ((ref($x))?$x:$x))/$self->{A},0)),"\n";
    return sort {$a<=>$b} grep {(1 > $_ || 1 eq $_) && ($_ > 0 || $_ eq 0)} &cubicformula($self->{BdA},$self->{CdA},($self->{D} - ((ref($x))?$x:$x))/$self->{A},1);
    }


# Adapted from our own CubicFormula.pm
# and optimized for use with cubic Beziers -
# just real roots, some values pre-computed,
# willing to make ugly if it's faster.
# Also computes first and second derivatives.
# Would also like to dev this in direction of
# knowing numerically unstable scenarios and fixing or avoiding them,
# without, perhaps, resorting to BigFloats as much, or at all.

sub solveXforTheta_newself {
    my $self = shift;
    my $x    = shift;
    my $doThetaPrime = @_?shift:0;

    my $C = ($self->{D} - $x)/$self->{A};
    my $bigpi = (ref($self->{cubA}) || ref($self->{cubB}) || ref($C))?1:0; # boolean, whether to use BigFloat version of pi
    my $thispi=$bigpi?$pi_big:$pi;
    my $thistwopi=$bigpi?$twopi_big:$twopi;
    my $thisfourpi=$bigpi?$fourpi_big:$fourpi;

    my @inrange;
    my @deriv;
    my @deriv2;

    my ($R, $D, $X1, $X2, $X3);

    # do the full R calc here, because if you precalc some of it's parts
    # you get rounding sometimes on the last digit of those intermediates
    # and the whole result doesn't come out the same as doing it all at once
    # all together like this.
    $R = (9.0*$self->{cubA}*$self->{cubB} - 27.0*$C - 2.0*$self->{cubA}**3)/54.0;

    $D = $self->{cubQ}**3 + $R**2;                  # polynomial discriminant

    # was just reading about rationalizing to avoid subtraction
    # since subtraction can be such a hungry hippo with significant digits
    # that Q can be negative, and then that's like R^2 - Q^3
    # could you rationalize to (R^2 - Q^3)(R^2 + Q^3) / (R^2 + Q^3)
    # (R^4 - Q^6) / (R^2 + Q^3)
    # sorta think that's no good, but maybe just maybe?
    # would need tests to see what it's effect is on the typical problem areas
    # in here. look at any other subtractions in here too. Maybe we can eliminate
    # resort to BigFloats?
    # Might be better to experiment on this with a fresh copy of standalone
    # cubic solver, and set up direct tests for each problem area in these calcs,
    # prob starting with that bit a bit below where you upgrade to bigfloat
    # when abs($R - $self->{cubsqrtofnegQcubed}) < 0.000000000001.
    # that might become (R^2 - (-Q^3)) / ($R + $self->{cubsqrtofnegQcubed})
    # hokum?

    my $allReal = 1;

    if ($D > 0 || $D eq 0) { # 1 real, and 2 complex or 2 real duplicate roots
        $allReal = 0;
        my $sqrtD=(ref($D)?bigsqrt($D):sqrt($D));
        my $preS=$R + $sqrtD;
        my $preT=$R - $sqrtD;
        my $S = (($preS < 0)?-1:1) * (abs($preS)**(1/3));
        my $T = (($preT < 0)?-1:1) * (abs($preT)**(1/3));
        my $SpT = $S + $T;
        $X1 = -$self->{cubAd3} + $SpT;              # real root
        push(@inrange, $X1) if (1 > $X1 || 1 eq $X1) && ($X1 > 0 || $X1 eq 0);

        my $deriv  = (1/(6  * $sqrtD*$self->{A}))
                   * (  ($preS/(($preS**2)**(1/3)))
                      - ($preT/(($preT**2)**(1/3)))
                     );

        my $deriv2 = (1/(12 * $D    *$self->{A}**2))
                   * (  ($preS/(($preS**2)**(1/3)))
                      * ( (1/3) - $R/($sqrtD) )
                      +
                        ($preT/(($preT**2)**(1/3)))
                      * ( (1/3) + $R/($sqrtD) )
                     )
        ;

        push @deriv , $deriv;
        push @deriv2, $deriv2;

        #if ((1 > $X1 || 1 eq $X1) && ($X1 > 0 || $X1 eq 0)) {
        #    $main::dodebug=1;
        #    }



        if ($D eq 0) { # $D==0 case, duplicate real roots
            $allReal=1;
            $X2 = -$self->{cubAd3} - $SpT/2; # real part of complex root (but complex part is 0 here, so real)
warn "yall need a deerivdiv herabouts too";
warn "this right??";
        push @deriv , $deriv /2, $deriv /2;
        push @deriv2, $deriv2/2, $deriv2/2;
die "die here until ready to figure out";

            # sort as you go
            if ((1 > $X2 || 1 eq $X2) && ($X2 > 0 || $X2 eq 0)) {
                if ($X2 >= $X1) {push( @inrange, $X2, $X2);}
                else {unshift(@inrange, $X2,$X2);}
                }
            }
        }
    else {                                          # 3 distinct real roots
        # if only we could do it in one line like this
        #my $th = acos($R/sqrt(-1 * ($Q**3)));
        # but numerical instability lurks in here.

        my $toAcos;
        my $sqrtofnegQcubed = (
            ref($self->{cubQ})
            # this check is key to fixing rare numerical instability bug
            || abs($R - $self->{cubsqrtofnegQcubed}) < 0.000000000001
            )
            ? $self->cubsqrtofnegQcubed_big()
            : $self->{cubsqrtofnegQcubed};

        # upgrade R if nec. so / operator overload works with R on left, setting Big mode
        $R = Math::BigFloat->new(''.$R) if (ref($sqrtofnegQcubed) && ! ref($R));

        $toAcos = $R / $sqrtofnegQcubed;

        my $th;

        if (ref($toAcos)) {
            # this stuff with "$noom" is from debug of rare numerical instabliltiy.
            # if tempted to simplify any of this, do that in context of making all the math here robust.
            my $noom=$toAcos->copy();
            #if ($noom eq 'NaN') {die "noom is NaN after big copy from toAcos\n";}
            $noom->bpow(2)->bmul(-1.0)->badd(1.0);
            #if ($noom eq 'NaN') {die "noom is NaN after series of pow, mul, add operations\n";}

            if ($noom < 0) {
                # snap neg to zero - not sure that's good, but maybe -
                # guessing this should not ever go negative
                $noom = Math::BigFloat->bzero;
                }
            my $newnoom;
            if ($noom <= 0) {$newnoom=Math::BigFloat->bzero;}
            else {$newnoom=bigsqrt($noom);}
            #if ($newnoom eq 'NaN') {die "noom is NaN after bigsqrt. before, noom was: ",$noom->bstr,"\n";}
            #else {
                $noom = $newnoom;
            #    }
            # This arctan is not the same as atan2 or even atan - it's input is
            # limited to abs(x)<=1. Math::Big::arctan uses a slower-to-converge
            # Taylor series that is prone to not converging near 1, so I made my
            # own arctan based on some faster Euler series.
            my $abscmp=$toAcos->bacmp($noom);
            if ($abscmp < 0 || $abscmp eq 0) {
                my $quo=$toAcos->copy()->bdiv($noom);
                #if ($quo eq 'NaN') {die "pre arctan_Euler call 1: have a NaN derived from toAcos/noom : ",$toAcos->bstr," / ",$noom->bstr,"\n"}
                $th=Math::CubicFormula::arctan_Euler($quo,25);
                $th->bmul(-1)->badd($piovertwo_big);
                }
            else {
                my $quo=$noom->copy()->bdiv($toAcos);
                #if ($quo eq 'NaN') {die "pre arctan_Euler call 2: have a NaN derived from noom/toAcos : ",$noom->bstr," / ",$toAcos->bstr,"\n"}
                $th=Math::CubicFormula::arctan_Euler($quo,25);
                }
            $th = Math::CubicFormula::atanfix($th,$toAcos,$noom);
            }
        else {
            #   arccosine
            $th=atan2(sqrt(1 - $toAcos * $toAcos),$toAcos);

            # arctan_Euler() might be called for here sometimes too
            # could that make it robust while using all doubles maybe?



            }

        # think you can work out the $th range that works
        # and then grep the three ($th+C)/3 angles so you don't always have to
        # calculate all three roots below.Often a 60% speed up for this section, right?
        # Prob not. but should register.

        if ($bigpi || ref($th) || ref($sqrtofnegQcubed) || ref($R)) { # was testing $A too, but
                                                        # bigpi should be true if ref($A) is true
                                                        # - probably true for (some of) the others too, but check
            if (ref($sqrtofnegQcubed)) {
                my $cubtwotimessqrtofnegQ_big = $self->cubtwotimessqrtofnegQ_big();
                $X1 = $cubtwotimessqrtofnegQ_big * Math::Big::cos( $th               /3,40) - $self->{cubAd3};
                $X2 = $cubtwotimessqrtofnegQ_big * Math::Big::cos(($th + $thistwopi )/3,40) - $self->{cubAd3};
                $X3 = $cubtwotimessqrtofnegQ_big * Math::Big::cos(($th + $thisfourpi)/3,40) - $self->{cubAd3};
                push @inrange, grep {(1 > $_ || 1 eq $_) && ($_ > 0 || $_ eq 0)} ($X1, $X2, $X3);
                }
            else {
                $X1 = $self->{cubtwotimessqrtofnegQ}     * Math::Big::cos( $th               /3,40) - $self->{cubAd3};
                $X2 = $self->{cubtwotimessqrtofnegQ}     * Math::Big::cos(($th + $thistwopi )/3,40) - $self->{cubAd3};
                $X3 = $self->{cubtwotimessqrtofnegQ}     * Math::Big::cos(($th + $thisfourpi)/3,40) - $self->{cubAd3};
                push @inrange, grep {(1 > $_ || 1 eq $_) && ($_ > 0 || $_ eq 0)} ($X1, $X2, $X3);
                }
            }

            # OPTIMIZE this else, then copy approach to above two versions of same


        else {
            #my $twotimessqrtofnegQ=2*sqrt(-1 * $Q);

            # if it's more often not in range, this speeds things up
            # if it's usually in range this slows it down

            # The range you should really precalc is the input x range
            # like the C range, I guess.
            # Maybe you work back to that?
            # Find upstream ranges from the range you now have.
            # Go a step back in theta calc, and figure corresponding ranges for that.
            #toAcos might be first point for that

            #my $cos1=cos($th/3);
            #push(@inrange, $self->{cubtwotimessqrtofnegQ} * $cos1 - $self->{cubAd3}) if ($cos1 <= $self->{cubtheta_slope_upper_bound} && $cos1 >= $self->{cubtheta_slope_lower_bound});
            #my $cos2=cos(($th + $thistwopi)/3);
            #push(@inrange, $self->{cubtwotimessqrtofnegQ} * $cos2 - $self->{cubAd3}) if ($cos2 <= $self->{cubtheta_slope_upper_bound} && $cos2 >= $self->{cubtheta_slope_lower_bound});
            #my $cos3=cos(($th + $thisfourpi)/3);
            #push(@inrange, $self->{cubtwotimessqrtofnegQ} * $cos3 - $self->{cubAd3}) if ($cos3 <= $self->{cubtheta_slope_upper_bound} && $cos3 >= $self->{cubtheta_slope_lower_bound});
=cut
            my $cos1=cos($th/3);
            if ($cos1 <= $self->{cubtheta_slope_upper_bound} && $cos1 >= $self->{cubtheta_slope_lower_bound}) {
                my $root = $self->{cubtwotimessqrtofnegQ} * $cos1 - $self->{cubAd3};
                push(@inrange, $root) if (($root > 0 || $root eq 0) && ($root < 1 || $root eq 1));

                }
            my $cos2=cos(($th + $thistwopi)/3);
            if ($cos2 <= $self->{cubtheta_slope_upper_bound} && $cos2 >= $self->{cubtheta_slope_lower_bound}) {
                my $root = $self->{cubtwotimessqrtofnegQ} * $cos2 - $self->{cubAd3};
                push(@inrange, $root) if (($root > 0 || $root eq 0) && ($root < 1 || $root eq 1));
                }
            my $cos3=cos(($th + $thisfourpi)/3);
            if ($cos3 <= $self->{cubtheta_slope_upper_bound} && $cos3 >= $self->{cubtheta_slope_lower_bound}) {
                my $root = $self->{cubtwotimessqrtofnegQ} * $cos3 - $self->{cubAd3};
                push(@inrange,$root) if (($root > 0 || $root eq 0) && ($root < 1 || $root eq 1));
                }
=cut



            # ifs as blocks instead of at end of push line does not slow things
            # down significantly
            # if conditions hopefully are permissive around the edges - checking
            # if theta is (not outside) the limits rather than (inside || on)
            # Even when checking inside or on, you get false positives that you have to weed out
            # with the if on the push line. So if that has to be there anyway
            # avoid the >||eq slowness, and don't worry about a few false positives getting in.


            my $tocos1=$th;
            #warn "$self->{cubtheta_upper_bound1} > $tocos1 : ",($self->{cubtwotimessqrtofnegQ} * cos($tocos1/3) - $self->{cubAd3})," > $self->{cubtheta_lower_bound1}\n";
            if (1 || !($tocos1 > $self->{cubtheta_upper_bound1}) && !($tocos1 < $self->{cubtheta_lower_bound1})) {
                my $root = $self->{cubtwotimessqrtofnegQ} * cos($tocos1/3) - $self->{cubAd3};
                push(@inrange, $root) if (($root > 0 || $root eq 0) && ($root < 1 || $root eq 1));


                my $deriv =
                    ($self->{cubsqrtofnegQ}/(3*$self->{A}))
                    * (sin($tocos1/3) / (sqrt(-$D)))
                ;
                my $deriv2 =
                    -1 * # because I didn't keep good enough track of the square root of -1 while deriving this
                    ($self->{cubsqrtofnegQ}/(6 * $D * $self->{A}**2))
                    *
                    ( (sin($tocos1/3) / sqrt(-$D))  * $R - cos($tocos1/3)/3 )
                ;



                push(@deriv, $deriv) if $doThetaPrime && (($root > 0 || $root eq 0) && ($root < 1 || $root eq 1));
                push(@deriv2, $deriv2) if $doThetaPrime==2 && (($root > 0 || $root eq 0) && ($root < 1 || $root eq 1));
                }
            my $tocos2=($th + $thistwopi);
            #warn "$self->{cubtheta_upper_bound2} > $tocos2 : ",($self->{cubtwotimessqrtofnegQ} * cos($tocos2/3) - $self->{cubAd3})," > $self->{cubtheta_lower_bound2}\n";
            if (1 || !($tocos2 > $self->{cubtheta_upper_bound2}) && !($tocos2 < $self->{cubtheta_lower_bound2})) {
                my $root = $self->{cubtwotimessqrtofnegQ} * cos($tocos2/3) - $self->{cubAd3};
                push(@inrange, $root) if (($root > 0 || $root eq 0) && ($root < 1 || $root eq 1));


                my $deriv =
                    ($self->{cubsqrtofnegQ}/(3*$self->{A}))
                    * (sin($tocos2/3) / (sqrt(-$D)))
                ;
                my $deriv2 =
                    -1 * # because I didn't keep good enough track of the square root of -1 while deriving this
                    ($self->{cubsqrtofnegQ}/(6 * $D * $self->{A}**2))
                    *
                    ( (sin($tocos2/3) / sqrt(-$D))  * $R - cos($tocos2/3)/3 )
                ;


                push(@deriv, $deriv) if $doThetaPrime && (($root > 0 || $root eq 0) && ($root < 1 || $root eq 1));
                push(@deriv2, $deriv2) if $doThetaPrime==2 && (($root > 0 || $root eq 0) && ($root < 1 || $root eq 1));
                }
            my $tocos3=($th + $thisfourpi);
            #warn "$self->{cubtheta_upper_bound3} > $tocos3 : ",($self->{cubtwotimessqrtofnegQ} * cos($tocos3/3) - $self->{cubAd3})," > $self->{cubtheta_lower_bound3}\n";
            if (1 || !($tocos3 > $self->{cubtheta_upper_bound3}) && !($tocos3 < $self->{cubtheta_lower_bound3})) {
                my $root = $self->{cubtwotimessqrtofnegQ} * cos($tocos3/3) - $self->{cubAd3};
                push(@inrange, $root) if (($root > 0 || $root eq 0) && ($root < 1 || $root eq 1));

                my $deriv =
                    ($self->{cubsqrtofnegQ}/(3*$self->{A}))
                    * (sin($tocos3/3) / (sqrt(-$D)))
                ;
                my $deriv2 =
                    -1 * # because I didn't keep good enough track of the square root of -1 while deriving this
                    ($self->{cubsqrtofnegQ}/(6 * $D * $self->{A}**2))
                    *
                    ( (sin($tocos3/3) / sqrt(-$D))  * $R - cos($tocos3/3)/3 )
                ;
                #warn "\n root deriv: $deriv [$self->{cubtwotimessqrtofnegQ} , $R / $sqrtofnegQcubed , $self->{A}]\n";
                push(@deriv, $deriv) if $doThetaPrime && (($root > 0 || $root eq 0) && ($root < 1 || $root eq 1));
                push(@deriv2, $deriv2) if $doThetaPrime==2 && (($root > 0 || $root eq 0) && ($root < 1 || $root eq 1));
                }

            }
#die;

        }


    if ($doThetaPrime) {

        # returns list of array refs instead of scalars

        my @ret_inds;
        # three item sort, because we think "sort {...} @ret" is expensive
        # hmm for this one - is this still faster than sort? maybe
        splice(@ret_inds,$#ret_inds==-1?0:($inrange[$_]<$inrange[$ret_inds[0]]?0:($inrange[$_]>$inrange[$ret_inds[-1]]?$#ret_inds+1:$#ret_inds)),0,$_)
            for (0 .. $#inrange);
        if ($doThetaPrime == 2) {
            return map {[$inrange[$_], $deriv[$_], $deriv2[$_]]} @ret_inds;
            }
        else {
            return map {[$inrange[$_], $deriv[$_]]             } @ret_inds;
            }
        }
    elsif ($allReal) { # up to three real roots
        my @ret;
        # three item sort, because we think "sort {...} @ret" is expensive
        splice(@ret,$#ret==-1?0:($_<$ret[0]?0:($_>$ret[-1]?$#ret+1:$#ret)),0,$_)
            for @inrange;
        return @ret;
       }
    else { # one real root, when D > 0
        return ($inrange[0]);
        }

    }


sub solveXforThetaBig {
    my $self = shift;
    my $x    = shift;
    # x = At^3 +  Bt^2 +     Ct +     D
    # 0 = At^3 +  Bt^2 +     Ct +     D -x
    # 0 =  t^3 + (B/A)t^2 + (C/A)t + (D-x)/A
    #print "for BIG cubic formula: $self->{BdA_Big}  ,  $self->{CdA_Big}, ",$self->{D_Big}->copy()->bsub($x)->bdiv($self->{A_Big})->bstr,"\n";
    #print "pre filter cubic formula          returns:",join(" , ",&cubicformula($self->{BdA}    ,$self->{CdA}    ,($self->{D}     - ((ref($x))?$x:$x))/$self->{A}    ,0)),"\n";
    #print "pre filter cubic formula with BIG returns:",join(" , ",&cubicformula($self->{BdA_Big},$self->{CdA_Big},($self->{D_Big} - ((ref($x))?$x:$x))/$self->{A_Big},0)),"\n";
    #print "going BIGGGGGGG with x: $x\n";
    my $cubic_solve_C = $self->{D_Big}->copy();
    $cubic_solve_C->bsub($x)->bdiv($self->{A_Big});
    return sort {$a<=>$b} grep {warn "big: ",$_,"\n"; (1 > $_ || 1 eq $_) && ($_ > 0 || $_ eq 0)} &cubicformula($self->{BdA_Big},$self->{CdA_Big},$cubic_solve_C,1);
    }
sub solveYforTheta {
    my $self = shift;
    my $y    = shift;
    # y = Et^3 +  Ft^2 +     Gt +     H
    # 0 = Et^3 +  Ft^2 +     Gt +     H - y
    # 0 =  t^3 + (F/E)t^2 + (G/E)t + (H - y)/E
    return sort {$a<=>$b} grep {(1 > $_ || 1 eq $_) && ($_ > 0 || $_ eq 0)} &cubicformula($self->{FdE},$self->{GdE},($self->{H} - ((ref($y))?$y:$y))/$self->{E},1);
    }
sub solveYforThetaBig {
    my $self = shift;
    my $y    = shift;
    # y = Et^3 +  Ft^2 +     Gt +     H
    # 0 = Et^3 +  Ft^2 +     Gt +     H - y
    # 0 =  t^3 + (F/E)t^2 + (G/E)t + (H - y)/E
    #my @ts = &cubicformula($self->{FdE_Big},$self->{GdE_Big},($self->{H_Big} - ((ref($y))?$y:$y))/$self->{E_Big},1);
    #print " ts: ",join(", ",@ts),"\n";
    return sort {$a<=>$b} grep {(1 > $_ || 1 eq $_) && ($_ > 0 || $_ eq 0)} &cubicformula($self->{FdE_Big},$self->{GdE_Big},($self->{H_Big} - ((ref($y))?$y:$y))/$self->{E_Big},1);
    }
sub solveXPrimeforTheta { # x's slope related to theta (not y)
    my $self = shift;
    my $xp    = shift;
    # x'= 3At^2+  2Bt  +     C
    # 0 = 3At^2+  2Bt  +     C - x'
    return sort {$a<=>$b} grep {(1 > $_ || 1 eq $_) && ($_ > 0 || $_ eq 0)} &quadraticformula($self->{Am3},$self->{Bm2},$self->{C} - ((ref($xp))?$xp:$xp),1);
    }
sub solveXPrimeforThetaBig { # x's slope related to theta (not y)
    my $self = shift;
    my $xp    = shift;
    return sort {$a<=>$b} grep {(1 > $_ || 1 eq $_) && ($_ > 0 || $_ eq 0)} &quadraticformula($self->{Am3_Big},$self->{Bm2_Big},$self->{C_Big} - ((ref($xp))?$xp:$xp),1);
    }

# same as above two functions, but allow results outside 0 to 1 range
sub solveXPrimeforTheta_noFilter {
    my $self = shift;
    my $xp    = shift;
    # x'= 3At^2+  2Bt  +     C
    # 0 = 3At^2+  2Bt  +     C - x'
    return sort {$a<=>$b} &quadraticformula($self->{Am3},$self->{Bm2},$self->{C} - ((ref($xp))?$xp:$xp),1);
    }
sub solveXPrimeforThetaBig_noFilter {
    my $self = shift;
    my $xp    = shift;
    return sort {$a<=>$b} &quadraticformula($self->{Am3_Big},$self->{Bm2_Big},$self->{C_Big} - ((ref($xp))?$xp:$xp),1);
    }


sub solveYPrimeforTheta {
    my $self = shift;
    my $yp    = shift;
    # y'= 3Et^2+  2Ft  +     G
    # 0 = 3Et^2+  2Ft  +     G - y'
    return sort {$a<=>$b} grep {(1 > $_ || 1 eq $_) && ($_ > 0 || $_ eq 0)} &quadraticformula($self->{Em3},$self->{Fm2},$self->{G} - ((ref($yp))?$yp:$yp),1);
    }
sub solveYPrimeforThetaBig {
    my $self = shift;
    my $yp    = shift;
    return sort {$a<=>$b} grep {(1 > $_ || 1 eq $_) && ($_ > 0 || $_ eq 0)} &quadraticformula($self->{Em3_Big},$self->{Fm2_Big},$self->{G_Big} - ((ref($yp))?$yp:$yp),1);
    }
sub solvefprimefortheta { #obsolete? since this is just slopeTangent(x) ?
    my $self = shift;
    my $fp    = shift;
    # f'= y'/x' = (3Et^2 + 2Ft + G) / (3At^2 + 2Bt + C)
    # 0 = 3(E-f'A)t^2 + 2(F-f'B)t + (G-f'C)
    return sort {$a<=>$b} grep {(1 > $_ || 1 eq $_) && ($_ > 0 || $_ eq 0)} &quadraticformula($self->{Em3} - $fp * $self->{Am3},$self->{Fm2} - $fp * $self->{Bm2},$self->{G} - $fp * $self->{C},1);
    }

our $BigFloatOneHalf = Math::BigFloat->new('0.5');
our $BigFloatTen     = Math::BigFloat->new('10');
sub bigsqrt {
    #because the sqrt and root functions in Math::BigFloat sometimes fail
    #Wikipedia says:
    #sqrt(x) = 10**(1/2 * log_10(x))
    return $BigFloatTen->copy()->bpow($BigFloatOneHalf->copy()->bmul($_[0]->copy()->blog(10)),25);
    }

sub dimensionalStepFromTheta {
    my $self=shift;

    my $dim=shift;
    my $theta=shift;
    my $direction=scalar(@_)?shift:1; # 1 or 0

    my $findnexttheta = sub {
        my $ret;
        my $pt_last = $self->point($theta);
        if (ref($pt_last->[0])) {$pt_last->[0]=0 + sprintf("%.20f",$pt_last->[0]->bstr);}
        if (ref($pt_last->[1])) {$pt_last->[1]=0 + sprintf("%.20f",$pt_last->[1]->bstr);}
        if (!ref($_[0])) {
            #warn "\nNORMALLLLLLLLLLLLLLLLLLLLLLLLLLLLLL    \n";
            my $pt_new  = $self->point($_[0]);
            if (ref($pt_new->[0])) {$pt_new->[0]=0 + sprintf("%.20f",$pt_new->[0]->bstr);}
            if (ref($pt_new->[1])) {$pt_new->[1]=0 + sprintf("%.20f",$pt_new->[1]->bstr);}

            $ret = $dim - CORE::sqrt(($pt_new->[0] - $pt_last->[0])**2 + ($pt_new->[1] - $pt_last->[1])**2);
            #print "$ret = $dim - CORE::sqrt(($pt_new->[0] - $pt_last->[0])**2 + ($pt_new->[1] - $pt_last->[1])**2) = $ret\n";
            }
        else {
            warn "I don't think you want to be here - not sure if this mess is debugged.\n";
            my $pt_new  = $self->point($_[0]);
            print "using BigFloat with trial theta = $_[0]\n";
        #    my $dx=(ref($pt_new->[0]))?$pt_new->[0]->copy()->bsub($pt_last->[0]):$pt_new->[0] - $pt_last->[0];
        #    my $dxsqrd=(ref($dx))?$dx->copy()->bpow(2):$dx**2;
        #    my $dy=(ref($pt_new->[1]))?$pt_new->[1]->copy()->bsub($pt_last->[1]):$pt_new->[1] - $pt_last->[1];
        #    my $dysqrd=(ref($dy))?$dy->copy()->bpow(2):$dy**2;
        #    my $distsqrd=(ref($dysqrd))?$dysqrd->copy()->badd($dxsqrd):$dysqrd + $dxsqrd;
        #    my $dist = (ref($distsqrd))?bigsqrt($distsqrd):sqrt($distsqrd);

            my $dxsqrd = ref($pt_new->[0]) ? $pt_new->[0]->copy() : $pt_new->[0];
            if (ref($dxsqrd)) {$dxsqrd->bsub($pt_last->[0]);} else {$dxsqrd -= $pt_last->[0];}
            if (ref($dxsqrd)) {$dxsqrd->bpow(2);} else {$dxsqrd=$dxsqrd**2;}
            my $dysqrd = ref($pt_new->[1]) ? $pt_new->[1]->copy() : $pt_new->[1];
            if (ref($dysqrd)) {$dysqrd->bsub($pt_last->[1]);} else {$dysqrd -= $pt_last->[1];}
            if (ref($dysqrd)) {$dysqrd->bpow(2);} else {$dysqrd=$dysqrd**2;}
            my $distsqrd = ref($dysqrd) ? $dysqrd->copy() : $dysqrd;
            if (ref($distsqrd)) {$distsqrd->badd($dxsqrd)} else {$distsqrd += $dxsqrd;}
            my $dist = (ref($distsqrd))?bigsqrt($distsqrd):sqrt($distsqrd);

            if (ref($dim)) {
                #$ret= $dim - $dist;
                $ret= (0 + sprintf("%.20f",$dim->bstr)) -  (ref($dist) ? (0 + sprintf("%.20f",$dist->bstr)) : $dist);
                }
            else {
                my $dimB = Math::BigFloat->new(''.$dim);
                $ret= $dimB - $dist;
                }
            }
        return $ret;
        };

    my $newtheta;
    my $er;
    #warn "$dim , $theta , $direction , $self->{precision}";
    if (!defined($theta)) {warn "theta undefined";}
    ($newtheta,$er) = FalsePosition($findnexttheta,($direction ? [$theta,1]:[0,$theta]),$self->{precision},($direction ? ($theta + (1-$theta)/2):($theta/2)),'dimensionalStepFromTheta for Bezier segment');
    #print " \nAFTR ROOT FIND: $newtheta , REF?:",(ref($newtheta) ? 'yes':'no'),"\n";
    if (length($newtheta)>17) {$newtheta=sprintf("%.15f",$newtheta);} #FalsePosition will return number as string if it ended up with BigFloat and none of inputs were BigFloat - so eval() here to turn that big string into a perl-sized number
    #print " dim step result ($newtheta,$er)\n";
    if (defined($er)) {
        warn "dimstep er: $er";
        #probably just reached the end
        if (abs(&{$findnexttheta}(($direction ? 1:0))) < $dim) {
            warn "snapping to end, I think, in dimensionalStepFromTheta";
            $newtheta=($direction ? 1:0);
            }
        #otherwise the error might be real
        else {warn "error in dimensionalStepFromTheta???? $er"}
        }
    my $r=$self->point($newtheta);
    if (ref($r->[0])) {$r->[0]=0 + sprintf("%.20",$r->[0]->bstr);}
    if (ref($r->[1])) {$r->[1]=0 + sprintf("%.20",$r->[1]->bstr);}
    return ($r,$newtheta);
    }

sub asin {
    #Based on Wolfram MathWorld
    #http://mathworld.wolfram.com/InverseSine.html
    #but it's supposed to use complex numbers, hmmm
    if    ($_[0] eq -1) {return $pi/2 * -1} #eq 2
    elsif ($_[0] eq  1) {return $pi/2     } #eq 4
    elsif ($_[0] eq  0) {return 0;        } #eq 3
    else {                                  #eq 15
# SPEED ATTEMPTS - why was this eval() here?
#        if (eval(abs($_[0]))>1) {warn("giving something bigger than +/-1 to your asin:",$_[0],"\n");}
        if (abs($_[0])>1) {warn("giving something bigger than +/-1 to your asin:",$_[0],"\n");}
        return atan2($_[0],sqrt(1-$_[0]**2));
        }
    }
sub acos {return ($pi/2) - asin($_[0]);}

}
1;
