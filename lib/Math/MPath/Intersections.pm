####################################################################################
###      Math::MPath::Intersections    #############################################
package Math::MPath::Intersections;
{
use Math::MPath::Function::Root qw(BrentsMethod FalsePosition);

our $pi = 4 * atan2(1,1);

sub bez_bez_intersect {

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
            # 

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
                                
                my @int1 = bez_bez_intersect($bezA,$bezB,
                                             [[$spanA->[0],[$spanx->[0]  , $span_split_x], $spanA->[3]?$sub_t_span_2_A:$sub_t_span_1_A, $spanA->[3]]],
                                             [[$spanB->[0],[$spanx->[0]  , $span_split_x], $spanB->[3]?$sub_t_span_2_B:$sub_t_span_1_B, $spanB->[3]]]
                                             );
                
                my @int2 = bez_bez_intersect($bezA,$bezB,
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


sub bezoff_bezoff_intersect {

    # The first attempt at adapting bez_bez_intersect() for offset beziers

    # returns list of intersections and the corresponding Bezier parameters
    # for each input Bezier offset curve.

    # lutA and lutB optional - used to re-call this function in the intersection

    my ($bezA, $bezB, $offA, $offB, $lutA, $lutB) = @_;

    my $XtoTLUT_A = $lutA ? $lutA : [map [$_->[0],[$_->[1]->[0],$_->[1]->[-1]],[$_->[2]->[0],$_->[2]->[-1]],$_->[3]], @{$bezA->{XtoTLUT}}];
    my $XtoTLUT_B = $lutB ? $lutB : [map [$_->[0],[$_->[1]->[0],$_->[1]->[-1]],[$_->[2]->[0],$_->[2]->[-1]],$_->[3]], @{$bezB->{XtoTLUT}}];

    if (!defined($lutA)) { # assume LUTs already contain offset Xs in their x range entries
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
                # 

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
                    
                    my $tA = $bezA->t_from_xoff($_[0],$offA,[$spanA->[2]->[0],$spanA->[2]->[-1]],$spanA->[0]->[1],$spanA->[3]);
                    my $tB = $bezB->t_from_xoff($_[0],$offB,[$spanB->[2]->[0],$spanB->[2]->[-1]],$spanB->[0]->[1],$spanB->[3]);
    
                    my $nA = $bezA->F_prime(undef,$tA);
                    my $nB = $bezB->F_prime(undef,$tB);
    
                    my $nA_prime = $bezA->F_2prime(undef,$tA);
                    my $nB_prime = $bezB->F_2prime(undef,$tB);
    
                    my $YPrimeA = $bezA->bezierEvalYPrimeofT($tA);
                    my $YPrimeB = $bezB->bezierEvalYPrimeofT($tB);
    
                    $yoffset_primeA = -($offA/2.0) * 1/(sqrt($nA**2 + 1)**3) * 2*$nA * $nA_prime;
                    $yoffset_primeB = -($offB/2.0) * 1/(sqrt($nB**2 + 1)**3) * 2*$nB * $nB_prime;
    
                    my $YoffPrimeA = ($YPrimeA + $yoffset_primeA);
                    my $YoffPrimeB = ($YPrimeB + $yoffset_primeB);
                    
                    my $ret = $YoffPrimeB - $YoffPrimeA;
    
                    return $ret;
    
                };
    
                my $at_start = $y_diff_prime->($spanx->[0]);
                my $at_end   = $y_diff_prime->($spanx->[1]);
    
                #warn "ydiffprime start, end: [$spanx->[0]] $at_start, [$spanx->[1]] $at_end\n";
    
                if ($at_start > 0 && $at_end < 0 || $at_start < 0 && $at_end > 0) {
    
                    my $bounds_xoffA = [$spanx->[0] + 1.00001 * $e,$spanx->[1] - 1.00001 * $e];
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
                    my @int1 = bezoff_bezoff_intersect($bezA,$bezB,$offA,$offB,
                                                 [[$spanA->[0],[$spanx->[0]  , $span_split_x], $spanA->[3]?$sub_t_span_2_A:$sub_t_span_1_A, $spanA->[3]]],
                                                 [[$spanB->[0],[$spanx->[0]  , $span_split_x], $spanB->[3]?$sub_t_span_2_B:$sub_t_span_1_B, $spanB->[3]]]
                                                 );
                    #warn "re-call 2\n";
                    #warn "  [[$spanA->[0],[$span_split_x, $spanx->[1]  ], $spanA->[3]?$sub_t_span_1_A:$sub_t_span_2_A, $spanA->[3]]]\n";
                    #warn "  [[$spanB->[0],[$span_split_x, $spanx->[1]  ], $spanB->[3]?$sub_t_span_1_B:$sub_t_span_2_B, $spanB->[3]]]\n";
    
                    my @int2 = bezoff_bezoff_intersect($bezA,$bezB,$offA,$offB,
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

} # end package
1;