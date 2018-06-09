####################################################################################
###      Math::MPath::MoveTo        #############################################
package Math::MPath::MoveTo;
{
push @Math::MPath::MoveTo::ISA, 'Math::MPath::LineSegment';

sub getLength {return 0;} #//per SVG spec
sub getFeet {return ();} #//per SVG spec
sub getIntersections {return ();} #//per SVG spec

}
####################################################################################
###      Math::MPath::ClosePath     #############################################
package Math::MPath::ClosePath;
{
push @Math::MPath::ClosePath::ISA, 'Math::MPath::LineSegment';
}
####################################################################################
###      Math::MPath::LineSegment   #############################################
package Math::MPath::LineSegment;
{
use Math::BigFloat;
use Math::MPath::Function::Root qw(FalsePosition); #get rid of this when/if you rewrite the dimensionalStepFromTheta sub
use Carp qw(cluck);


sub new {
    my $class = shift;
    my $self={};
    bless $self,$class;
    $self->{p1} = shift;
    $self->{p2} = shift;
    $self->{precision} = shift;
    $self->{isLite} = @_?shift:0;
    ($self->{minx},$self->{maxx}) = sort {$a<=>$b} ($self->{p1}->[0],$self->{p2}->[0]);
    ($self->{miny},$self->{maxy}) = sort {$a<=>$b} ($self->{p1}->[1],$self->{p2}->[1]);
    $self->{maxdiaglength} = sqrt(($self->{maxx} - $self->{minx})**2 + ($self->{maxy} - $self->{miny})**2);
    $self->{m}  = ($self->{p2}->[0] - $self->{p1}->[0] == 0)?(($self->{p2}->[1] < $self->{p1}->[1])?'-':'').'inf':($self->{p2}->[1] - $self->{p1}->[1])/($self->{p2}->[0] - $self->{p1}->[0]);
    $self->{b}  = $self->{p1}->[1] - $self->{m} * $self->{p1}->[0];
    $self->{slopeTangent}=$self->{m};
    $self->{slopeNormal}= ($self->{slopeTangent} == 0)? (($self->{p2}->[0]<$self->{p1}->[0])?'-':'').'inf':-1/$self->{slopeTangent};
    $self->{dx}=$self->{p2}->[0] - $self->{p1}->[0];
    $self->{dy}=$self->{p2}->[1] - $self->{p1}->[1];
    $self->{angleTangent}=atan2($self->{dy},$self->{dx});
    $self->{angleNormal}=atan2(-$self->{dx},$self->{dy});
    $self->{length} = sqrt(($self->{p2}->[0] - $self->{p1}->[0])**2 + ($self->{p2}->[1] - $self->{p1}->[1])**2);
    $self->{length_big} = undef; # BigFloat used in point(theta) if theta is a ref. Create there only if needed.
    return $self;
}

sub getInfiniteSlopeThetas {
    my $self = shift;
    if ($self->{dx}==0) {return wantarray ? (0,1):0;}
    else {return;}
}

sub getLength {
    my $self = shift;
    my $res = shift;
    my $start_theta = shift;
    my $end_theta = shift;
    if (!defined($res)) {$res=1000;}
    if (!defined($start_theta)) {$start_theta=0;}
    if (!defined($end_theta)) {$end_theta=1;}
    return $self->{length} * abs($end_theta - $start_theta);
}
sub precision {
    my $self = shift;
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
    if (
          defined($coords->[0]) &&
          (
             $self->{minx} < $coords->[0] ||
             $self->{minx} eq $coords->[0]
          ) &&
          (
             $self->{maxx} > $coords->[0] ||
             $self->{maxx} eq $coords->[0]
          )
        ) {$xok=1;}
    if (
          defined($coords->[1]) &&
          (
             $self->{miny} < $coords->[1] ||
             $self->{miny} eq $coords->[1]
          ) &&
          (
             $self->{maxy} > $coords->[1] ||
             $self->{maxy} eq $coords->[1]
          )
        ) {$yok=1;}
    return $xok,$yok;
}
sub onSegment {
    my $self = shift;
    my $point= shift;
    if ($point->[1] - ($self->{m} * $point->[0] + $self->{p1}->[0]) < $self->{precision}) {return 1;}
    else {return 0;}
}
sub getFeet {
    my $self = shift;
    my $x = shift;
    my $y = shift;
    my @feet=();
    if ($self->{m} eq 0) {
        if (($self->inRange([$x,undef]))[0]) { push(@feet, [$x,$self->{p1}->[1], $self->solveXforTheta($x)] );}
    }
    elsif ($self->{m} =~ /inf/) {
        if (($self->inRange([undef,$y]))[1]) { push(@feet, [$self->{p1}->[0],$y, $self->solveYforTheta($y)] );}
    }
    else {
        my $intersect_x = (($self->{m}*$self->{p1}->[0])-($self->{p1}->[1])+((1/$self->{m})*$x)+($y))/($self->{m} + (1/$self->{m}));
        my $bx;
        my $by;
        ($bx,$by) = $self->inRange([$intersect_x,undef]);
        if ($bx) {
            my $intersect_y = $self->f($intersect_x);
            ($bx,$by) = $self->inRange([undef,$intersect_y]);
            if ($by) {
                push(@feet,[$intersect_x,$intersect_y, $self->solveXforTheta($intersect_x)]);
            }
        }
    }
    return @feet;
}
sub f {
    my $self = shift;
    my $x = shift;
    if ($self->{minx} > $x || $self->{maxx} < $x) {return;}
    elsif ($self->{m}=~/(\-?)inf/) {
        if (!wantarray) {return $self->{maxy};} #fixes my immediate problem, maybe, but what's the best policy for this whole inf slope thing?
        my $n = (0 + ($1.'1'));
        #this gives too many, and it's hard to figure a good general purpose way to come up with a number here
        #my $stepcount = abs(($self->{p2}->[1] - $self->{p1}->[1])/($self->{precision} * 10));
        #print "stepcount: $stepcount = abs(($self->{p2}->[1] - $self->{p1}->[1])/($self->{precision} * 10));\n";
        #so lets say 20 for now and revisit next time we're stuck on this situation
        my $stepcount=20;
        my $wholestepcount=int($stepcount);
        my $remainder = $stepcount - $wholestepcount;
        my $step = $n * $self->{precision};
        my @ret=();
        foreach (0 .. $wholestepcount) {push(@ret,$self->{p1}->[1] + $step * $_);}
        push(@ret,$self->{p1}->[1] + $step * $wholestepcount + $step * $remainder);
        return @ret;
    }
    elsif ($self->{m} eq 0) {return $self->{p1}->[1];}
    elsif ($x eq $self->{p1}->[0]) {return $self->{p1}->[1];}
    elsif ($x eq $self->{p2}->[0]) {return $self->{p2}->[1];}
    else {return $self->{m} * $x + $self->{b};}
}
sub F {
    my $self = shift;
    my $y = shift;
    if ($self->{miny} > $y || $self->{maxy} < $y) {return;}
    if ($self->{m} eq 0) {
        if (!wantarray) {return $self->{maxx};} #fixes my immediate problem, but what's the best poicy for this whole inf slope thing?
        my $n = ($self->{p2}->[0] > $self->{p1}->[0])?1:-1;
        my $stepcount = abs(($self->{p2}->[0] - $self->{p1}->[0])/10); #hmm
        my $wholestepcount=int($stepcount);
        my $remainder = $stepcount - $wholestepcount;
        my $step = $n * $self->{precision};
        my @ret=();
        foreach (0 .. $wholestepcount) {push(@ret,$self->{p1}->[0] + $step * $_);}
        push(@ret,$self->{p1}->[0] + $step * $wholestepcount + $step * $remainder);
        return @ret;
    }
    elsif ($self->{m}=~/(\-?)inf/) {return $self->{p1}->[0];}
    else {return ($y - $self->{b})/$self->{m};}
}
sub solveXforTheta {
    my $self=shift;
    my $x=shift;
    my $th = abs(($x - $self->{p1}->[0])/($self->{p2}->[0] - $self->{p1}->[0]));
    if ((1 > $th || 1 eq $th) && ($th > 0 || $th eq 0)) {return ($th);}
    else {return ();}
}
sub solveYforTheta {
    my $self=shift;
    my $y=shift;
    if ($self->{p2}->[1] - $self->{p1}->[1] eq 0) {
        # TODO: handle
        warn "[$self->{p1}->[0] , $self->{p1}->[1]]\n[$self->{p2}->[0] , $self->{p2}->[1]]";
    }
    my $th = abs(($y - $self->{p1}->[1])/($self->{p2}->[1] - $self->{p1}->[1]));
    if ((1 > $th || 1 eq $th) && ($th > 0 || $th eq 0)) {return ($th);}
    else {return ();}
}
sub point {
    my $self = shift;
    my $theta = shift;
    my ($x,$y);
    if ($theta==0)    { $x=$self->{p1}->[0];$y=$self->{p1}->[1];}
    elsif ($theta==1) { $x=$self->{p2}->[0];$y=$self->{p2}->[1];}
    else {
        # At some point I started doing all line segment point lookups
        # with BigFloats. It didn't hurt too much at the time, and it
        # allowed the 0-to-1 paramaterization to work with longer
        # lines when doing fine-grained root finding.
        # But later I was processing lots of short lines
        # and this became a bottleneck. So, back to normal PERL-sized
        # numbers.
        #
        # What's needed is a good way to choose when to use Bigs.
        #
        # For now I'll say, if you pass in a BigFloat theta, you'll
        # get BigFloat processing and a BigFloat response.
        #
        if (!ref($theta)) { # if theta is not a reference to a BigFloat, just do the normal simple math with normal numbers
            $x=($self->{dx} * $theta) + $self->{p1}->[0];
            $y=($self->{dy} * $theta) + $self->{p1}->[1];
        }
        else {
            if (!defined $self->{length_big}) {$self->{length_big} = new Math::BigFloat ''.$self->{length};}
            my $partdist=$self->{length_big}->copy()->bmul($theta);
            $x=$partdist->copy()->bmul(CORE::cos($self->{angleTangent}));
            $y=$partdist->copy()->bmul(CORE::sin($self->{angleTangent}));
            $x->badd($self->{p1}->[0]);
            $y->badd($self->{p1}->[1]);
            if (!ref($theta)) {
                if (ref($x)) {$x=0 + sprintf("%.20f",$x->bstr());}
                if (ref($y)) {$y=0 + sprintf("%.20f",$y->bstr());}
            }
        }
    }
    return [$x,$y];
}
sub point_offset {
    my ($self, $t, $distance) = @_;
    my ($x,$y) = @{$self->point($t)};
    $x += -($self->{dy}/$self->{length}) * $distance;
    $y +=  ($self->{dx}/$self->{length}) * $distance;
    return [$x,$y];
}

sub secondDerivative {return 0;}
sub slopeTangent {return $_[0]->{slopeTangent};}
sub slopeTangent_byTheta {return $_[0]->{slopeTangent};}
sub fDoublePrime {return 0};
sub slopeNormal  {return $_[0]->{slopeNormal};}
sub slopeNormal_byTheta  {return $_[0]->{slopeNormal};}
sub angleTangent {return $_[0]->{angleTangent};}
sub angleNormal  {return $_[0]->{angleNormal};}
sub angleTangent_byTheta {return $_[0]->{angleTangent};}
sub angleNormal_byTheta  {return $_[0]->{angleNormal};}

our $BigFloatOneHalf = Math::BigFloat->new('0.5');
our $BigFloatTen     = Math::BigFloat->new('10');
sub bigsqrt {
    #because the sqrt and root functions in Math::BigFloat sometimes fail
    #Wikipedia says:
    #sqrt(x) = 10**(1/2 * log_10(x))
    return $BigFloatTen->copy()->bpow($BigFloatOneHalf->copy()->bmul($_[0]->copy()->blog(10)),25);
}

sub dimensionalStepFromTheta {


    # Oh, things don't run as fast as you would like?
    # Well, why don't you rewrite this for the simple line segment case, instead of using
    # the copy-pasted code from CubicBezier - duh.
    # (actually this is pretty fast - not likely a bottleneck - still, rewrite it sometime)
    # and it looks like I did this first for all the javascript versions
    # including arc, so if it's good enough for javascript, maybe it's
    # good enough for Perl


    my $self = shift;

    my $dim = shift;
    my $theta = shift;
    my $direction = scalar(@_) ? shift : 1;

    my $pt_last = $self->point($theta);

    my $findnexttheta = sub {
        my $ret;
        if (!ref($_[0])) {
            my $pt_new  = $self->point($_[0]);
            $ret = $dim - CORE::sqrt(($pt_new->[0] - $pt_last->[0])**2 + ($pt_new->[1] - $pt_last->[1])**2);
        }
        else {
            #warn "I don't think you want to be here - not sure if this mess is debugged.\n";
            my $pt_new  = $self->point($_[0]);
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
    ($newtheta,$er) = FalsePosition($findnexttheta,($direction ? [$theta,1]:[0,$theta]),$self->{precision},($theta+($direction ? 1:0))/2,'dimensionalStepFromTheta for Line segment');
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
