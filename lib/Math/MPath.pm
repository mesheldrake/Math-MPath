package Math::MPath;
{

use Math::MPath::BezierCubicSegment;
use Math::MPath::BezierQuadraticSegment;
use Math::MPath::LineSegment;
use Math::MPath::EllipticalArc;

######################### start copy-paste of old MikePath guts
### Please edit and improve drastically.
### Keep old MikePath API when convenient, but don't be afraid to impove, change, extend.

use Math::MPath::CubicFormula; # Eliminate. Used once in intersection stuff. Smarter Math::MPath::BezierCubicSegment should make this obsolete.
use Math::MPath::QuadraticFormula; #Eliminate. Used once for line-ellipse intersection. make minimal optimized version and stuff in here, or relegate to EllipticalArc code.
use vars qw($VERSION);
use strict;
use warnings;
#use LWP::Simple;
use Carp qw(cluck croak);
use Math::MPath::Function::Root qw(BrentsMethod FalsePosition);
use POSIX qw();
$VERSION = '0.01';

our $enableCarefulFofy=1;
our $enableCarefulfofx=1;

sub mergePaths {
	my @segs;
	push(@segs,$_[0]->{pathSegmentSpecs}->[0]);
	foreach (@_) {
		push (@segs,@{$_->{pathSegmentSpecs}}[1 .. scalar(@{$_->{pathSegmentSpecs}}) - 1]);
		}
	my $ret = new Math::MPath(join('',@segs),$_[0]->{resolution},$_[0]->{precision});
	#foreach (@_) {undef($_);}
	return $ret;
	}
sub new {

	my $class = shift;
	my $self={};
	bless $self,$class;
#	if ($_[0] =~ /^((file|http|ftp)\:.*)/) {
#		$_=get($1);die "Couldn't get file: $1" unless defined $_;
#		s/^#.*\r?\n//g;s/\r?\n//g;
#		push(@_,$_);
#		}
#	else {
		$self->{pathspec} = shift;
#		}
	#print "  making : ",$self->{pathspec},"\n";
	$self->{resolution} = shift;
	$self->{precision} = @_?shift:$self->{resolution}/1000;
	$self->constructSegments($self->{pathspec});
	# not used anymore? at least at this composit level $self->{thetaprecision} = $self->{precision}/$self->{xmax};
	#$self->{length} = getLength($self);
	#my $toff=0;
	#$self->{pathThetaToCompThetaRanges}=[];
	#for (my $i=0;$i<scalar(@{$self->{pathSegments}});$i++) {
	#	$self->{pathThetaToCompThetaRanges}->[scalar(@{$self->{pathThetaToCompThetaRanges}})] = [$toff,$toff + ($self->{pathSegments}->[$i]->{length}/$self->{length})];
	#	$toff += ($self->{pathSegments}->[$i]->{length}/$self->{length});
	#	#now the indexes of those ranges should match the indexes (indeces) of the components
	#	}
	#$self->{pathThetaToCompThetaRanges}->[scalar(@{$self->{pathThetaToCompThetaRanges}})-1]->[1]=1;
	return $self;
	}
sub newlite {
	my $class = shift;
	my $self={};
	bless $self,$class;
	$self->{isLite}=1;
#	if ($_[0] =~ /^((file|http|ftp)\:.*)/) {
#		$_=get($1);die "Couldn't get file: $1" unless defined $_;
#		s/^#.*\r?\n//g;s/\r?\n//g;
#		push(@_,$_);
#		}
#	else {
		$self->{pathspec} = shift;
#		}
	#print "  making : ",$self->{pathspec},"\n";
	$self->{resolution} = shift;
	$self->{precision} = @_?shift:$self->{resolution}/1000;
	$self->constructSegments($self->{pathspec});
	my $toff=0;
	return $self;
	}
sub constructSegments {
	my $self = shift;
	$self->{pathspec} = shift;
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

# set these up so they get handled like linetos
#warn " at $i $thisSegSpec\n";
if ($thisSegSpec=~/^H/i) {
    $thisSegSpec=~s/\s//g;
    if ($i==1) {
        $lastM=~/M\s*?[0-9\-\.eE]+\s*?[\s,]\s*?([0-9\-\.eE]+)/;
        $thisSegSpec.=','.$1;
        }
    else {$thisSegSpec.=','.$self->{pathSegments}->[$i - 2]->{p2}->[1];}
    $thisSegSpec=~s/^H/L/;
#print "h: $thisSegSpec\n";
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
#print "v: $thisSegSpec\n";
    }


		$self->{pathSegments}->[$i - 1] = $self->constructSegment($thisSegSpec,$lastSegSpec);
		$lastSegSpec=$thisSegSpec;
		}
	$self->{minx} = (        sort {$a<=>$b} map {$_->{minx}} @{$self->{pathSegments}})[0];
	$self->{maxx} = (reverse sort {$a<=>$b} map {$_->{maxx}} @{$self->{pathSegments}})[0];
	$self->{miny} = (        sort {$a<=>$b} map {$_->{miny}} @{$self->{pathSegments}})[0];
	$self->{maxy} = (reverse sort {$a<=>$b} map {$_->{maxy}} @{$self->{pathSegments}})[0];

	$self->{pathspecparts}=[($self->{pathspec}=~/([MmZzLlHhVvCcSsQqTtAa]|[0-9\-\.e]+[, ]+[0-9\-\.e]+)[, ]*/g)];#might be useful for an interpolation function, to avoid some regex's
	#print "pathspec parts: ",join("  ,  ",@{$self->{pathspecparts}}) , "\n";
	$self->{pathspecparts} = [(map {if ($_=~/([0-9\-\.e]+)[, ]+([0-9\-\.e]+)/) {[0,[$1,$2]]} else {[1,$_]}} @{$self->{pathspecparts}})];
	#print "pathspec parts: ",join("  ,  ",map {ref($_->[1])?"$_->[1]->[0],$_->[1]->[1]":$_->[1]} @{$self->{pathspecparts}}) , "\n";

	#print "minx:$self->{minx},maxx:$self->{maxx},miny:$self->{miny},maxy$self->{maxy}\n";
	}
sub parameterize {
	my $self=shift;
	my $newlength=0;
	my $toff=0;
	foreach (@{$self->{pathSegments}}) {$_->{length}=$_->getLength();$newlength+=$_->{length};}
	$self->{length} = $newlength;
    if ($newlength == 0) {croak "zero length for path: ",$self->{pathspec};}
	$self->{pathThetaToCompThetaRanges}=[];
	for (my $i=0;$i<scalar(@{$self->{pathSegments}});$i++) {
		$self->{pathThetaToCompThetaRanges}->[scalar(@{$self->{pathThetaToCompThetaRanges}})] = [$toff,$toff + ($self->{pathSegments}->[$i]->{length}/$self->{length})];
		$toff += ($self->{pathSegments}->[$i]->{length}/$self->{length});
		#now the indeces of those ranges should match the indeces of the components
		}
	$self->{pathThetaToCompThetaRanges}->[scalar(@{$self->{pathThetaToCompThetaRanges}})-1]->[1]=1;
	}
sub getSegsInRange {
	my $self = shift;
	my $point= shift;
	return grep {my ($bx,$by) = $_->inRange($point);$bx || $by} @{$self->{pathSegments}};
	}
sub precision {
	my $self=shift;
	my $np=@_?shift:undef;
	if (defined($np)) {
		$self->{precision}=$np;
		#$self->{thetaprecision} = $self->{precision}/$self->{xmax};
		foreach (@{$self->{pathSegments}}) {
			$_->{precision}=$self->{precision};
			#$_->{thetaprecision}=$self->{thetaprecision};
			}
		}
	return $self->{precision};
	}

sub getLength {

	# NOTE - you now support paths with "M" moveto commands
	# per SVG spec, MoveTo does not contribute to path length
	# new MoveTo package ISA LineSegment, but it's getLength
	# is overridden to return zero.

	my $self=shift;
	my $res = shift;
	if (!defined($res)) {$res=1000;}
	my $length=0;
	foreach (@{$self->{pathSegments}}) {
		$length += $_->getLength($res);
		}
	return $length;
	}

sub dimensionalCurve {
	my $self     = shift;
	my $step     = shift; # constant (scalar) or reference to function of distance along curve - f(t), t:[0.0 to 1.0]
	my $startstep= @_?shift:0; # fraction of first step to offset first point from start
	my $endstep  = @_?shift:undef; # fraction of last step to offset last point from end
	                               # for now, only used if $step is a constant
	my $segrange = @_?shift:[0,scalar(@{$self->{pathSegments}}) - 1]; #2 element array of start and stop indeces of path segments to work over

	# developed in and copied from javascript
    # so keep that synced with any changes here

	my @ret=(); # list of points like [x,y,normal,theta] (theta is (normalized, i.e. fraction of distance) theta on whole path, not from segs)

	my $segslength=0;
	for (my $i=$segrange->[0];$i<=$segrange->[1];$i++) {
        # default theta divisor (sample count) is 1000,
        # but that's not enough it looks like at board scale -
        # determined while working on rail piece design
	    $segslength+=$self->{pathSegments}->[$i]->getLength(10000);
	    }
	my $lengthsum=0;
    if ($segslength < 1/10000) {return ();} # quick way to avoid divide by zero
	my $lengthfraction=$lengthsum/$segslength;
	my $stepfunc=(ref($step))?$step:sub {return $step};

	my $stepval=&{$stepfunc}($lengthfraction);

    # default is first step will be zero-length, giving start point of start segment
    # TODO otherwise, seems like if stepfunc is not constant, any $startstep fraction
    # needs to take into account non-linear step growth... maybe a root find thing and or weighted value thing
    # or maybe the factor multiplies the variable of the step function, and not the span of two evals of that function
    # so see theres stuff to figure here.
	my $firststep = defined($startstep) ? $startstep*$stepval : 0;
	my $laststep;# undef is significant, because "0" needs to be a true last step request - actually we don't make that test anywhere below though

    # Check if the whole path is too short to subdivide by the requested dimension.
	if ($segslength<$firststep) {warn "path too short ($segslength) to divide by first step $firststep\n";return ();}

    # If the step size is fixed, and we want to hit a specified end point,
    # adjust step size to be smaller than requested,
    # to hit that endpoint in a whole number of steps.
	if (!ref($step) && defined $endstep) {
		$laststep = (1-$endstep)*$stepval;
		my $fitlength=$segslength - ($firststep + $laststep);
		my $d = $fitlength/$step;
		if ($d<1) { # Path too short to support both the requested first step offset and last step offset
		    $stepfunc = sub {return $step}; # so this will overshoot on the first step, and result will just have the first step point
		    }
		else {
            my $floored_d = POSIX::floor($d);
			my $r = $d - $floored_d;
            my $newstep = $step;
            if ($r != 0) {
                my $adj = ($step*(1.0 - $r))/($floored_d + 1.0); # because we're shrinking step, one more will fit in
                $newstep = $step - $adj;
                }
			$stepfunc = sub {return $newstep};



            #// a*adjstd + b*adjstd + c*adjstd = segslength, c is number of steps we'll end up with
            #// (a + b + c) * adjstd = segslength
            #// (AB + c) * adjstd = segslength
            #// a and b are known
            #// c can be one more than what you get when you figure this with unadjusted step size
            #// call that C
            #// so,
            #// adjstd = segslength/(A+B+C)
            #// yep. looks reasonable.

            my $C = (POSIX::floor($d) + 1);
            my $ABC = ( $startstep + $endstep + $C );
            my $adjustedstep = $segslength / $ABC;
            $stepfunc = sub {return $adjustedstep};
            #// console.log('adjusted: ' + step + ' to ' + adjustedstep);
            $stepval = &{$stepfunc}($lengthfraction);
            $firststep = $startstep * $adjustedstep;
            $laststep = (1 - $endstep) * $adjustedstep; #// laststep isn't used, is it?


			}
		}
    # Otherwise, if the step size is fixed, and no end point offset is specified
    # use exactly requested step size, and don't care where the last point lands,
    # as long as it's less than stepsize from end of last segment.

    # TODO: (this todo is duplicated in both perl and javascript versions)
    # variable (assumed) step size, so need to smart-figure what first step is, and for last step
    # that depends both on first step and step function, so weirder to figure.
    # maybe analogous to fixed step case, you can find a tiny number to subtract from each input to the eval of the step function
    # to shrink the resulting steps (unless the function's weird and that makes them grow? detect that?)
    # If you pull that off, your rail designer will work like you want, more, sorta.

    # The key feature of that would be that you could do a 50% step on that first step
    # and then a mirror version of the path, mirrored through start point, would
    # have a "natural" looking space between it's first point and the orig first point.
    # "Natural" in view of the way the step function is varying steps in that region.

	#elsif (ref($step) && (defined $endstep || defined $endstep)) {}

	my $zeropt= [@{$self->{pathSegments}->[$segrange->[0]]->point(0)},
	            $self->{pathSegments}->[$segrange->[0]]->angleNormal_byTheta(0),
	            $lengthfraction
	            ];

	if ($stepval == $firststep || (!defined($startstep) || $startstep eq '0')) {push(@ret,$zeropt);}

	my $endpt; # set to each seg endpoint in for loop,
	           # and we might use the last one, after the for loop

	for (my $i=$segrange->[0];$i<=$segrange->[1];$i++) {
		my $prevtheta=0; # seg theta, not path theta
		my $firstpt=[@{$self->{pathSegments}->[$i]->point($prevtheta)},
		            $self->{pathSegments}->[$i]->angleNormal_byTheta($prevtheta),
		            $lengthfraction
		            ];

		my @thispointandtheta=($firstpt,0);
		my @prevptnth=@thispointandtheta;
		$endpt=[@{$self->{pathSegments}->[$i]->point(1)},
		       $self->{pathSegments}->[$i]->angleNormal_byTheta(1),
		       1
		       ];

        # downgrade any BigFloats
#        $_ = sprintf("%.15f",$firstpt->[0]->bstr())
        $_ = sprintf("%.15f",$_->bstr())
            for grep ref($_),
                (@$firstpt,@$endpt);

        my $dist_to_end = sqrt(
                       ($endpt->[1] - $thispointandtheta[0]->[1])**2
                     + ($endpt->[0] - $thispointandtheta[0]->[0])**2
        );

        my $first_while_loop = 1;

		while ( $dist_to_end > $stepval ) {

			my $thisstep=$stepval;
			if ($firststep) {$thisstep=$firststep;$firststep=0;}

            # downgrade any BigFloats
            $_ = sprintf("%.15f",$firstpt->[0]->bstr())
                for grep ref($_),
                    (@{$thispointandtheta[0]},$thispointandtheta[1]);

			@prevptnth=@thispointandtheta;

            warn("ref in dimensionalCurve") for grep ref($_), ($thisstep,$prevtheta);

			@thispointandtheta=$self->{pathSegments}->[$i]->dimensionalStepFromTheta($thisstep,$prevtheta,1);

            # downgrade any BigFloats
#            $_ = sprintf("%.15f",$firstpt->[0]->bstr())
            $_ = sprintf("%.15f",$_->bstr())
                for map {warn "saw Big\n";$_}
                    grep ref($_),
                    (@{$thispointandtheta[0]},$thispointandtheta[1]);

			$lengthsum+=sqrt(($thispointandtheta[0]->[0]-$prevptnth[0]->[0])**2 + ($thispointandtheta[0]->[1]-$prevptnth[0]->[1])**2);
			$lengthfraction=$lengthsum/$segslength;
			$stepval=&{$stepfunc}($lengthfraction);

			push(@{$thispointandtheta[0]},$self->{pathSegments}->[$i]->angleNormal_byTheta($thispointandtheta[1]),$lengthfraction);

            # downgrade any BigFloats
            if (ref($thispointandtheta[0]->[2])) {
                $thispointandtheta[0]->[2]=sprintf("%.15f",$thispointandtheta[0]->[2]->bstr());
                }

			push(@ret,$thispointandtheta[0]);

            # avoid infinite loop when theta isn't changing
            if ($prevtheta == $thispointandtheta[1]
                # A same-theta result of 0 can be valid on the first run through
                # this while loop, if our previous point left us very very close
                # to where this step (which would be a $firststep) being right
                # on this new seg's start point. So that's what this test is for:
                && !($prevtheta == 0 && $first_while_loop)
               ) {
                warn "hit LAST to avoid infinite loop when theta isn't changing\nprob should be die?\n";
                last;
                }

			$prevtheta = $thispointandtheta[1];

            $dist_to_end = sqrt(
                           ($endpt->[1] - $thispointandtheta[0]->[1])**2
                         + ($endpt->[0] - $thispointandtheta[0]->[0])**2
            );

            $first_while_loop = 0;
			}

        # TODO: actually hitting every segment endpoint might be a usful
        #       option to stick right here. You _don't_ want that for what
        #       you're using this for now, mostly, like getting nicely spaced
        #       locations for radial rail pieces. But it would be good for
        #       other things where you need to capture inflection points where
        #       segments meet along a path.

		# While loop ended because next step would have stepped beyond
		# end of current segment. Now set first step for next segment,
		# accounting for whatever partial step distance was taken up by
		# the end of this segment.
		$firststep=$stepval - $dist_to_end;

		}

    # The while loop should _not_ get the last step in this case.
    # (though rarely it might? giving duplicate to this point?)
    if (defined($endstep) && $endstep == 0) {
        push(@ret,$endpt);
        }

	return @ret;
	}

sub f {
	my $self=shift;
	my $x=shift;
	my @res;
	foreach ($self->getSegsInRange([$x,undef])) {
		if (!wantarray) {my $this=$_->f($x);if ($this || $this eq '0') {push(@res,$this);last;} } # call segment function in scalar context, and avoid collecting points after the first one is found
		else {my @these = $_->f($x);push(@res,@these);}
		}
	my %dupsieve;
	@res = grep {!$dupsieve{$_}++} @res;
	return wantarray ? @res : $res[0];
	}
sub F {
	my $self=shift;
	my $y=shift;
	#my $scalefactor = @_?shift:1; #line seems misguided, obsolete
	my @res;
	foreach ($self->getSegsInRange([undef,$y])) {
		my @these;
		if (!wantarray) {$these[0]=$_->F($y);if ($these[0] || $these[0] eq 0) {push(@res,@these);last;} } # call segment function in scalar context, and avoid collecting points after the first one is found
		else {@these = $_->F($y);push(@res,@these);}
		push(@res,@these);
		}
	my %dupsieve;
	@res = grep {!$dupsieve{$_}++} @res;
	return wantarray ? @res : $res[0];
	}
sub point {
	my $self = shift;
	my $theta = shift;
	my @seg;
	if (!defined($self->{pathThetaToCompThetaRanges})) {$self->parameterize();} #delay this until you need it
	@seg=$self->getSegThetaIndexAtPathTheta($theta);
	#print "path->point($theta) find seg?: $seg[0], $seg[1], $seg[2]\n";
	return $seg[0]->point($seg[1]);
	}
sub getSegThetaIndexAtPathTheta {
# how do you map segment's non-uniform theta spacing to your path-wide pseudo-theta indexing?
# here I'm tring to find the segments native theta that correspondes to a given fraction of it's length
# Expensive, but should be doable. Make it work, then decide if it's too expensive.
	my $self = shift;
	my $theta = shift;
	if (!defined($self->{pathThetaToCompThetaRanges})) {$self->parameterize();} #delay this until you need it
	for (my $i=0;$i < scalar(@{$self->{pathThetaToCompThetaRanges}});$i++) {
		#warn "test:  (($self->{pathThetaToCompThetaRanges}->[$i]->[0] < $theta || $self->{pathThetaToCompThetaRanges}->[$i]->[0] eq $theta)  &&  ($self->{pathThetaToCompThetaRanges}->[$i]->[1] > $theta || $self->{pathThetaToCompThetaRanges}->[$i]->[1] eq $theta)) ? ",((($self->{pathThetaToCompThetaRanges}->[$i]->[0] < $theta || $self->{pathThetaToCompThetaRanges}->[$i]->[0] eq $theta) && ($self->{pathThetaToCompThetaRanges}->[$i]->[1] > $theta || $self->{pathThetaToCompThetaRanges}->[$i]->[1] eq $theta))?'yes':'no'),"\n";
		if (
			($self->{pathThetaToCompThetaRanges}->[$i]->[0] < $theta || $self->{pathThetaToCompThetaRanges}->[$i]->[0] eq $theta)
			&&
			($self->{pathThetaToCompThetaRanges}->[$i]->[1] > $theta || $self->{pathThetaToCompThetaRanges}->[$i]->[1] eq $theta)
			) {
			my $segoffsetlengthfraction = ($theta - $self->{pathThetaToCompThetaRanges}->[$i]->[0]) / ($self->{pathThetaToCompThetaRanges}->[$i]->[1] - $self->{pathThetaToCompThetaRanges}->[$i]->[0]);
			my $segtheta=$segoffsetlengthfraction; #might that be good enough? w/o worrying about real length issues?
			return ($self->{pathSegments}->[$i],$segtheta,$i);
			}
		}
	}
sub getPathThetaAtSegTheta {
	my $self = shift;
	my $segindex  = shift;
	my $segtheta = shift;
	if (!defined($self->{pathThetaToCompThetaRanges})) {$self->parameterize();} #delay this until you need it
	my $theta = ( $segtheta * ($self->{pathThetaToCompThetaRanges}->[$segindex]->[1] - $self->{pathThetaToCompThetaRanges}->[$segindex]->[0]) ) + $self->{pathThetaToCompThetaRanges}->[$segindex]->[0];
	#print "getPathThetaAtSegTheta:\n$theta = ( $segtheta * (",$self->{pathThetaToCompThetaRanges}->[$segindex]->[1]," - ",$self->{pathThetaToCompThetaRanges}->[$segindex]->[0],") ) + ",$self->{pathThetaToCompThetaRanges}->[$segindex]->[0],";\n";
	return $theta;
	}
sub curve {
	my $self = shift;
	my $cnt  = @_?shift:20;
	my @ret;
	#my $l=$self->getLength();
	#my $inc=$l/$cnt;
	my $inc=1/$cnt;
	#the sprintf is attmpting to overcome increment overshoot at "1" due to inherent rounding errors
	for (my $i=0;$i<1;$i+=$inc) {push(@ret,$self->point($i));}
	push(@ret,$self->point(1));
	#my $inc = abs($self->{maxx} - $self->{minx})/$cnt;
	#if ($inc>0) {
	#	foreach (my $i=$self->{minx};$i<=$self->{maxx};$i+=$inc) {
	#		#print "f($i) = ",$self->f($i),"\n";
	#		push(@ret,[$i,$self->f($i)]);
	#		if ($i + $inc > $self->{maxx}) {push(@ret,[$self->{maxx},$self->f($self->{maxx})]);last;}
	#		}
	#	}
	#elsif ($inc eq 0) {
	#	$inc = abs($self->{maxy} - $self->{miny})/$cnt;
	#	#print "in MPath curve(), min and max x were same, so using min ($self->{miny}) and max ($self->{maxy}) y, and F()";
	#	foreach (my $i=$self->{miny};$i<=$self->{maxy};$i+=$inc) {
	#		#print "F($i) = ",$self->F($i),"\n";
	#		push(@ret,[$self->F($i),$i]);
	#		if ($i + $inc > $self->{maxy}) {push(@ret,[$self->F($self->{maxy}),$self->{maxy}]);last;}
	#		}
	#	}
	#else {
	#	foreach (my $i=$self->{maxx};$i<=$self->{minx};$i+=$inc) {
	#		#print "f($i) = ",$self->f($i),"\n";
	#		push(@ret,[$i,$self->f($i)]);
	#		if ($i + $inc < $self->{minx}) {push(@ret,[$self->{minx},$self->f($self->{minx})]);}
	#		}
	#	}
	return @ret;
	}
sub secondDerivative {
	my $self=shift;
	my $x=shift;
	my $y=@_?shift:undef;
	my @res;
	foreach ($self->getSegsInRange([$x,$y])) {push(@res,$_->secondDerivative($x,$y));}
	return wantarray ? @res : $res[0];
	}
sub slopeTangent {
	my $self=shift;
	my $x=shift;
	my $y=@_?shift:undef;
	my @res;
	foreach ($self->getSegsInRange([$x,$y])) {push(@res,$_->slopeTangent($x,$y));}
	return wantarray ? @res : $res[0];
	}
sub slopeNormal  {
	my $self=shift;
	my $x=shift;
	my $y=@_?shift:undef;
	my @res;
	foreach ($self->getSegsInRange([$x,$y])) {push(@res,$_->slopeNormal($x,$y));}
	return wantarray ? @res : $res[0];
	}
sub slopeTangent_byTheta {
	my $self=shift;
	my $theta=shift;
	my @seg=$self->getSegThetaIndexAtPathTheta($theta);
	return $seg[0]->slopeTangent_byTheta($seg[1]);
	}
sub slopeNormal_byTheta  {
	my $self=shift;
	my $theta=shift;
	my @seg=$self->getSegThetaIndexAtPathTheta($theta);
	return $seg[0]->slopeNormal_byTheta($seg[1]);
	}
sub angleTangent {
	my $self=shift;
	my $x=shift;
	my $y=@_?shift:undef;
	my @res;
	foreach ($self->getSegsInRange([$x,$y])) {push(@res,$_->angleTangent($x,$y));}
	return wantarray ? @res : $res[0];
	}
sub angleNormal  {
	my $self=shift;
	my $x=shift;
	my $y=@_?shift:undef;
	my @res;
	foreach ($self->getSegsInRange([$x,$y])) {push(@res,$_->angleNormal($x,$y));}
	return wantarray ? @res : $res[0];
	}
sub angleTangent_byTheta {
	my $self=shift;
	my $theta=shift;
	my @seg=$self->getSegThetaIndexAtPathTheta($theta);
	return $seg[0]->angleTangent_byTheta($seg[1]);
	}
sub angleNormal_byTheta  {
	my $self=shift;
	my $theta=shift;
	my @seg=$self->getSegThetaIndexAtPathTheta($theta);
	return $seg[0]->angleNormal_byTheta($seg[1]);
	}

sub solveXforTheta {
	my $self=shift;
	my $x=shift;
	my @thetas;
	my $toff=0;
	for (my $i=0;$i<scalar(@{$self->{pathSegments}});$i++) {
		my @segthetas;
		my ($bx,$by) = $self->{pathSegments}->[$i]->inRange([$x,undef]);
		if ($bx || $by) {
			push(@segthetas,$self->{pathSegments}->[$i]->solveXforTheta($x));
			}
		foreach (@segthetas) {push(@thetas,$toff + (  $self->{pathSegments}->[$i]->getLength(1000,0,$_) / $self->{length}  ));} #uh, I hope
		$toff+=$self->{pathSegments}->[$i]->{length}/$self->{length};
		#print "toff: $toff after i: $i last length: $self->{pathSegments}->[$i]->{length} self length: $self->{length}\n";
		}
	my %dupsieve;
	@thetas = grep {!$dupsieve{$_}++} @thetas;
	return wantarray ? @thetas:$thetas[0];
	}

sub solveYforTheta {
	my $self=shift;
	my $y=shift;
	my @thetas;
	my $toff=0;
	for (my $i=0;$i<scalar(@{$self->{pathSegments}});$i++) {
		my @segthetas;
		my ($bx,$by) = $self->{pathSegments}->[$i]->inRange([undef,$y]);
		if ($bx || $by) {
			push(@segthetas,$self->{pathSegments}->[$i]->solveYforTheta($y));
			}
		foreach (@segthetas) {push(@thetas,$toff + (  $self->{pathSegments}->[$i]->getLength(100,0,$_) / $self->{length}  ));} #uh, I hope
		$toff+=$self->{pathSegments}->[$i]->{length}/$self->{length};
		}
	my %dupsieve;
	@thetas = grep {!$dupsieve{$_}++} @thetas;
	return wantarray ? @thetas:$thetas[0];
	}

sub normalizeY {
	my $self=shift;
	my $nymin = shift;
	my $nymax = shift;
	my $constrain = @_?shift:0;
	my $scalefactor=1;
	my $ydiff=$self->{maxy} - $self->{miny};
	my $oldx=$self->{minx};
	my $untranslatex=$self->{pathSegments}->[0]->{p1}->[0];
	my $untranslatey=$self->{pathSegments}->[0]->{p1}->[1];
	#$self->translate(-$self->{minx},-$self->{miny});
	$self->translate(-$self->{pathSegments}->[0]->{p1}->[0],-$self->{pathSegments}->[0]->{p1}->[1]);
	if (abs($ydiff) > 0.0000000000001) {$scalefactor = ($nymax - $nymin)/$ydiff;}
	my $newspec='';
	#print "yoldspec:",$self->{pathspec},"\n";
	foreach (@{$self->{pathSegmentSpecs}}) {
		my $segTypeLetter = substr($_,0,1);
		my @thesepoints= $self->extractPointsFromPathSpec($_);
		$newspec.=$segTypeLetter;
		my $first=1;
		foreach my $point (@thesepoints) {
			if (ref($point) eq 'ARRAY') {
				#my $newy=$point->[1] * $scalefactor * ($first && $segTypeLetter eq 'A'?(($self->{maxy} - $self->{miny})/(2*$point->[1])):1);
				my $newy=$point->[1] * $scalefactor;
				my $newx=$point->[0];
				if ($constrain && !($first && $segTypeLetter eq 'A')) {
					$newx*=$scalefactor;
					}
# SPEED ATTEMPTS
#				if (ref($newy) && !ref($point->[1])) {$newy=eval(substr($newy->bstr,0,25));}
#				if (ref($newx) && !ref($point->[0])) {$newx=eval(substr($newx->bstr,0,25));}
				if (ref($newy) && !ref($point->[1])) {$newy=0 + sprintf("%.20f",$newy->bstr);}
				if (ref($newx) && !ref($point->[0])) {$newx=0 + sprintf("%.20f",$newx->bstr);}
				$newspec.=$newx.','.$newy.',';
				}
			else {$newspec.=$point.',';}
			$first=0;
			}
		$newspec=substr($newspec,0,-1);

		##if (!$constrain) {$newspec.=$segTypeLetter . join(',',map {(ref($_) eq 'ARRAY')? $_->[0].','.eval(substr(''.($_->[1] * $scalefactor),0,25)) : $_} @thesepoints);}
		#if (!$constrain) {$newspec.=$segTypeLetter . join(',',map {(ref($_) eq 'ARRAY')? $_->[0].','.($_->[1] * $scalefactor) : $_} @thesepoints);}
		##else             {$newspec.=$segTypeLetter . join(',',map {(ref($_) eq 'ARRAY')?eval(substr(($_->[0] * $scalefactor),0,25)).','.eval(substr(($_->[1] * $scalefactor),0,25)) : $_} @thesepoints);}
		#else             {$newspec.=$segTypeLetter . join(',',map {(ref($_) eq 'ARRAY')?($_->[0] * $scalefactor).','.($_->[1] * $scalefactor) : $_} @thesepoints);}
		}
	#print "ynewspec:$newspec\n";
	$self->constructSegments($newspec);
	my $offset = $nymin - $self->{miny};
	#$self->translate($oldx,$offset);
	$self->translate($untranslatex,$offset);
	if ($nymin ne $self->{miny}) {
	#	print "UH OH y: tried to set miny to $nymin, but it came out as $self->{miny}\n";
#		print "try to fix that by doing it with bigfloat (or should you have just rounded ?)\n";
#		my $bignewymin=Math::BigFloat->new(''.$nymin);
#		my $bigymin=Math::BigFloat->new(''.$self->{miny});
#		my $smalleroffset = $bignewymin - $bigymin;
#		if ($smalleroffset) {$self->translate(undef,$smalleroffset);}
#		if ($bignewymin ne $bigymin) {
#			print "     same old problem: $bignewymin ne $bigymin\n"
#			}
		}
	#print "ynewspec:",$self->{pathspec},"\n";
	}
sub normalizeX {
	my $self=shift;
	my $nxmin = shift;
	my $nxmax = shift;
	my $constrain = @_?shift:0;
	my $scalefactor=1;
	#print "norm X xoldspec:",$self->{pathspec},"\n";
	#$self->translate(-$self->{minx},-$self->{miny});
	#lets do this instead: assume the first point of the path is the reference point
	my $untranslatex=$self->{pathSegments}->[0]->{p1}->[0];
	my $untranslatey=$self->{pathSegments}->[0]->{p1}->[1];

	$self->translate(-$self->{pathSegments}->[0]->{p1}->[0],-$self->{pathSegments}->[0]->{p1}->[1]);
	my $xdiff=$self->{maxx} - $self->{minx};
	if (abs($xdiff) > 0.0000000000001) {
		$scalefactor = ($nxmax - $nxmin)/$xdiff;
		}
	my $newspec='';
	#print "xoldspec:$self->{pathspec}\n";
	foreach (@{$self->{pathSegmentSpecs}}) {
		#print "pathsegspec:",$_,"\n";
		my $segTypeLetter = substr($_,0,1);
		my @thesepoints= $self->extractPointsFromPathSpec($_);
		#if (!$constrain) {$newspec.=$segTypeLetter . join(',',map {(ref($_) eq 'ARRAY')?eval(substr(($_->[0] * $scalefactor),0,25)).','.$_->[1]:$_} @thesepoints);}

		$newspec.=$segTypeLetter;
		my $first=1;
		foreach my $point (@thesepoints) {
			if (ref($point) eq 'ARRAY') {
				#my $newx=$point->[0] * $scalefactor * ($first && $segTypeLetter eq 'A'?(($self->{maxx} - $self->{minx})/(2*$point->[0])):1);
				my $newx=$point->[0] * $scalefactor ;
				my $newy=$point->[1];

				#don't recal reason for $first, and now it's causing problem
				#so I'll comment it out, and recreat a problem to be rediscoverd done the road
				#so that maybe I can get what I'm working on today working
				#if ($constrain && !($first && $segTypeLetter eq 'A')) {
				if ($constrain) {
					$newy*=$scalefactor;
					}
# SPEED ATTEMPTS
#				if (ref($newx) && !ref($point->[0])) {$newx=eval(substr($newx->bstr,0,25));}
#				if (ref($newy) && !ref($point->[1])) {$newy=eval(substr($newy->bstr,0,25));}
				if (ref($newx) && !ref($point->[0])) {$newx=0 + sprintf("%.20f",$newx->bstr);}
				if (ref($newy) && !ref($point->[1])) {$newy=0 + sprintf("%.20f",$newy->bstr);}
				$newspec.=$newx.','.$newy.',';
				}
			else {$newspec.=$point.',';}
			$first=0;
			}
		$newspec=substr($newspec,0,-1);#chop off last comma


		#if (!$constrain) {$newspec.=$segTypeLetter . join(',',map {(ref($_) eq 'ARRAY')?($_->[0] * $scalefactor).','.$_->[1]:$_} @thesepoints);}
		##else             {$newspec.=$segTypeLetter . join(',',map {(ref($_) eq 'ARRAY')?eval(substr(($_->[0] * $scalefactor),0,25)).','.eval(substr(($_->[1] * $scalefactor),0,25)):$_} @thesepoints);}
		#else             {$newspec.=$segTypeLetter . join(',',map {(ref($_) eq 'ARRAY')?($_->[0] * $scalefactor).','.($_->[1] * $scalefactor):$_} @thesepoints);}
		}
	#print "xnewspec:$newspec\n";
	$self->constructSegments($newspec);
	my $offset = $nxmin - $self->{minx};

	$self->translate($offset,$untranslatey);

	if ($nxmin ne $self->{minx}) {
	#	print "UH OH x: tried to set minx to $nxmin, but it came out as $self->{minx}\n";
#		print "try to fix that by doing it with bigfloat (or should you have just rounded ?)\n";
#		my $bignxmin=Math::BigFloat->new(''.$nxmin);
#		my $smalleroffset = $bignxmin - $self->{minx};
#		if ($smalleroffset) {$self->translate($smalleroffset,0);}
#		if ($nxmin ne $self->{minx}) {
#			print "     same old problem: $nxmin ne $self->{minx}\n"
#			}
		}
	#print "norm X newspec:",$self->{pathspec},"\n";

	}
sub translate {
	my $self = shift;
	my $xoff = shift;
	my $yoff = shift;
	$xoff ||= 0;
	$yoff ||= 0;
	if (!$xoff && !$yoff) {return;}
	my $xoffbig=$xoff;
	my $yoffbig=$yoff;
# SPEED ATTEMPTS
#	if (ref($xoffbig)) {$xoff = eval(substr($xoffbig->bstr(),0,25));}
#	if (ref($yoffbig)) {$yoff = eval(substr($yoffbig->bstr(),0,25));}
	if (ref($xoffbig)) {$xoff = 0 + sprintf("%.20f",$xoffbig->bstr());}
	if (ref($yoffbig)) {$yoff = 0 + sprintf("%.20f",$yoffbig->bstr());}

	#print "translate:xoff,yoff:$xoff,$yoff\ntranslate:oldspec:$self->{pathspec}\n";
	my $newspec='';
	foreach (@{$self->{pathSegmentSpecs}}) {
		if (substr($_,0,1) eq 'A') {
			my @pts=$self->extractPointsFromPathSpec($_);
			$newspec.=substr($_,0,1) . join(',',map {(ref($_) eq 'ARRAY')?($_->[0]).','.($_->[1]) : $_} @pts[0 .. $#pts - 1]);
# SPEED ATTEMPTS
#			$newspec.=','.eval(sprintf("%.9f",$pts[$#pts]->[0] + $xoff)).','.eval(sprintf("%.9f",$pts[$#pts]->[1] + $yoff));
			$newspec.=','.sprintf("%.9f",$pts[$#pts]->[0] + $xoff).','.sprintf("%.9f",$pts[$#pts]->[1] + $yoff);
			}
		else {
			#$newspec.=substr($_,0,1) . join(',',map {eval(sprintf("%.9f",$_->[0] + $xoff)).','.eval(sprintf("%.9f",$_->[1] + $yoff))} $self->extractPointsFromPathSpec($_));
			#$newspec.=substr($_,0,1) . join(',',map {eval(substr($_->[0] + $xoff,0,25)).','.eval(substr($_->[1] + $yoff,0,25))} $self->extractPointsFromPathSpec($_));
			$newspec.=substr($_,0,1) . join(',',map {($_->[0] + ((ref($_->[0])&&ref($xoff))?$xoffbig:$xoff)).','.($_->[1] + ((ref($_->[1])&&ref($yoff))?$yoffbig:$yoff))} $self->extractPointsFromPathSpec($_));
			}
		}
	$self->constructSegments($newspec);
	#print "translate:newspec:",$self->{pathspec},"\n";
	}
sub scale { # new, untested, copied from javascript version where it seems to work
	my $self = shift;
	my $xscale = shift;
	my $yscale = shift;

    # this stuff was wrong in javascript, but inertly so
    # but now I'm fixing it there and here in a way that might introduce weirdness,
    # since I'm not running and testing this at the moment
	if (!$xscale && $yscale) {$xscale=1;}
	if (!$yscale && $xscale) {$yscale=$xscale;}

	my $osp = "scale:oldspec:\n".$self->{pathspec}."\n";
	my $newspec='';
	for (my $i=0;$i<scalar(@{$self->{pathSegmentSpecs}});$i++) {
		my @pts=$self->extractPointsFromPathSpec($self->{pathSegmentSpecs}->[$i]);
		if (substr($self->{pathSegmentSpecs}->[$i],0,1) eq 'A') {
			$newspec.=substr($self->{pathSegmentSpecs}->[$i],0,1);

            # TODO: copy this expanded, partly fixed arc scale stuff back to js version

            # two radii (these scale, but should remain positive if either scale is negative)
            my $radii = shift @pts;
            # phi, large arc flag, sweep flag
            my $phi=shift @pts;
            my $lrgarc=shift @pts;
            my $sweep=shift @pts;
            # is this right? when one or other scale is negative (but not both)
            # flip sweep flag?
            if (($yscale < 0 && $xscale > 0) || ($yscale < 0 && $xscale > 0)) {
                if ($sweep) {$sweep=0;}
                else        {$sweep=1;}
                }

            # TODO: seems you should take phi into account here
            if ($phi ne 0) {die "not sure you want to simply scale radii when phi in effect - come address this";}

            $newspec.= ($radii->[0]*abs($xscale)) . ',' . ($radii->[1]*abs($yscale)) . ',';
            $newspec.= $phi.',';
            $newspec.= $lrgarc.',';
            $newspec.= $sweep;
			# end of arc point coords
			$newspec.=','.( ($pts[$#pts]->[0]) * $xscale).','.( ($pts[$#pts]->[1]) * $yscale);
			}
		else {
			$newspec.=substr($self->{pathSegmentSpecs}->[$i],0,1) . join(',',map { ($_->[0] * $xscale).','.($_->[1] * $yscale) } @pts);}
			}
	$self->constructSegments($newspec);
	print $osp,"scale:newspec:\n",$newspec,"\n";
	}
sub constructSegment {
	my $self = shift;
	my $segspec = shift;
	my $segspecprevious = shift;
	my $segTypeLetter = substr($segspec,0,1);
	my @lastpoints = $self->extractPointsFromPathSpec($segspecprevious);
	my @thesepoints= $self->extractPointsFromPathSpec($segspec);
	if    ($segTypeLetter eq 'C') {
		if ( # check for degenerate curve
			($lastpoints[$#lastpoints]->[0] eq $thesepoints[0]->[0] && $lastpoints[$#lastpoints]->[0] eq $thesepoints[1]->[0] && $lastpoints[$#lastpoints]->[0] eq $thesepoints[2]->[0]) ||
			($lastpoints[$#lastpoints]->[1] eq $thesepoints[0]->[1] && $lastpoints[$#lastpoints]->[1] eq $thesepoints[1]->[1] && $lastpoints[$#lastpoints]->[1] eq $thesepoints[2]->[1])
			) {
			my $ret=Math::MPath::LineSegment->new($lastpoints[$#lastpoints],$thesepoints[2],$self->{precision},$self->{isLite});
			return $ret;
			}
		else {
			return Math::MPath::BezierCubicSegment->new($lastpoints[$#lastpoints],@thesepoints,$self->{precision},$self->{isLite});
			}
		}
	elsif ($segTypeLetter eq 'L') {return Math::MPath::LineSegment->new($lastpoints[$#lastpoints],@thesepoints,$self->{precision},$self->{isLite});}
	elsif ($segTypeLetter eq 'H') {return Math::MPath::LineSegment->new($lastpoints[$#lastpoints],@thesepoints,$self->{precision},$self->{isLite});}
	elsif ($segTypeLetter eq 'V') {return Math::MPath::LineSegment->new($lastpoints[$#lastpoints],@thesepoints,$self->{precision},$self->{isLite});}
   	elsif ($segTypeLetter eq 'M') {return Math::MPath::MoveTo->new(     $lastpoints[$#lastpoints],@thesepoints,$self->{precision},$self->{isLite});}
	elsif ($segTypeLetter eq 'Z') {return Math::MPath::ClosePath->new(  $lastpoints[$#lastpoints],@thesepoints,$self->{precision},$self->{isLite});}
	elsif ($segTypeLetter eq 'z') {return Math::MPath::ClosePath->new(  $lastpoints[$#lastpoints],@thesepoints,$self->{precision},$self->{isLite});}
	elsif ($segTypeLetter eq 'A') {return Math::MPath::EllipticalArc->new($lastpoints[$#lastpoints],@thesepoints,$self->{precision},$self->{isLite});}

	}
sub extractPointsFromPathSpec {
	my $self=shift;
	my $segspec = shift;
	#print "EXTRACTING FROM:$segspec\n";
	my @ret;
	my $segTypeLetter = substr($segspec,0,1);
	if    ($segTypeLetter eq 'M') {push @ret , [split(/[ ,]+/,substr($segspec,1))]}
	elsif ($segTypeLetter eq 'L') {push @ret , [split(/[ ,]+/,substr($segspec,1))]}
	elsif ($segTypeLetter eq 'H') {push @ret , [split(/[ ,]+/,substr($segspec,1))]}
	elsif ($segTypeLetter eq 'V') {push @ret , [split(/[ ,]+/,substr($segspec,1))]}
	#coords parsed here for Z command were injected upstream, but shouldn't appear in representations of SVG path specs
	elsif ($segTypeLetter eq 'Z') {push @ret , [split(/[ ,]+/,substr($segspec,1))]}
	elsif ($segTypeLetter eq 'z') {push @ret , [split(/[ ,]+/,substr($segspec,1))]}
	elsif ($segTypeLetter eq 'C') {my ($cp1x,$cp1y,$cp2x,$cp2y,$px,$py) = split(/[ ,]+/,substr($segspec,1));push @ret , ([$cp1x,$cp1y],[$cp2x,$cp2y],[$px,$py])}
	elsif ($segTypeLetter eq 'A') {my ($rx,$ry,$phi,$lrgarc,$sweep,$px,$py) = split(/[ ,]+/,substr($segspec,1));push @ret , ([$rx,$ry],$phi,$lrgarc,$sweep,[$px,$py])}
	for (my $i=0;$i<@ret;$i++) {
		if (ref($ret[$i]) eq 'ARRAY') {
			foreach my $num (@{$ret[$i]}) {
				my $sigdigits=0;
				if ($num =~ /^-?([0-9]+\.[0-9]+[1-9])0*$/i) {$sigdigits=length($1) - 1;}
				if ($num =~ /^-?0\.0*([1-9][0-9]+?[1-9])0*$/) {$sigdigits=length($1);} #small numbers that might fit in perl float
				if ($num =~ /^-?([1-9][0-9]+?[1-9])0*\.?0*$/) {$sigdigits=length($1);} #big numbers that might fit in perl float
				}
			}
		}
	return @ret;
	}

sub getFeet {
	my $self = shift;
	my $x= shift;
	my $y= shift;
	my @feet;
	for (my $i=0;$i< scalar(@{$self->{pathSegments}});$i++) {
		my @f=$self->{pathSegments}->[$i]->getFeet($x,$y);
		if (scalar(@f)) {
			for (my $j=0;$j<scalar(@f);$j++) {
				$f[$j]->[2] = $self->getPathThetaAtSegTheta($i,$f[$j]->[2]);
				}
			push(@feet,@f);
			}
		}
	#if (!scalar(@feet)) {
	#	#print "NOOOO FEET for [$x,$y]\n";
	#	}
	return @feet;
	}
sub getIntersections {
	my $self = shift;
	my $other= shift;
	my $wantThetas = scalar(@_) ? shift:0;
	my $wantNativeThetas = scalar(@_) ? shift:0;
	my @intersects=();
	for (my $i=0;$i<@{$self->{pathSegments}};$i++) {
		foreach my $otherseg (@{$other->{pathSegments}}) {
			push(@intersects,map {
			 ($wantThetas && !$wantNativeThetas)
			 ? $self->getPathThetaAtSegTheta($i,$_)
			 : $_
			 } getSegSegIntersects($self->{pathSegments}->[$i],$otherseg,$wantThetas));
			}
		}
	return @intersects;
	}
sub getSegSegIntersects {
	my $seg1=shift;
	my $seg2=shift;
	my $wantThetas = scalar(@_) ? shift:0;
	my $refstrings = ref($seg1).'--'.ref($seg2);
	my @ret;

	if (($refstrings=~/LineSegment/ || $refstrings=~/ClosePath/) && $refstrings=~/BezierCubicSegment/) {
		my $line;
		my $curve;
		my $lineIsSelf;
		if ($refstrings=~/LineSegment.*?--/ || $refstrings=~/ClosePath.*?--/) {$lineIsSelf=1;($line,$curve) = ($seg1,$seg2);}
		else                                 {$lineIsSelf=0;($line,$curve) = ($seg2,$seg1);}

		# t^3 + [(F-mB)/(E-mA)]t^2 + [(G-mC)/(E-mA)]t + [(H-mD-x0)/(E-mA)] = 0 , I hope

        # Yeah, but you should have noted why. How'd you work that up again?
        # You took the equation of a line y = mx + b
        # Used the Bezier's y(t) and x(t) for y and x, used the line's m and b,
        # then arranged it like this: y - mx - b = 0
        # and the left side comes out as a cubic polynomial in t that we can
        # use the good old cubic solver for. Sounds right.

        # First some special cases for vertical and horizontal lines.
		my @thetas;
		if ($line->{m} eq 'inf' || $line->{m} eq '-inf') {
			#then I hope this is right
			@thetas = $curve->solveXforTheta($line->{maxx});
			foreach my $t (@thetas) {
				my $y = $curve->bezierEvalYofT($t);
				#print "an y: $y\n";
				#print "line maxy: $line->{maxy}\n";
				#if ($y <= $line->{maxy} + $line->{precision} &&
				#	$y >= $line->{miny} - $line->{precision}) {
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
			#warn "zero slope, snap special";
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

			#@thetas = &cubicformula(($curve->{F_Big}-$line->{m}*$curve->{B_Big})/($curve->{E_Big}-$line->{m}*$curve->{A_Big}),($curve->{G_Big}-$line->{m}*$curve->{C_Big})/($curve->{E_Big}-$line->{m}*$curve->{A_Big}),($curve->{H_Big}-$line->{m}*$curve->{D_Big}-$line->{b})/($curve->{E_Big}-$line->{m}*$curve->{A_Big}),1);
			#print "uh thetas:",join(",",@thetas),"  and line slope:$line->{m}\n";
			@thetas = sort {$a<=>$b} grep {(1 > $_ || 1 eq $_) && ($_ > 0 || $_ eq 0)}  @thetas;
			#print "uh _Big thetas:",join(",",@thetas),"\n";
			#@thetas = sort {$a<=>$b} grep {(1 > $_ || 1 eq $_) && ($_ > 0 || $_ eq 0)} @thetas;
			foreach my $t (@thetas) {
				my $x = $curve->bezierEvalXofT($t);
				#print "an x: $x\n";
				#print "line maxx: $line->{maxx}\n";
				#if ($x <= $line->{maxx} + $line->{precision} &&
				#	$x >= $line->{minx} - $line->{precision}) {
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
		}
	elsif ($refstrings=~/(LineSegment|ClosePath).*?--/ && $refstrings=~/--.*?(LineSegment|ClosePath)/) {
#		# m1x + b1 = m2x + b2
#		# (m1-m2)x = b2-b1
#		# x = (b2-b1)/(m1-m2)
#		my $x;
#		my $dm = ($seg1->{m}-$seg2->{m});
#		if    (($seg1->{m} eq 'inf' || $seg1->{m} eq '-inf') && $seg2->{m} ne 'inf' && $seg2->{m} ne '-inf') {$x = $seg1->{p1}->[0];}
#		elsif (($seg2->{m} eq 'inf' || $seg2->{m} eq '-inf') && $seg1->{m} ne 'inf' && $seg1->{m} ne '-inf') {$x = $seg2->{p1}->[0];}
#		elsif (abs($dm) > 0.000000000001) {$x=($seg2->{b}-$seg1->{b})/$dm;}
#		#my $mostprec = ($seg1->{precision}<$seg2->{precision})?$seg1->{precision}:$seg2->{precision};
#		if (defined($x) &&
#			#$x <= $seg1->{maxx} + $mostprec &&
#			#$x >= $seg1->{minx} - $mostprec &&
#			#$x <= $seg2->{maxx} + $mostprec &&
#			#$x >= $seg2->{minx} - $mostprec) {
#			($x < $seg1->{maxx} || $x eq $seg1->{maxx} ) &&
#			($x > $seg1->{minx} || $x eq $seg1->{minx} ) &&
#			($x < $seg2->{maxx} || $x eq $seg2->{maxx} ) &&
#			($x > $seg2->{minx} || $x eq $seg2->{minx} ) ) {
#			if ($wantThetas) {push(@ret,($seg1->{m} eq 'inf')?$seg1->solveYforTheta($seg2->f($x)):$seg1->solveXforTheta($x));}
#			else {push(@ret,[$x,($seg1->{m} eq 'inf' || $seg1->{m} eq '-inf')?$seg2->f($x):$seg1->f($x)]);}
#			}



## start new stuff adapted from best working seg_line_intersection() function from elsewhere
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
		#print "M1:$m1 , M2:$m2, DM:$dm, XI: $xi\n";
		my @lowhiu=($u2>$u1)?($u1,$u2):($u2,$u1);
		if ($m1 ne 'Inf') {
			my @lowhix=($x2>$x1)?($x1,$x2):($x2,$x1);
			if ($m2 eq 'Inf' &&   ($u2<$lowhix[0] || $u2>$lowhix[1]) ) {
				#nothing - because no intersection
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
				if ($m2 eq 'Inf' && # in this case we set $xi above even though thre might not be an intersection. If $y not in range of other seg's y extremes, no intersection
					($y<$lowhiv[0] || $y>$lowhiv[1])
					) {
					#nothing
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
### end new stuff


		}
	elsif (($refstrings=~/LineSegment/ || $refstrings=~/ClosePath/) && $refstrings=~/EllipticalArc/) {
		my $line;
		my $arc;
		my $lineIsSelf;
		if ($refstrings=~/LineSegment.*?--/ || $refstrings=~/ClosePath.*?--/) {$lineIsSelf=1;($line,$arc)=($seg1,$seg2);}
		else                                 {$lineIsSelf=0;($line,$arc)=($seg2,$seg1);}
		my @intersections;
		my $x1=$line->{p1}->[0];
		my $y1=$line->{p1}->[1];
		my $x2=$line->{p2}->[0];
		my $y2=$line->{p2}->[1];
		#print "line-arc intercepts:\n    line : [$x1,$y1] , [$x2,$y2]\n";
		$x1-=$arc->{cx};
		$y1-=$arc->{cy};
		$x2-=$arc->{cx};
		$y2-=$arc->{cy};
		#print "if ($line->{maxx}>$arc->{minx} && $line->{minx}<$arc->{maxx} && $line->{maxy}>$arc->{miny} && $line->{miny}<$arc->{maxy}) {\n";
		#if ($line->{maxx}>$arc->{minx} && $line->{minx}<$arc->{maxx} && $line->{maxy}>$arc->{miny} && $line->{miny}<$arc->{maxy}) {
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

				#print "phi : $arc->{phi_radians}\nsin(phi)/cos(phi):",sin($arc->{phi_radians}),"/",cos($arc->{phi_radians}),"\n";
				#print "rot line: [$rot_line_p1->[0],$rot_line_p1->[1]],[$rot_line_p2->[0],$rot_line_p2->[1]]\n";

				my $a = (($rot_line_slope)**2/$arc->{ry}**2) + 1/$arc->{rx}**2;
				my $b = ( 2 * ($rot_line_slope) * ($rot_line_p1->[1] - ($rot_line_slope)*$rot_line_p1->[0]))/$arc->{ry}**2;
				my $c =(($rot_line_p1->[1] - ($rot_line_slope)*$rot_line_p1->[0])**2 / $arc->{ry}**2 ) - 1;
				#print "quad coeffs: $a, $b, $c\n";
				my @xs = &quadraticformula($a,$b,$c,1);
				#print "solution(s) from quad form:",join(",",@xs),"\n";
				for (my $i=0;$i<@xs;$i++) {
					my $y=$rot_line_slope * $xs[$i] + ($rot_line_p1->[1] - $rot_line_slope * $rot_line_p1->[0]); #line formula
					push(@intersections,[$xs[$i],$y]);
					}
				}
			else {#vertical line - use ellipse formula to get points
				#print "VERT (min/max x: $line->{minx},$line->{maxx}) line [$line->{p1}->[0],$line->{p1}->[1] , $line->{p2}->[0],$line->{p2}->[1]] doin ok\n";
				my $y=sqrt($arc->{ry}**2 * (1 - ($x1**2)/($arc->{rx}**2)));#vertical line. use ellipse formula to get the +/- y vals
				push(@intersections,[$x1,$y],[$x1,-$y]);
				#print "   ints even: [ $x1 , +/- $y]\n";
				}
			#for (my $i=0;$i<@intersections;$i++) {
			#	print "centered, unrotated intersection $i: [$intersections[$i]->[0],$intersections[$i]->[1]]\n";
			#	}
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
				#print "add thing: $intersections[$i]->[0]+=$arc->{cx} = ";
				$intersections[$i]->[0]+=$arc->{cx};
				#print "$intersections[$i]->[0]  then  ";
				my $sigcnt3=length(($intersections[$i]->[0]=~/\.([0-9]*)/g)[0]);
				$intersections[$i]->[0]=0 + sprintf("%.".(sort {$b<=>$a} ($sigcnt1,$sigcnt2))[0]."f",$intersections[$i]->[0]);
				#print "$intersections[$i]->[0]\n";
				$intersections[$i]->[1]+=$arc->{cy};
				}

			#Now check to see of those intersections are within bounds and within sweep

#if ( 0 &&   #degug thing
#	! grep {
#					($_->[0] < $line->{maxx} || abs($_->[0] - $line->{maxx}) < $line->{precision}) &&
#					($_->[0] > $line->{minx} || abs($_->[0] - $line->{minx}) < $line->{precision}) &&
#					($_->[1] < $line->{maxy} || abs($_->[1] - $line->{maxy}) < $line->{precision}) &&
#					($_->[1] > $line->{miny} || abs($_->[1] - $line->{miny}) < $line->{precision})
#
#				} @intersections
#	) {
#	map {
#	print "			($_->[0] < $line->{maxx} || abs($_->[0] - $line->{maxx}) < $line->{precision}) && \n";
#	print "			($_->[0] > $line->{minx} || abs($_->[0] - $line->{minx}) < $line->{precision}) && \n";
#	print "			($_->[1] < $line->{maxy} || abs($_->[1] - $line->{maxy}) < $line->{precision}) && \n";
#	print "			($_->[1] > $line->{miny} || abs($_->[1] - $line->{miny}) < $line->{precision})\n";
#	 } @intersections;
#
#	}
				@intersections = grep {
				#$_->[0] <= $line->{maxx} + $line->{precision} &&
				#$_->[0] >= $line->{minx} - $line->{precision} &&
				#$_->[1] <= $line->{maxy} + $line->{precision} &&
				#$_->[1] >= $line->{miny} - $line->{precision}

				#($_->[0] < $line->{maxx} || $_->[0] eq $line->{maxx}) &&
				#($_->[0] > $line->{minx} || $_->[0] eq $line->{minx}) &&
				#($_->[1] < $line->{maxy} || $_->[1] eq $line->{maxy}) &&
				#($_->[1] > $line->{miny} || $_->[1] eq $line->{miny})

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
                        #warn "CIR--LINE allarcthetas: ",join(',',@allArcThetas),"\n";
						foreach my $t (@allArcThetas) {
							my $tp=$arc->point($t);
							if (abs($tp->[1] - $int->[1]) < 0.0000000001) {push(@ret,$t);}
							#man.
							}
						}
					}
				}
			else {
				push(@ret,@intersections);
				}
			#print "got ",scalar(@intersections)," line-arc intersections\n";
			#print "  it/they is/are: \n  ",join("\n  ",map {"[$_->[0],$_->[1]]"} @intersections),"\n";

			}
		}
    # circle-circle special (but common) case
	elsif (   $refstrings=~/EllipticalArc.*?--/ && $refstrings=~/--.*?EllipticalArc/
	       && $seg1->{rx} eq $seg1->{ry}
	       && $seg2->{rx} eq $seg2->{ry}
	      ) {

		my $arc1 = $seg1;
		my $arc2 = $seg2;

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

        #warn "  circ--circ pre filter\n  ",join("\n  ",map {"[$_->[0],$_->[1]]"} @intersections),"\n";

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

            #warn "  circ--circ post arc1 filter\n  ",join("\n  ",map {"[$_->[0],$_->[1]]"} @intersections),"\n";

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


			#warn "got ",scalar(@intersections)," circle_arc-circle_arc intersections\n";
			#warn "  it/they is/are: \n  ",join("\n  ",map {"[$_->[0],$_->[1]]"} @intersections),"\n";

			if ($wantThetas) {
				foreach my $int (@intersections) {
					my @allArcThetas=$arc1->solveXforTheta($int->[0]);
                    #warn "CIR--CIR allarcthetas: ",join(',',@allArcThetas),"\n";
					foreach my $t (@allArcThetas) {
						my $tp=$arc1->point($t);
						if (abs($tp->[1] - $int->[1]) < 0.0000000001) {push(@ret,$t);}
						#man. # why man? were you complaining about doing a tolerance thing here? whatever man!
						}
					}
				}
			else {
				push(@ret,@intersections);
				}
			#print "got ",scalar(@intersections)," circle_arc-circle_arc intersections\n";
			#print "  it/they is/are: \n  ",join("\n  ",map {"[$_->[0],$_->[1]]"} @intersections),"\n";
            }

        }

    # general ellipse-ellipse case
	elsif (   $refstrings=~/EllipticalArc.*?--/ && $refstrings=~/--.*?EllipticalArc/ ) {
        # will be kinda hairy. probably needs quartic solver, unless
        # you can find shortcut due to working with arcs and not full ellipses, in general
        # or unless you cook up a quick (to code) rootfinding appoach
        die "elliptical arc--elliptical arc intersection not handled yet (when both aren't circular arcs)";
        }

    # yeah bez--bez isn't here! you thought you had it but you don't
    # but now you might be able to figure it, based on breaking bezs down to y(t(x)) sub functions
    # -- you'd solve intersection of those. Gotta work all that up though to see if that works.
    # You'll also be in position to do offset bez and offset intersections then, I think.
    # But need serious work sessions for that. Better space. More peace, for longer.

	return @ret;
	}

# this bigsqrt() was the result of a bug killing mission related to dimensionalStepFromTheta()
# useful. Just hiding behind greater precision, but gets the job done.
our $BigFloatOneHalf = Math::BigFloat->new('0.5');
our $BigFloatTen     = Math::BigFloat->new('10');
sub bigsqrt {
	#because the sqrt and root functions in Math::BigFloat sometimes fail
	#Wikipedia reminds us that:
	#sqrt(x) = 10**(1/2 * log_10(x))
    # doing similar thing in QuadraticFormula.pm. see comments there
	return $BigFloatTen->copy()->bpow($BigFloatOneHalf->copy()->bmul($_[0]->copy()->blog(10)),25);
	}

sub dimensionalStepFromTheta {
	my $self=shift;
	#print " in dimstep ref:",ref($self),"\n";

	my $dim=shift;
	my $theta=shift;
	my $direction=scalar(@_)?shift:1; # 1 or 0

	my $findnexttheta = sub {
		my $ret;
		#print " in sub dimstep ref:",ref($self),"\n";
		my $pt_last = $self->point($theta);
		if (!ref($_[0])) {
			my $pt_new  = $self->point($_[0]);
			$ret = $dim - CORE::sqrt(($pt_new->[0] - $pt_last->[0])**2 + ($pt_new->[1] - $pt_last->[1])**2);
			#warn "$ret = $dim - CORE::sqrt(($pt_new->[0] - $pt_last->[0])**2 + ($pt_new->[1] - $pt_last->[1])**2) = $ret\n";
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
	($newtheta,$er) = FalsePosition($findnexttheta,($direction ? [$theta,1]:[0,$theta]),$self->{resolution}/10,($direction ? ($theta + (1-$theta)/2):($theta/2)),'dimensionalStepFromTheta');
	#warn " dim step result ($newtheta,$er)\n";
	if (defined($er)) {
        warn "dimstep er: $er";
		#probably just reached the end
		if (abs(&{$findnexttheta}(($direction ? 1:0))) < $dim) {
			$newtheta=($direction ? 1:0);
			}
		#otherwise the error might be real
		}
	return ($self->point($newtheta),$newtheta);
	}
sub getInfiniteSlopeThetas {
	my $self=shift;
	my $bounds=shift;
	if (!defined($bounds)) {$bounds=[0,1];}
	my @infthetas;
	my ($seg,$segtheta,$segind);
	($seg,$segtheta,$segind)=$self->getSegThetaIndexAtPathTheta($bounds->[0]);
	for (my $i=$segind;$i<scalar(@{$self->{pathSegments}});$i++) {
		if ($self->{pathThetaToCompThetaRanges}->[$segind]->[1] > $bounds->[1]) {last;}
		push(@infthetas,$self->{pathSegments}->[$segind]->getInfiniteSlopeThetas());
		}
	@infthetas = grep {($_ > $bounds->[0] || $_ eq $bounds->[0]) && ($_ < $bounds->[1] || $_ eq $bounds->[1])} @infthetas;
	return @infthetas;
	}

sub getKeyPointsOnLine { ### Not quite right yet, but maybe usable
	my $self=shift;
	my @ret=();
	push(@ret,
		[$self->{pathSegments}->[0]->{maxx},$self->{pathSegments}->[0]->f($self->{pathSegments}->[0]->{maxx}),ref($self->{pathSegments}->[0])=~/bezier/i?$self->{pathSegments}->[0]->solveXforTheta($self->{pathSegments}->[0]->{maxx}):$self->{pathSegments}->[0]->{maxx}],
		[$self->{pathSegments}->[0]->{minx},$self->{pathSegments}->[0]->f($self->{pathSegments}->[0]->{minx}),ref($self->{pathSegments}->[0])=~/bezier/i?$self->{pathSegments}->[0]->solveXforTheta($self->{pathSegments}->[0]->{minx}):$self->{pathSegments}->[0]->{minx}],
		(map {[$_,$self->{pathSegments}->[0]->{maxy},ref($self->{pathSegments}->[0])=~/bezier/i?$self->{pathSegments}->[0]->solveYforTheta($_):$_] }   $self->{pathSegments}->[0]->F($self->{pathSegments}->[0]->{maxy})),
		(map {[$_,$self->{pathSegments}->[0]->{miny},ref($self->{pathSegments}->[0])=~/bezier/i?$self->{pathSegments}->[0]->solveYforTheta($_):$_] }   $self->{pathSegments}->[0]->F($self->{pathSegments}->[0]->{miny})),
		[$self->{pathSegments}->[0]->{p1}->[0],$self->{pathSegments}->[0]->{p1}->[1],ref($self->{pathSegments}->[0])=~/bezier/i?$self->{pathSegments}->[0]->solveXforTheta($self->{pathSegments}->[0]->{p1}->[0]):$self->{pathSegments}->[0]->{p1}->[0]]
		);
	if (ref($self->{pathSegments}->[0])=~/bezier/i) {print "sorting bezier points by thetas1\n";@ret = sort {$a->[2]<=>$b->[2]} @ret;}
	elsif (ref($self->{pathSegments}->[0])=~/line/i) {
		if ($self->{pathSegments}->[0]->[0]<$self->{pathSegments}->[0]->{p2}->[0]) {@ret = sort {$a->[0]<=>$b->[0]} @ret;}
		else {@ret = sort {$b->[0]<=>$a->[0]} @ret;}
		}
#	print "first sort numbers : " , join(",",map {$_->[2]} @ret),"\n";
	@ret = map {[$_->[0],$_->[1]]} @ret;

	push(@ret,[111,111]);#testing

	for (my $i=1;$i<scalar(@{$self->{pathSegments}});$i++) {
		# so far, assumes all segments have ->{p1} and ->{p2} start and end points
		my @thisgroup=(
			[$self->{pathSegments}->[$i]->{maxx},$self->{pathSegments}->[$i]->f($self->{pathSegments}->[$i]->{maxx}), ref($self->{pathSegments}->[$i])=~/bezier/i?$self->{pathSegments}->[$i]->solveXforTheta($self->{pathSegments}->[$i]->{maxx}):$self->{pathSegments}->[$i]->{maxx}],
			[$self->{pathSegments}->[$i]->{minx},$self->{pathSegments}->[$i]->f($self->{pathSegments}->[$i]->{minx}), ref($self->{pathSegments}->[$i])=~/bezier/i?$self->{pathSegments}->[$i]->solveXforTheta($self->{pathSegments}->[$i]->{minx}):$self->{pathSegments}->[$i]->{minx}],
			(map {[$_,$self->{pathSegments}->[$i]->{maxy}, ref($self->{pathSegments}->[$i])=~/bezier/i?$self->{pathSegments}->[$i]->solveYforTheta($_):$_]}   $self->{pathSegments}->[$i]->F($self->{pathSegments}->[$i]->{maxy})),
			(map {[$_,$self->{pathSegments}->[$i]->{miny}, ref($self->{pathSegments}->[$i])=~/bezier/i?$self->{pathSegments}->[$i]->solveYforTheta($_):$_]}   $self->{pathSegments}->[$i]->F($self->{pathSegments}->[$i]->{miny})),
			[$self->{pathSegments}->[$i-1]->{p2}->[0],$self->{pathSegments}->[$i-1]->{p2}->[1], ref($self->{pathSegments}->[$i-1])=~/bezier/i?$self->{pathSegments}->[$i-1]->solveXforTheta($self->{pathSegments}->[$i-1]->{p2}->[0]):$self->{pathSegments}->[$i-1]->{p2}->[0]],
			[$self->{pathSegments}->[$i]->{p1}->[0],  $self->{pathSegments}->[$i]->{p1}->[1],   ref($self->{pathSegments}->[$i])=~/bezier/i?$self->{pathSegments}->[$i]->solveXforTheta($self->{pathSegments}->[$i]->{p2}->[0]):$self->{pathSegments}->[$i]->{p2}->[0]]);
		my %dupsieve={};
		@thisgroup = grep {!$dupsieve{$_->[0].','.$_->[1]}++} @thisgroup;
		if (ref($self->{pathSegments}->[$i])=~/bezier/i) {
#			print "sorting bezier points by thetas2\n";
			@thisgroup = sort {$a->[2]<=>$b->[2]} @thisgroup;
#			print "in loop sort numbers $i: " , join(",",map {$_->[2]} @thisgroup),"\n";
			}
		elsif (ref($self->{pathSegments}->[$i])=~/line/i) {
			if ($self->{pathSegments}->[$i]->[0]<$self->{pathSegments}->[$i]->{p2}->[0]) {@thisgroup = sort {$a->[0]<=>$b->[0]} @thisgroup;}
			else {@thisgroup = sort {$b->[0]<=>$a->[0]} @thisgroup;}
			}
		push(@ret,map {[$_->[0],$_->[1]]} @thisgroup);
		push(@ret,[222,222]);#testing
		}
	push(@ret,$self->{pathSegments}->[scalar(@{$self->{pathSegments}})-1]->{p2});
	my %dupsieve={};
	@ret = grep {!$dupsieve{$_->[0].','.$_->[1]}++} @ret;
	return @ret;
	}


sub _rotate2d {
	my ($origin,$point,$angle) = @_;
	my $dx=($point->[0]-$origin->[0]);
	my $dy=($point->[1]-$origin->[1]);
	#{a c-b d, a d+b c}
	return [$origin->[0] + ($dx*cos($angle) - $dy*sin($angle)),$origin->[1] + ($dx*sin($angle) + $dy*cos($angle))];
	}

######################### end copy-paste old mikepath guts


} # end package Math::MPath


__END__
=pod

=head1 NAME

Math::MPath - An SVG-like Path, with enhancements

=head1 VERSION

version 0.01

=head1 SYNOPSIS

    use Math::MPath qw(); etc.

=head1 ABSTRACT

This module...

=head1 METHODS

=head2 method_one

Does stuff

=for Pod::Coverage new

=head1 ACKNOWLEDGEMENTS

Thanks to ...

=head1 AUTHOR

somebdy <bdy@some.net>

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 4002 by So And So.

What, if any, licence though?

=cut

