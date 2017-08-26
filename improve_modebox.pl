#!/usr/bin/env perl
#

use strict;
use warnings;

use constant PI    => 4 * atan2(1, 1);
use POSIX;

use File::Slurp;
use File::Copy;
use List::Util qw(sum);

use Math::GSL::Matrix;
use Math::Trig qw(acos_real);

my $magic_matrix = compute_magic_matrix();
my $magic_matrix_inverse = $magic_matrix->transpose();


for (@ARGV) {
    improve_modebox($_);
}

sub improve_modebox {
    my $file = shift;
    my $summary_file = "${file}_Summary.txt";
    my $summary_file2 = "${file}_Summary2.txt";
    my $fh;
    my @names = qw(id bx by bz mode obs_cluster obs_-cluster
		   box_cluster box_-cluster
		   density_mode sum_dens_cluster sum_dens_-cluster
		   vol_cluster vol_-cluster
		   mean_dist_mode max_dist_mode
		   peak_cluster
		   mean_dist_all_mode max_dist_all_mode
		   peak_all
		 );
    my @lines = read_file($summary_file, {chomp => 1});
    if (@lines != 1) {
	# We only attemt to do something when a single cluster has been found.
	copy($summary_file, $summary_file2);
	return;
    }
    my %fields;
    @fields{@names} = split /\t/, $lines[0];

    my @obs = read_observations($file);
    my ($bx, $by, $bz, $best) = find_best(\%fields, \@obs);
    if ($bx == $fields{bx} && $by == $fields{by} && $bz == $fields{bz}) {
	# No change, no need to do anything.
	copy($summary_file, $summary_file2);
	return;
    }

    # {
    # 	my $actual = compute_distance(\@obs, map {idx_to_coord($_)} @fields{qw|bx by bz|})/@obs;
    # 	printf STDERR "(%d,%d,%d,%f,%f) => (%d,%d,%d,%f) %s\n", @fields{qw|bx by bz|}, $fields{mean_dist_all_mode},
    # 	    $actual, $bx, $by, $bz, $best, $actual < $best ? "!" : "";
    # }

    # Update the bx,by,bz,mean_dist_all_mode fields and clear out the
    # fields which are certainly meaningless when we modify the 'mode'
    # box.
    @fields{qw|bx by bz|} = ($bx, $by, $bz);
    $fields{mean_dist_all_mode} = $best;

    $fields{mode} = 0;
    $fields{density_mode} = 0;
    $fields{mean_dist_mode} = 0;
    $fields{max_dist_mode} = 0;
    $fields{max_dist_all_mode} = 0;
    write_file($summary_file2, join("\t", @fields{@names}) . "\n");
}

sub read_observations {
    my $file = shift;
    open(my $fh, '<', $file)
	or die "cannot open $file: $!";
    my @obs;
    while (<$fh>) {
	chomp;
	my ($x, $y, $z, $angle) = split;
	push @obs, do_magic_rotation($x, $y, $z, $angle);
    }
    close($fh);
    return @obs;
}

sub find_best {
    my $fields = shift;
    my $obs = shift;
    my @b = @{$fields}{qw|bx by bz|};
    my $dist = compute_squared_distance($obs, map {idx_to_coord($_)} @b);
    my %visited;

    while (1) {
	my @best = @b;
	my $best_dist = $dist;
	my $found_better = 0;
	my @current;
	# printf "%d,%d,%d,%f\n", @best, $best_dist/@$obs;
	for my $dx (-2 ... 2) {
	    $current[0] = $b[0] + $dx;
	    for my $dy (-2 ... 2) {
		$current[1] = $b[1] + $dy;
		for my $dz (-2 ... 2) {
		    $current[2] = $b[2] + $dz;
		    my $s = "$current[0],$current[1],$current[2]";
		    next if exists $visited{$s};
		    my @coord = map {idx_to_coord($_)} @current;
		    # If we've gone outside the sphere of radius pi,
		    # ignore this 'box'.
		    if (sqrt(sum(map {$_*$_} @coord)) >= PI) {
			# print STDERR "skipping @current...\n";
			$visited{$s} = 1000000;
			next;
		    }
		    my $d = compute_squared_distance($obs, @coord);

		    if ($d < $best_dist) {
			@best = @current;
			$best_dist = $d;
			$found_better = 1;
			# printf "%d,%d,%d,%f\n", @best, $best_dist/@$obs;
		    }
		    $visited{$s} = $d;
		}
	    }
	}
	last if (!$found_better);
	@b = @best;
	$dist = $best_dist;
    }
    return @b, compute_distance($obs, map {idx_to_coord($_)} @b)/scalar @$obs;
}

sub compute_distance {
    my ($obs, $x, $y, $z) = @_;
    my $dist = 0;
    my $mat = EV_to_SO3($x, $y, $z);
    for (@$obs) {
	$dist += SO3_distance($mat, $_);
    }
    return $dist;
}
sub compute_squared_distance {
    my ($obs, $x, $y, $z) = @_;
    my $dist = 0;
    my $mat = EV_to_SO3($x, $y, $z);
    for (@$obs) {
	$dist += SO3_distance($mat, $_)**2;
    }
    return $dist;
}

sub coord_to_idx {
    # The boxes are numbered 1..81, with box i covering [-pi +
    # 2pi*(i-1)/81, -pi + 2pi*i/81]. Hence the index corresponding to
    # a given coordinate w is ceil((w+pi)*81/2pi).
    my $w = shift;
    # The int() is to force perl to treat the scalar as an integer; we
    # need its stringification to be e.g. "39", not "39.0" or some
    # other floating point representation. Yes, it's fragile as
    # hell...
    return int(ceil( ($w + PI)*81 / (2*PI) ));
}

sub idx_to_coord {
    # We want the "middle" of the box. According to the above, this
    # means we should return this:
    my $w = shift;
    return -(PI) + (2*PI*($w - 0.5))/81;
}


sub Math::GSL::Matrix::trace {
    my ($mat) = @_;
    my $t = 0;
    return sum(map {$mat->get_elem($_,$_)} (0..2));
}

sub SO3_distance {
    my ($A, $B) = @_;
    my $t = ($A->transpose() * $B)->trace();
    # $t = 3 if ($t > 3);
    # $t = -1 if ($t < -1);
    return acos_real(($t - 1)/2);
}

sub EV_to_SO3 {
    my ($x, $y, $z, $angle);

    if (@_ == 4) {
	($x, $y, $z, $angle) = @_;
	my $n = sqrt($x*$x + $y*$y + $z*$z);
	if ($n) {
	    $x /= $n;
	    $y /= $n;
	    $z /= $n;
	}
    } elsif (@_ == 3) {
	($x, $y, $z) = @_;
	$angle = sqrt($x*$x + $y*$y + $z*$z);
	if ($angle > 0) {
	    $x /= $angle;
	    $y /= $angle;
	    $z /= $angle;
	} else {
	    $x = 1;
	    $y = 0;
	    $z = 0;
	}
    }

    my $K = Math::GSL::Matrix->new(3, 3);
    $K->set_row(0, [0, -$z, $y]);
    $K->set_row(1, [$z, 0, -$x]);
    $K->set_row(2, [-$y, $x, 0]);
    my $K2 = $K*$K;
    return Math::GSL::Matrix->new(3, 3)->identity() + sin($angle)*$K + (1-cos($angle))*$K2;
}

sub SO3_to_EV {
    my $M = shift;
    my $t = $M->trace();
    # $t = 3 if ($t > 3);
    # $t = -1 if ($t < -1);
    my $angle = acos_real(($t - 1)/2);
    # Kludge #1:
    return (0, 0, 0) if $angle == 0;

    if (sin($angle) > 0.0001) {
	# The normal case: Angle is away from 0 and pi.
	my $x = $M->get_elem(2, 1) - $M->get_elem(1, 2);
	my $y = $M->get_elem(0, 2) - $M->get_elem(2, 0);
	my $z = $M->get_elem(1, 0) - $M->get_elem(0, 1);
	my $n = sqrt($x*$x + $y*$y + $z*$z);
	if (abs($n-2*sin($angle)) > 0.001) {
	    printf STDERR "Unexpected large diff (%f) between |v| and 2sin(angle)\n", abs($n-2*sin($angle));
	}
	# $x /= $n;
	# $y /= $n;
	# $z /= $n;
	# return ($x, $y, $z, $angle);
	return ($x, $y, $z);
    }

    if ($t < 0) {
	# Angle is close to pi.
	my $B = ($M + Math::GSL::Matrix->new(3,3)->identity()) * .5;
	my ($x, $y, $z) = (sqrt($B->get_elem(0,0)),
			   sqrt($B->get_elem(1,1)),
			   sqrt($B->get_elem(2,2)));
	if ($x > 0) {
	    $y = -$y if ($M->get_elem(0, 1) < 0);
	    $z = -$z if ($M->get_elem(0, 2) < 0);
	} else {
	    $z = -$z if ($M->get_elem(1, 2) < 0);
	}
	# my $n = sqrt($x*$x + $y*$y + $z*$z);
	# $x /= $n;
	# $y /= $n;
	# $z /= $n;
	# return ($x, $y, $z, $angle);
	return ($x, $y, $z);
    }

}

sub do_magic_rotation {
    my ($x, $y, $z, $angle) = @_;
    $angle //= 1.0;
    # Normalize $x,$y,$z
    my $n = sqrt($x*$x + $y*$y + $z*$z);
    $x /= $n;
    $y /= $n;
    $z /= $n;
    $angle *= $n;
    my $mat = EV_to_SO3($x, $y, $z, $angle);
    # $mat is a Math::GSL::Matrix instance. Can multiply directly by $magic_matrix.
    $mat = $magic_matrix * $mat;

    return $mat;

    ## The caller may also need the matrix, so let him apply
    ## SO3_to_EV, to avoid having to apply EV_to_SO3 to the returned
    ## Euler vector...

    # return SO3_to_EV($mat);
}

sub undo_magic_rotation {
    my ($x, $y, $z, $angle) = @_;
    $angle //= 1.0;
    # Normalize $x,$y,$z
    my $n = sqrt($x*$x + $y*$y + $z*$z);
    if ($n) {
	$x /= $n;
	$y /= $n;
	$z /= $n;
    }
    $angle *= $n;
    my $mat = EV_to_SO3($x, $y, $z, $angle);
    # $mat is a Math::GSL::Matrix instance. Can multiply directly by $magic_matrix.
    $mat = $magic_matrix_inverse * $mat;

    return $mat;

    ## The caller may also need the matrix, so let him apply
    ## SO3_to_EV, to avoid having to apply EV_to_SO3 to the returned
    ## Euler vector...

    # return SO3_to_EV($mat);
}

sub compute_magic_matrix {
    my $delta = 2*PI/81;
    my $xt = 32;
    my $yt = 70;
    my $zt = 31;
    $xt = -(PI) + ($xt - .5)*$delta;
    $yt = -(PI) + ($yt - .5)*$delta;
    $zt = -(PI) + ($zt - .5)*$delta;
    return EV_to_SO3($xt, $yt, $zt)->transpose(); # The inverse is the same as the transpose...
}
