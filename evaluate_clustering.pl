#!/usr/bin/env perl
#
# File: evaluate_clustering.pl
# Time-stamp: <2017-05-01 13:54:43 yuki>
# Author: Rasmus Villemoes

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use List::Util qw(sum);
use File::Slurp;

use constant PI    => 4 * atan2(1, 1);
use POSIX;

use Math::GSL::Matrix;
use Math::Trig qw(acos_real);

use List::Util qw(sum);

use v5.10.1;

my $magic_matrix = compute_magic_matrix();
my $magic_matrix_inverse = $magic_matrix->transpose();


# We accept a pattern as good if:
#
# (1) It resulted in just one cluster
# (2) The peakness of that cluster is at most $max_peak
# (3) The ratio #(minus cluster)/#(plus cluster) is at most $max_minus_plus_ratio
# (4) At least $min_obs points were clustered
#
# If less than $min_obs points were clustered, we give up; that is,
# note that when we reach this pattern, we cannot hope to determine
# the rotation from pattern information alone.
#

my $output_fh = \*STDOUT;

my $max_minus_plus_ratio = 0.1;
my $max_peak = 0.6;
my $min_obs = 50;
my $super_majority = 0.95;
my $too_few = undef;
my $L_max_minus_plus_ratio = undef;
my $L_max_peak = undef;
my $help = 0;
my $directory;

my $flp_option_file;
my $flp_out_dir;

my $output_dir = ".";

my %total;
my %count;

sub print_usage {
    print "$0 [--mp-ratio=X] [--max-peak=Y] [--min-obs=N] <dir with _Summary files>";
    exit(0);
}

sub get_options_from_file {
    my (undef, $file) = @_;
    open(my $fh, '<', $file)
	or die "cannot open $file: $!";
    while (<$fh>) {
	chomp;
	next unless m/^# @ ([a-zA-Z0-9_-]+)\s+(.*)$/;
	my $opt = $1;
	my $val = $2;
	given ($opt) {
	    when ('max-mp-ratio') { $max_minus_plus_ratio = $val; }
	    when ('max-peak') {$max_peak = $val; }
	    when ('L-max-mp-ratio') { $L_max_minus_plus_ratio = $val; }
	    when ('L-max-peak') { $L_max_peak = $val; }
	    when ('min-obs') { $min_obs = $val; }
	    when ('majority') { $super_majority = $val; }
	    when ('too-few') { $too_few = $val; }
	}
    }
    close($fh);
}

GetOptions("max-mp-ratio=f" => \$max_minus_plus_ratio,
           "max-peak=f" => \$max_peak,
	   "L-max-mp-ratio=f" => \$L_max_minus_plus_ratio,
           "L-max-peak=f" => \$L_max_peak,
           "min-obs=i"  => \$min_obs,
	   "too-few=i" => \$too_few,
	   "majority=f" => \$super_majority,
	   "option-file=s" => \&get_options_from_file,
	   "flp-option-file=s" => \$flp_option_file,
	   "flp-out-dir=s" => \$flp_out_dir,
	   "output-dir=s" => \$output_dir,
           "help|h"     => \$help)
    or print_usage();

if ($help) {
    print_usage();
}


# If --too-few not given, default to value of $min_obs
$too_few //= $min_obs;
# For "long-range" Hbonds, we allow a different set of parameters -
# for example, we don't and can't expect these to be quite as 'peaked'
# as shorter bonds. But if not given, we default to the ordinary
# parameters.
$L_max_minus_plus_ratio //= $max_minus_plus_ratio;
$L_max_peak //= $max_peak;

# if (!defined $flp_option_file) {
#     die "either --flp-option-file or --flp-out-dir must be given"
# 	unless defined $flp_out_dir;
#     $flp_option_file = $flp_out_dir . '/' . basename($directory) . '_opts';
# }

if (!defined $flp_out_dir) {
    die "--flp-out-dir must be given";
}

die "only one file allowed" if (@ARGV > 1);
my $inputfile = $ARGV[0];

process_input($inputfile);

sub process_input {
    my $file = shift;
    my @lines;
    my $current_step = '';
    my %processed_steps;

    open(my $ifh, '<', $file)
	or die "cannot open $file: $!";

    while (<$ifh>) {
	my ($summary, $rest) = split /\t/, $_;
	my ($step, $fn) = split /\//, $summary;
	if ($step eq $current_step) {
	    push @lines, $_;
	    next;
	}
	if (@lines) {
	    die "lines belonging to step $current_step not contiguous" if exists $processed_steps{$current_step};
	    process_one_step($current_step, \@lines);
	    $processed_steps{$current_step} = 1;
	    @lines = ();
	}
	$current_step = $step;
	push @lines, $_;
    }
    close($ifh);
    if (@lines) {
	die "lines belonging to step $current_step not contiguous" if exists $processed_steps{$current_step};
	process_one_step($current_step, \@lines);
	$processed_steps{$current_step} = 1;
    }

}

sub process_one_step {
    my $step = shift;
    my $lines = shift;
    my %per_file;
    foreach (@$lines) {
	my ($fn, $rest) = split /\t/, $_, 2;
	push @{$per_file{$fn}}, $rest;
    }
    my @summaries;
    for my $k (keys %per_file) {
	push @summaries, create_summary($k, $per_file{$k});
    }

    $flp_option_file = $flp_out_dir . '/' . $step . '_opts';
    open(my $outfh, '>', $output_dir . '/' . $step . '_assess');
    my $oldfh = $output_fh;
    $output_fh = $outfh;

    v("step %s", $step);
    c("Using");
    v("%-12s %f", "max-mp-ratio", $max_minus_plus_ratio);
    v("%-12s %f", "max-peak", $max_peak);
    v("%-12s %f", "L-max-mp-ratio", $L_max_minus_plus_ratio);
    v("%-12s %f", "L-max-peak", $L_max_peak);
    v("%-12s %d", "min-obs", $min_obs);
    v("%-12s %f", "majority", $super_majority);
    v("%-12s %d", "too-few", $too_few);
    c("");
    v("%-12s %s", "pattern-options",
      join(" ", grep {! m/^#/ } read_file($flp_option_file, {chomp => 1})));
    c("");

    foreach (sort {$b->{total_obs} <=> $a->{total_obs}} @summaries) {
	process_summary($_);
    }

    $output_fh = $oldfh;
    close($outfh);
}


# die "only one directory allowed" if (@ARGV > 1);
# $directory = $ARGV[0];

# Let's try to create output which is both machine- and
# human-readable. The simplest is to just precede all free-form text
# by a #. A third kind of line contains "variable settings", so that
# other scripts reading this file can use the same parameters.

sub c {
    my $fmt = shift;
    my $s = sprintf($fmt, @_);
    $s =~ s/^/# /mg;
    $s .= "\n" unless $s =~ m/\n$/;
    print $output_fh $s;
}

sub v {
    my $fmt = shift;
    my $s = sprintf($fmt, @_);
    $s =~ s/^/# @ /mg;
    $s .= "\n" unless $s =~ m/\n$/;
    print $output_fh $s;
}


sub p {
    my $fmt = shift;
    my $s = sprintf($fmt, @_);
    $s .= "\n" unless $s =~ m/\n$/;
    print $output_fh $s;
}

# my @files = find_summary_files($directory);

# sub find_summary_files {
#     my $d = shift;
#     opendir(my $dirh, $d)
# 	or die "unable to open directory $d: $!";
#     my @files = grep(/_Summary\.txt$/, readdir($dirh));
#     closedir($dirh);
#     return @files;
# }

# v("directory %s", $directory);
# c("Using");
# v("%-12s %f", "max-mp-ratio", $max_minus_plus_ratio);
# v("%-12s %f", "max-peak", $max_peak);
# v("%-12s %f", "L-max-mp-ratio", $L_max_minus_plus_ratio);
# v("%-12s %f", "L-max-peak", $L_max_peak);
# v("%-12s %d", "min-obs", $min_obs);
# v("%-12s %f", "majority", $super_majority);
# v("%-12s %d", "too-few", $too_few);
# c("");
# v("%-12s %s", "pattern-options",
#   join(" ", grep {! m/^#/ } read_file($flp_option_file, {chomp => 1})));
# c("");

# my @summaries;
# foreach my $f (@files) {
#     push @summaries, read_summary("$directory/$f");
# }
# foreach (sort {$b->{total_obs} <=> $a->{total_obs}} @summaries) {
#     process_summary($_);
# }

sub create_summary {
    my %summary = ();
    my $fn = shift;
    my $lines = shift;

    my @names = qw(id bx by bz mode obs_cluster obs_-cluster
		 box_cluster box_-cluster
		 density_mode sum_dens_cluster sum_dens_-cluster
		 vol_cluster vol_-cluster
		 mean_dist_mode max_dist_mode
		 peak_cluster
		 mean_dist_all_mode max_dist_all_mode
		 peak_all
		 );
    $summary{fn} = $fn;
    $summary{total_obs} = 0;
    $summary{clusters} = [];

    # klusternr
    # bx
    # by
    # bz
    # antal i mode
    # antal i kluster
    # antal i -kluster
    # antal bokse i kluster
    # antal bokse i -kluster
    # density i mode
    # sum of densitites i kluster
    # sum of densitites i -kluster,
    # volumen af kluster
    # volumen af -kluster
    # mean distance to mode
    # maximum distance to mode
    # peakness for kluster,
    # mean distance (all) to mode
    # maximum distance (all) to mode
    # peakness (all) for kluster
    for (@$lines) {
	chomp;
	my @fields = split;
	my %h;
	@h{@names} = @fields;
	warn "weird: empty cluster: $fn" unless defined $h{'obs_cluster'} && $h{'obs_cluster'} > 0;
	$summary{total_obs} += $h{'obs_cluster'};
	$summary{total_obs} += $h{'obs_-cluster'};
	push @{$summary{clusters}}, \%h;
    }
    return \%summary;
}

sub cluster_ok {
    my $h = shift;
    return 0 unless $h->{peak_cluster} < $max_peak;
    return 0 unless $h->{'obs_cluster'} + $h->{'obs_-cluster'} >= $min_obs;
    return 0 unless $h->{'obs_-cluster'} / $h->{'obs_cluster'} < $max_minus_plus_ratio;
    return 1;
}
sub L_cluster_ok {
    my $h = shift;
    return 0 unless $h->{peak_cluster} < $L_max_peak;
    return 0 unless $h->{'obs_-cluster'} / $h->{'obs_cluster'} < $L_max_minus_plus_ratio;
    return 0 unless $h->{'obs_cluster'} + $h->{'obs_-cluster'} >= $min_obs;
    return 1;
}

# Apply bonuses/penalties to the base quality.
sub adjust_quality {
    my ($qual, $summary) = @_;
    my $primary = $summary->{clusters}[0];

    # # OK, let's throw away everything else and use JEA's scheme. This
    # # means we need to make the quality some sufficient low number
    # # when @{$summary->{clusters}} > 1, and when
    # # @{$summary->{clusters}} == 1 we make it $LARGE_NUMBER -
    # # $mean_dist_all_mode.
    # #
    # # Then, when we read in all the assess files, we can add a very
    # # small bonus for being 'later', say $step_number/1000. We can't
    # # do that here, because we don't know the step number when
    # # processing each individual file. Ugly hack, yes.

    # yk (20170406): Just look at primary cluster.
    # return 100 - $primary->{mean_dist_all_mode};

    if (@{$summary->{clusters}} == 1) {
    	return 100 - $primary->{mean_dist_all_mode};
    }
    return 10;


    # Subtract penalty for large average distance to mode. The average
    # is bounded above by pi. Dividing by 10 gives a penalty roughly
    # in [0.0, 0.3]. We only apply this penalty in the 1-cluster case.
    $qual -= $primary->{mean_dist_all_mode}/10
	if @{$summary->{clusters}} == 1;

    # If we have a major good cluster, give a bonus for being really
    # good (that is, if we use a limit of 95%, but the big cluster is
    # actually 99%, we want to distinguish that from another step
    # where it was only 96%).
    if ($primary->{'obs_cluster'} / $summary->{total_obs} > $super_majority) {
	my $num = $primary->{'obs_cluster'} / $summary->{total_obs} - $super_majority;
	my $den = 1.0 - $super_majority;
	$qual += .5*$num/$den;
    }

    # For the low quality cases, we also want to promote use of a
    # guess when the largest cluster constitutes a large percentage of
    # all bonds. Let's say add one 10th of the percentage the largest
    # cluster constitutes.

    $qual += .1 * $primary->{'obs_cluster'} / $summary->{total_obs}
	if @{$summary->{clusters}} > 1;


    # TODO

    return $qual;
}

sub process_summary {
    my $summary = shift;
    my $fn = $summary->{fn};
    my $base = basename($fn);
    $base =~ m/(.*)_Summary2?\.txt/
	or die "unexpected base filename $base";
    my $pattern = $1;
    my $length = $2;
    my $res = $3;
    my $ok_func = ($length eq 'L') ? \&L_cluster_ok : \&cluster_ok;
    my $total = $summary->{total_obs};
    my $ncluster = scalar @{$summary->{clusters}};
    c("Processing %s", $fn);
    c("pattern %s, length %s", $pattern, $length);
    c("total observations: %d", $total);
    c("#clusters: %d", $ncluster);
    c("\tid\t#obs\tpeakness\tmp-ratio\tmode");
    # GRRR: Why the F*CK can $_->{'obs_cluster'} be 0?
    for (@{$summary->{clusters}}) {
	c("\t%d\t%d\t%f\t%f\t%d,%d,%d", $_->{id}, $_->{'obs_cluster'} + $_->{'obs_-cluster'},
	 $_->{peak_cluster}, $_->{'obs_cluster'} ? $_->{'obs_-cluster'} / $_->{'obs_cluster'} : "nan",
	 $_->{bx}, $_->{by}, $_->{bz});
    }

    my ($whitelist, $reason);
    my ($base_quality, $quality);


    # Now it's time to write a single line determining our assessment
    # of this pattern.  The stuff we need to output, beyond columns
    # containing the pattern and length, is:
    #
    # (a) Whether to whitelist all bonds having this pattern
    # (b) A simple keyword explaining (a)
    # (c) In case we end up whitelisting everything with this pattern,
    # also the mode box of the primary cluster (as that is what our
    # 'generated algorithm' will end up guessing at)
    #
    # (0) If there are fewer than $min_obs across all (plus and minus)
    # clusters (that is, we started with fewer than $min_obs), we
    # stop, but note that we don't really believe in this.
    #
    # (1) If there is precisely one cluster, and that cluster
    # satisfies all criteria, we're happy.
    #
    # (2) If there is more than one cluster, and all satisfy the
    # criteria, we're almost happy, but still want to try to refine
    # the pattern (that may not help, so it might have been better to
    # stop at this point and say 'the rotation is one of these 3
    # possibilities; you need to figure out some other way to pick the
    # right one').
    #
    # (3) Otherwise, we have at least one cluster not satisfying the
    # criteria, so we (?)  XXX: This sentence was apparently never
    # finished.
    #
    # Anyway, we should probably switch away from giving a textual
    # "reason" and instead give a numeric "quality" (we may however
    # still want a short text string describing the "most important"
    # factor in the determination of the quality, which might well be
    # the current "reason"). We probably want to keep the "whitelist"
    # separately - it is still rather useful to be able to ignore a
    # huge number of Hbonds in subsequent runs, to cut down on the
    # number of patterns encountered (and thus the number of
    # clustering runs needed).

    c('');
    if ($summary->{total_obs} < $too_few) {
	$whitelist = 1;
	$base_quality = 1;
	$reason = 'too_few';
	c("Assessment: Too few observations.");
    } elsif ($ncluster == 1 and $ok_func->($summary->{clusters}[0])) {
	$whitelist = 1;
	$reason = 'all_good';
	$base_quality = 10;
	c("Assessment: One good cluster.");
    } elsif ($ok_func->($summary->{clusters}[0]) &&
	     $summary->{clusters}[0]{'obs_cluster'} / $summary->{total_obs} >= $super_majority) {
	$whitelist = 0;
	# $whitelist = 1;
	$reason = "major_good_cluster";
	$base_quality = 8;
	c("Assessment: Primary cluster good, dominates all other.");
    } elsif ($ncluster > 1 and
	     all(map {$ok_func->($_)} @{$summary->{clusters}})) {
	$whitelist = 0;
	$reason = 'ambiguous';
	$base_quality = 4;
	c("Assessment: Multiple good clusters.");
    } elsif ($ncluster > 1 and
    	     $ok_func->($summary->{clusters}[0])) {
    	$whitelist = 0;
    	$reason = 'ambiguous';
    	$base_quality = 4;
    	c("Assessment: Multiple clusters (primary good), continue.");
    } else {
	$whitelist = 0;
	$reason = 'bad_cluster';
	$base_quality = 2;
	c("Assessment: Primary cluster bad, continue.");
	# c("Assessment: At least one bad cluster.");
    }

    $quality = adjust_quality($base_quality, $summary);

    $total{"${whitelist}${reason}"} += $summary->{total_obs};
    $count{"${whitelist}${reason}"}++;

    my $primary = $summary->{clusters}[0];
    p("%s\t%s\t%d,%d,%d\t%.2f;%.2f,%.2f,%.2f\t%.4f" . ("\t%.2f"x6),
      "${pattern}_${length}${res}", $length,
      $primary->{bx}, $primary->{by}, $primary->{bz},
      unrotated_center($primary->{bx}, $primary->{by}, $primary->{bz}),
      $quality,
      map {$primary->{$_}} qw|mean_dist_mode     max_dist_mode     peak_cluster
			      mean_dist_all_mode max_dist_all_mode peak_all|);

    c('');
}

# foreach (sort keys %total) {
#     v("%s-bonds %d", $_, $total{$_});
#     v("%s-patterns %d", $_, $count{$_});
# }

sub all { $_ || return 0 for @_; 1 }

sub unrotated_center {
    my ($bx, $by, $bz) = @_;
    my $mat = undo_magic_rotation(map {idx_to_coord($_)} ($bx, $by, $bz));
    return SO3_to_EV($mat);
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
    return (0, 1, 0, 0) if $angle == 0;

    if (sin($angle) > 0.0001) {
	# The normal case: Angle is away from 0 and pi.
	my $x = $M->get_elem(2, 1) - $M->get_elem(1, 2);
	my $y = $M->get_elem(0, 2) - $M->get_elem(2, 0);
	my $z = $M->get_elem(1, 0) - $M->get_elem(0, 1);
	my $n = sqrt($x*$x + $y*$y + $z*$z);
	if (abs($n-2*sin($angle)) > 0.001) {
	    printf STDERR "Unexpected large diff (%f) between |v| and 2sin(angle)\n", abs($n-2*sin($angle));
	}
	$x /= $n;
	$y /= $n;
	$z /= $n;
	return ($angle, $x, $y, $z);
	# return ($x, $y, $z);
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
	my $n = sqrt($x*$x + $y*$y + $z*$z);
	$x /= $n;
	$y /= $n;
	$z /= $n;
	return ($angle, $x, $y, $z);
	# return ($x, $y, $z);
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

