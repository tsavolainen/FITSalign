#!/usr/bin/perl -s
#
################################################################
# FITSalign.pl v1.0 - A program to align VLBI images
# Copyright (C) 2017 Tuomas Savolainen
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License, see
#  <https://www.gnu.org/licenses/>.
#
################################################################
#
# usage: FITSalign.pl [-aco] <image1> <image2> [<x1> <y1> <x2> <y2>
#                     <xtemp1> <ytemp1> <xtemp2> <ytemp2>] 
#
# The purpose of this program is to align two FITS images 
# by normalized 2D cross correlation. The input images need to 
# have the same pixel size and their PSFs (or restoring beams) 
# should be identical. 
#
# The program first plots the input images and asks the user
# to select an appropriate area for the subsequent analysis.
# Then the user is prompted to select a feature in the image to be 
# matched. In case of aligning radio images of AGN jets taken at 
# different frequencies, it is important to select a feature 
# that is optically thin and does not show strong spectral gradients. 
# 
# The program calculates a normalized cross-correlation between the 
# first image and the selected feature of the second image using a 
# fast method described by J.P. Lewis (1995). The correlation plot is 
# shown along with useful information about the shift between the 
# input images (calculated from the maximum of the correlation plot).  
#
# The program also calculates two spectral index images: one 
# with the calculated shift applied and one without. Check the 
# results! The program gives a warning if the feature selected 
# for matching turns out to contain optically thick parts. Finally, it 
# is possible to write out the spectral index image (with the shift 
# applied) as well as the correlation plot as FITS files.
#
# Switches:
#
#   -o=<outfile>   If <outfile> is given, the output is written to 
#                  <outfile>.out and the spectral index image to 
#                  <outfile>.fits.  
#
#   -c=<icut>      Parameters <icut> and is the cut-off level used 
#                  in calculating a spectral index maps. They are 
#                  optional: if not specified by the user, a value 
#                  of 0.001 Jy/beam is used. 
#
#   -a             Automatic mode: user still needs to provide the 
#                  coordinates of the image region and template. No 
#                  output graphics.
#
#
# NOTE THAT BAD VALUE SUPPORT IN PDL IS NEEDED TO RUN THIS 
# PROGRAM.
#
# - Tuomas Savolainen 
#   tuomas.k.savolainen at aalto.fi
#
################################################################

use PDL;
use PDL::NiceSlice;
use PDL::IO::FITS;
use PDL::Image2D;
use PDL::ImageND;
use PDL::Graphics::PGPLOT::Window; 
use PDL::Graphics::LUT;
use PDL::Math;
use strict;

### Check that bad value support is installed ###
if ($PDL::Bad::Status == 0) {
    die "Bad value support is not installed!\n";
}


### Initialize hashes ###
our %plotopt;
our %factor;
our %output;

### Set up constants ###
my $pi = 4*atan2(1,1);
my $L = 10;

### Set up plot options ###
$plotopt{'units'} = 'mas';
$plotopt{'icut1'} = 0.001;
$plotopt{'icut2'} = 0.001;


### Set up conversion factor ###
$factor{'mas'}     = 3.6e6;  
$factor{'arcsec'}  = 3.6e3;
$factor{'arcmin'}  = 60;
$factor{'final'} = 1;

### Set up output files ###
$output{'result'} = 'alignment.out';
$output{'specim'} = '';
$output{'nccim'} = '';
$output{'automode'} = '';

### Parse switches  ###
our ($c, $o, $G, $a);
if ($c) {
    $plotopt{'icut1'} = $c;
    $plotopt{'icut2'} = $c;
}
if ($o) {
    $output{'result'} = "$o" .'.out';
    $output{'specim'} = "$o" .'.si.fits';
    $output{'nccim'} = "$o" .'.ncc.fits';
}
if ($a) {
    $output{'automode'} = $a;
    if (!$o) {
	$output{'specim'} = 'alignment.si.fits';
	$output{'nccim'} = 'alignment.ncc.fits';
    }	
    if (@ARGV < 10) {
	die "ERROR: Not enough arguments in the automatic mode\n";
    }
}


### Check the number of given arguments ### 
if (@ARGV < 2) {
    die "ERROR: Not enough arguments\n";
}


### Load FITS images ###
my $file1 = shift @ARGV; 
my $file2 = shift @ARGV;
print "Loading FITS images...\n";
my ($imag1,$freq1) = loadfits($file1); 
my ($imag2,$freq2) = loadfits($file2);
my $fqtit1 = int(rint($freq1/1e8))/10; # to GHz with one decimal
my $fqtit2 = int(rint($freq2/1e8))/10; # to GHz with one decimal
my $title1 = "$fqtit1 GHz";
my $title2 = "$fqtit2 GHz";

### Check pixelsize and beam size ###
my $h1 = $imag1->gethdr;
my $h2 = $imag2->gethdr;
my $degppix1 = $$h1{CDELT2};
my $npx1 = $$h1{NAXIS1};
my $bmaj1 = $$h1{BMAJ};
my $bmin1 = $$h1{BMIN};
my $bpa1 = $$h1{BPA};
my $degppix2 = $$h2{CDELT2};
my $npx2 = $$h2{NAXIS1};
my $bmaj2 = $$h2{BMAJ};
my $bmin2 = $$h2{BMIN};
my $bpa2 = $$h2{BPA};
if ($npx1 != $npx2) {
    die "ERROR: The input images do not have the same number of pixels\n";
}
if (abs($degppix1-$degppix2) > 1e-11) {
    die "ERROR: The pixel increment in mas is not the same in the two input figures\n";
}
if (abs($bmaj1-$bmaj2) > 1e-9 | abs($bmin1-$bmin2) > 1e-9 | abs($bpa1-$bpa2) > 0.01) {
    die "ERROR: The input figures do not have the same beam sizes\n";
}
my $beamx = $bmaj1*$bmin1/sqrt(($bmaj1*sin($bpa1*$pi/180-$pi/2))**2+($bmin1*cos($bpa1*$pi/180-$pi/2))**2);
my $beamy = $bmaj1*$bmin1/sqrt(($bmaj1*sin($bpa1*$pi/180-0))**2+($bmin1*cos($bpa1*$pi/180-0))**2);
print "...ok\n";

### Plot and crop the images ###
my ($win1, $win2, $win3, $win4, $x1, $x2, $y1, $y2);
if ($output{'automode'} > 0) {
    $x1 = shift @ARGV;
    $y1 = shift @ARGV;
    $x2 = shift @ARGV;
    $y2 = shift @ARGV;
    $imag1 = crop($imag1,$x1,$y1,$x2,$y2);
    $imag2 = crop($imag2,$x1,$y1,$x2,$y2);
} else {
    print "Select an image area to be used in the correlation analysis\n"; 
    $win1 = plot_im($imag1,$title1,0);
    $win1->hold();
    $win2 = plot_im($imag2,$title2,0);
    $win2->hold();
    ($x1,$y1,$x2,$y2) = sel_area($win1,$imag1);
    $x1 = rint($x1); $y1 = rint($y1);
    $x2 = rint($x2); $y2 = rint($y2);
    $imag1 = crop($imag1,$x1,$y1,$x2,$y2);
    $imag2 = crop($imag2,$x1,$y1,$x2,$y2);
    $win1->close();
    $win2->close();
    print "Coordinates of the selected area:\n";
    print "x1 = $x1 , y1 = $y1 , x2 = $x2 , y2 = $y2 \n";
}

### Select matched feature ###
my ($temp,$xt1, $xt2, $yt1, $yt2);
if ($output{'automode'} > 0) {
    $xt1 = shift @ARGV;
    $yt1 = shift @ARGV;
    $xt2 = shift @ARGV;
    $yt2 = shift @ARGV;
    $temp = crop($imag2,$xt1,$yt1,$xt2,$yt2);
} else {
    print "Select a feature to be matched\n"; 
    $win1 = plot_im($imag1,$title1,0);
    $win1->hold();
    $win2 = plot_im($imag2,$title2,0);
    $win2->hold();
    my $test = 0;
    while ($test == 0) {
	($xt1,$yt1,$xt2,$yt2) = sel_area($win2,$imag2);
	$xt1 = rint($xt1); $yt1 = rint($yt1);
	$xt2 = rint($xt2); $yt2 = rint($yt2);
	if ($xt1 > $L && $yt1 > $L && $x2-$x1-$xt2 > $L && $y2-$y1-$yt2 > $L) {
	    $test = 1;
	}
	if ($test == 0) {
	    print "Select a feature that is at least $L pixels \n";
	    print "away from the edges of the image.\n";
	}
    }
    $temp = crop($imag2,$xt1,$yt1,$xt2,$yt2);
    $win2->rect($xt1,$xt2,$yt1,$yt2,{FILLTYPE=>'Outline'});
    my $xt11 = $xt1+$x1; my $xt21 = $xt2+$x1;
    my $yt11 = $yt1+$y1; my $yt21 = $yt2+$y1;
    print "Coordinates of the selected feature:\n";
    print "x1 = $xt11 , y1 = $yt11 , x2 = $xt21 , y2 = $yt21 \n";
}
my $xtc = $xt1+floor(($xt2-$xt1)/2);
my $ytc = $yt1+floor(($yt2-$yt1)/2);


### Run correlation ###
print "Calculating cross-correlation... \n";
my ($ncc_im) = ncc2d($imag1,$temp);


### Show correlation plot ###
my $convert = $factor{'final'};
if (!$output{'automode'}) {
    $win2->close();
    $win1->close();
    my $title_ncc = "Normalized Cross Correlation between images at $title1 and $title2";
    $win1 = PDL::Graphics::PGPLOT::Window->new(Dev=>'/xserve'); 
    $win1->ctab(lut_data('rainbow4',0,'ramp'));
    $win1->imag($ncc_im, {Justify => 1, DrawWedge=>1});
    $win1->hold();
    $win1->label_axes({Font=>2,XTitle=>"Relative R.A. ($plotopt{'units'})", YTitle=>"Relative Declination ($plotopt{'units'})", Title => $title_ncc});
}

### Process results ###
my ($peak,$x0,$y0) = max2d_ind($ncc_im);
$peak = $peak->sclr; $x0 = $x0->sclr; $y0 = $y0->sclr;
$peak = rint($peak*1000)/1000;
my $dpix_x = $x0-$xtc; 
my $dpix_y = $y0-$ytc;
my $dRA = -$dpix_x*$convert; 
my $dDec = $dpix_y*$convert;
$beamx = $beamx*$factor{$plotopt{'units'}}; $beamy = $beamy*$factor{$plotopt{'units'}};
my $dRAbeam = $dRA/$beamx; my $dDecbeam = $dDec/$beamy;
$dRA = rint(1000*$dRA)/1000; $dDec = rint(1000*$dDec)/1000;
$dRAbeam = rint(1000*$dRAbeam)/1000; $dDecbeam = rint(1000*$dDecbeam)/1000;

### Print results ###
print "######################################################################\n";
print "Highest correlation (r = $peak ) is obtained with the following shift:\n";
print "dx = x_image1 - x_image2 = $dpix_x pixels\n";
print "dy = y_image1 - y_image2 = $dpix_y pixels\n";
print "or in relative sky coordinates:\n";
print "dRA = RA_image1 - RA_image2 = $dRA $plotopt{'units'}\n";
print "dDec = Dec_image1 - Dec_image2 = $dDec $plotopt{'units'}\n";
print "######################################################################\n";

### Write out results ###
open OUTPUT, ">", $output{'result'};
my $form = "%6.0f" . "%6.0f" . "%10.3f" . "%10.3f". "%10.3f" . "%10.3f";
select OUTPUT;
    print  "dx(pix) dy(pix) dRA(mas) dDec(mas) dRA(beam) dDec(beam)\n";
    printf "$form", $dpix_x, $dpix_y, $dRA, $dDec, $dRAbeam, $dDecbeam;
close OUTPUT;
select STDOUT;


### Load original images ###
print "Calculating spectral index images... ";
($imag1,$freq1) = loadfits($file1); 
($imag2,$freq2) = loadfits($file2);
$fqtit1 = int(rint($freq1/10**8))/10; # to GHz with one decimal
$fqtit2 = int(rint($freq2/10**8))/10; # to GHz with one decimal
$title1 = "$fqtit1 GHz";
$title2 = "$fqtit2 GHz";


### Shift image 2 ###
my $imag2_shift;
$imag2_shift = shift_im($imag2,-$dpix_x,-$dpix_y);


### Calculate spectral index maps ###
my $si_im1 = alpha_im($imag1,$imag2,$plotopt{'icut1'},$plotopt{'icut2'});
my $si_im2 = alpha_im($imag1,$imag2_shift,$plotopt{'icut1'},$plotopt{'icut2'});
print "ok\n";



if (!$output{'automode'}) {
    ### Plot spectral index map 1 ###
    my $title_si1 = "Spectral index image between $title1 and $title2 WITHOUT shift";
    $win3 = PDL::Graphics::PGPLOT::Window->new(Dev=>'/xserve'); 
    my $t = $win3->transform(dims($si_im1), {ImageCenter => 0,  Pixinc => $convert});
    $t(0) .= -$t(0); $t(1) .= -$t(1);
    $win3->ctab(lut_data('rainbow4',0,'ramp'));
    $win3->imag($si_im1, {Justify => 1, Transform => $t, DrawWedge=>1});
    $win3->hold();
    $win3->label_axes({Font=>2,XTitle=>"Relative R.A. ($plotopt{'units'})", YTitle=>"Relative Declination ($plotopt{'units'})", Title => $title_si1});
    ### Plot spectral index map 2 ###
    my $title_si2 = "Spectral index image between $title1 and $title2 WITH shift applied";
    $win4 = PDL::Graphics::PGPLOT::Window->new(Dev=>'/xserve'); 
    $t = $win4->transform(dims($si_im2), {ImageCenter => 0,  Pixinc => $convert});
    $t(0) .= -$t(0); $t(1) .= -$t(1);
    $win4->ctab(lut_data('rainbow4',0,'ramp'));
    $win4->imag($si_im2, {Justify => 1, Transform => $t, DrawWedge=>1});
    $win4->hold();
    $win4->label_axes({Font=>2,XTitle=>"Relative R.A. ($plotopt{'units'})", YTitle=>"Relative Declination ($plotopt{'units'})", Title => $title_si2});
}

### Check that the matched feature did not include any optically thick parts ###
my $xts1 = $x1+$xt1-$dpix_x; my $yts1 = $y1+$yt1-$dpix_y;
my $xts2 = $x1+$xt2-$dpix_x; my $yts2 = $y1+$yt2-$dpix_y;
if ($xts1 < 0) { $xts1 = 0; } 
if ($yts1 < 0) { $yts1 = 0; }
if ($xts2 > $npx1-1) { $xts2 = $npx1-1; }
if ($yts2 > $npx1-1) { $yts2 = $npx1-1; }
my $temp_chk = crop($si_im2,$xts1,$yts1,$xts2,$yts2);
$temp_chk = $temp_chk->badmask(-100);
if (max($temp_chk) > 0) {
    print "############################################################\n";
    print "WARNING!!! Matched feature contains optically thick parts!!!\n";
    print "You may have to re-run the matching.\n";
    print "############################################################\n";
} 


### Write out spectral index image ###
if ($output{'automode'} > 0) {
    print "Writing out FITS...\n";
    $si_im2->wfits($output{'specim'},-32); # 32-bit output
} else {
    my $ans;
    print "Do you want to write out the SHIFTED spectral index image in FITS format? (y/N) ";
    chomp($ans = <STDIN>);
    if ($ans eq "y" | $ans eq "Y" | $ans eq "yes" | $ans eq "YES") {
	if ($output{'specim'}) {
	    print "Writing out FITS...\n";
	    $si_im2->wfits($output{'specim'},-32); # 32-bit output
	} else {
	    print "Give a filename: ";
	    my $outfile = <STDIN>;
	    chomp($outfile);
	    print "Writing out FITS...\n";
	    $si_im2->wfits("$outfile",-32); # 32-bit output
	}
    }
}


### Write out correlation image ###
if ($output{'automode'} > 0) {
    print "Writing out FITS...\n";
    $ncc_im->wfits($output{'nccim'},-32); # 32-bit output
} else {
    my $ans;
    print "Do you want to write out the correlation image in FITS format? (y/N) ";
    chomp($ans = <STDIN>);
    if ($ans eq "y" | $ans eq "Y" | $ans eq "yes" | $ans eq "YES") {
	if ($output{'nccim'}) {
	    print "Writing out FITS...\n";
	    $ncc_im->wfits($output{'nccim'},-32); # 32-bit output
	} else {
	    print "Give a filename: ";
	    my $outfile = <STDIN>;
	    chomp($outfile);
	    print "Writing out FITS...\n";
	    $ncc_im->wfits("$outfile",-32); # 32-bit output
	}
    }
}

if (!$a) {
    $win1->close();
    $win3->close();
    $win4->close();
}

####################
### Subroutines ####
####################

### Loads FITS file ###
sub loadfits{
    my $file  = $_[0];
    my $im = rfits("$file");
    my $head = $im->gethdr;

    if ($$head{CDELT2} eq undef | $$head{CRVAL3} eq undef) {
	die "Missing header information!\n"
    }

    my $degpp = $$head{CDELT2};
    my $freq = $$head{CRVAL3};

    # Determine plot size units based on cellspacing
    if ($degpp < 1.0/360000.0)  {$plotopt{'units'}='mas';}  
    elsif ($degpp < 1.0/6000.0) {$plotopt{'units'}='arcsec';}  
    else  {$plotopt{'units'}='arcmin';}                         

    # Determine conversion
    my $convert = $degpp*$factor{$plotopt{'units'}};
    $factor{'final'} = $convert;

    return ($im, $freq);
}

### Plots image to a new PGPLOT windows ###
sub plot_im {
    my $im  = $_[0];
    my $title = $_[1];
    my $scl = $_[2];
    my $contours = $_[3];
    my $win = PDL::Graphics::PGPLOT::Window->new(Dev=>'/xserve');  
    $win->ctab(lut_data('rainbow4',0,'ramp'));
    if ($scl == 0) {
	$win->imag($im, {Justify => 1, ITF=>'LOG', RANGE=>[0,max($im)], DrawWedge=>1});
	$win->hold();
	$win->label_axes({Font=>2,XTitle=>"X (pixels)", YTitle=>"Y (pixels)", Title =>$title});
    } elsif ($scl == 1) {
	my $convert = $factor{'final'};
	my $t = $win->transform(dims($im), {ImageCenter => 0,  Pixinc => $convert});
	$t(0) .= -$t(0); $t(1) .= -$t(1);
	$win->imag($im, {Justify => 1, Transform => $t, ITF=>'LOG', RANGE=>[0,max($im)], DrawWedge=>1});
	$win->hold();
	$win->label_axes({Font=>2,XTitle=>"Relative R.A. ($plotopt{'units'})", YTitle=>"Relative Declination ($plotopt{'units'})", Title =>$title});
    } elsif ($scl == 2) {
	my $convert = $factor{'final'};
	my $t = $win->transform(dims($im), {ImageCenter => 0,  Pixinc => $convert});
	$t(0) .= -$t(0); $t(1) .= -$t(1);
	$win->imag($im, {Justify => 1, Transform => $t, ITF=>'LOG', RANGE=>[0,max($im)], DrawWedge=>1});
	$win->hold();
	$win->cont($im, {CONTOURS=>$contours, Transform => $t});
	$win->label_axes({Font=>2,XTitle=>"Relative R.A. ($plotopt{'units'})", YTitle=>"Relative Declination ($plotopt{'units'})", Title =>$title});
    }
    return $win;
}

### Normalized 2D cross correlation ###
### Uses the sum table method of Lewis (1995) ###
sub ncc2d {
    my $im = $_[0]; 
    my $temp = $_[1];

    # Calculate the normalized template and its energy
    $temp = $temp-avg($temp);
    my $temp_var = sum(($temp)**2);

    # Calculate cross correlation using fast method
    my $corr = convolveND($im, $temp(-1:0,-1:0));

    # Calculate sums of the image
    my $m_im = $im->dim(0);
    my $n_im = $im->dim(1);
    my $m_temp = $temp->dim(0);
    my $n_temp = $temp->dim(1);
    my $mc_temp = floor(($m_temp-1)/2);
    my $nc_temp = floor(($n_temp-1)/2);
    my $s = zeroes($m_im,$n_im);
    my $s2 = zeroes($m_im,$n_im);
    my ($u,$v);

    for (my $u = 0; $u < $m_im; $u++) {
	for (my $v = 0; $v < $n_im; $v++) {
	    if ($u-1 < 0 && $v-1 < 0) {
		$s($u,$v) .= $im($u,$v);
	    } elsif ($u-1 < 0) {
		$s($u,$v) .= $im($u,$v)+$s($u,$v-1);
	    } elsif ($v-1 < 0) {
		$s($u,$v) .= $im($u,$v)+$s($u-1,$v);
	    } else {
		$s($u,$v) .= $im($u,$v)+$s($u-1,$v)+$s($u,$v-1)-$s($u-1,$v-1);
	    }
	}
    }

    for (my $u = 0; $u < $m_im; $u++) {
	for (my $v = 0; $v < $n_im; $v++) {
	    if ($u-1 < 0 && $v-1 < 0) {
		$s2($u,$v) .= $im($u,$v)**2;
	    } elsif ($u-1 < 0) {
		$s2($u,$v) .= $im($u,$v)**2+$s2($u,$v-1);
	    } elsif ($v-1 < 0) {
		$s2($u,$v) .= $im($u,$v)**2+$s2($u-1,$v);
	    } else {
		$s2($u,$v) .= $im($u,$v)**2+$s2($u-1,$v)+$s2($u,$v-1)-$s2($u-1,$v-1);
	    }
	}
    }

    my $S_f = zeroes($m_im,$n_im);
    my $e_f = zeroes($m_im,$n_im);
    my $u_min = $mc_temp;
    my $v_min = $nc_temp;
    my $u_max = $m_im-($m_temp-$mc_temp);
    my $v_max = $n_im-($n_temp-$nc_temp);
    
    foreach($u_min..$u_max) {
	my $u = $_;
	foreach($v_min..$v_max) {
	    my $v = $_;
	    if ($u-$mc_temp-1 < 0 && $v-$nc_temp-1 < 0) {
		$S_f($u,$v) .= $s($u-$mc_temp+$m_temp-1,$v-$nc_temp+$n_temp-1);
		$e_f($u,$v) .= $s2($u-$mc_temp+$m_temp-1,$v-$nc_temp+$n_temp-1);
	    } elsif ($u-$mc_temp-1 < 0) {
		$S_f($u,$v) .= $s($u-$mc_temp+$m_temp-1,$v-$nc_temp+$n_temp-1)-$s($u-$mc_temp+$m_temp-1,$v-$nc_temp-1);
		$e_f($u,$v) .= $s2($u-$mc_temp+$m_temp-1,$v-$nc_temp+$n_temp-1)-$s2($u-$mc_temp+$m_temp-1,$v-$nc_temp-1);
	    } elsif ($v-$nc_temp-1 < 0) {
		$S_f($u,$v) .= $s($u-$mc_temp+$m_temp-1,$v-$nc_temp+$n_temp-1)-$s($u-$mc_temp-1,$v-$nc_temp+$n_temp-1);
		$e_f($u,$v) .= $s2($u-$mc_temp+$m_temp-1,$v-$nc_temp+$n_temp-1)-$s2($u-$mc_temp-1,$v-$nc_temp+$n_temp-1);
	    } else {
		$S_f($u,$v) .= $s($u-$mc_temp+$m_temp-1,$v-$nc_temp+$n_temp-1)-$s($u-$mc_temp-1,$v-$nc_temp+$n_temp-1)-$s($u-$mc_temp+$m_temp-1,$v-$nc_temp-1)+$s($u-$mc_temp-1,$v-$nc_temp-1);
		$e_f($u,$v) .= $s2($u-$mc_temp+$m_temp-1,$v-$nc_temp+$n_temp-1)-$s2($u-$mc_temp-1,$v-$nc_temp+$n_temp-1)-$s2($u-$mc_temp+$m_temp-1,$v-$nc_temp-1)+$s2($u-$mc_temp-1,$v-$nc_temp-1);
	    }
	}
    }
    
    # Normalize cross correlation
    my $S_f_avg = $S_f/($m_temp*$n_temp);
    my $E_f = $e_f-2*$S_f_avg*$S_f+$m_temp*$n_temp*$S_f_avg**2;
    my $norm = sqrt($E_f*$temp_var);
    my $ncc = $corr/$norm;
    $ncc = $ncc->badmask(0);

    return ($ncc, $corr, $norm);

}

### Subroutine to select an area from the image with a cursor ###
sub sel_area {
    my $win = $_[0];
    my $im = $_[1];
    my @dim = $im->dims;
    my $test = 0;
    my ($xc,$yc,$ch,$xref,$yref);
    while ($test == 0) {
	($xref,$yref,$ch) = $win->cursor();
	($xc,$yc,$ch,$xref,$yref) = $win->cursor({XRef=>$xref,YRef=>$yref,Type=>'Rectangle'});
	if ($xc >= -0.5 && $xc < $dim[0]-0.5 && $yc >= -0.5 && $yc < $dim[1]-0.5 && $xref >= -0.5 && $xref < $dim[0]-0.5 && $yref >= -0.5 && $yref < $dim[1]-0.5) {
	    $test = 1
	}
    }
    if ($xc > $xref) {
	my $x1 = $xref;
	my $x2 = $xc;
	$xc = $x1;
	$xref = $x2;
    }
    if ($yc > $yref) {
	my $y1 = $yref;
	my $y2 = $yc;
	$yc = $y1;
	$yref = $y2;
    }

    return ($xc,$yc,$xref,$yref);
}

### Crops image ###
sub crop {
    my $im = $_[0];
    my $x1 = $_[1];
    my $y1 = $_[2];
    my $x2 = $_[3];
    my $y2 = $_[4];
    my $xmin = 0;
    my $ymin = 0;

    # Crop the image
    if ($x1 > $x2) {
	$xmin = $x2;
	if ($y1 > $y2) {
	    $im = $im($x2:$x1,$y2:$y1);
	    $ymin = $y2;
	} else {
	    $ymin = $y1;
	    $im = $im($x2:$x1,$y1:$y2);
	}
    } else {
	$xmin = $x1;
	if ($y1 > $y2) {
	    $im = $im($x1:$x2,$y2:$y1);
	    $ymin = $y2;
	} else {
	    $im = $im($x1:$x2,$y1:$y2);
	    $ymin = $y1;
	}
    }

    # Fix the header
    my @dim = $im->dims;
    my $head = $im->gethdr;
    my $centx = $$head{CRPIX1}; 
    my $centy = $$head{CRPIX2};       
    $$head{NAXIS1} = $dim[0]; 
    $$head{NAXIS2} = $dim[1];
    $$head{CRPIX1} = $centx-$xmin;
    $$head{CRPIX2} = $centy-$ymin;
    $im->sethdr($head);
	
    # Return cropped image
    return $im;
}

### Subroutine to linearly shift image. Uses interpolation. ###
### Adapted from Dan Homan's FITSplot.pl                    ###
sub shift_im {
    my $im = $_[0];
    my $xshift = $_[1]; #X-shift in pixels
    my $yshift = $_[2]; #Y-shift in pixels

    my $head = $im->gethdr;

    my $px = pdl([$xshift,1.0],[0.0,0.0]);
    my $py = pdl([$yshift,0.0],[1.0,0.0]);
    my $temp = warp2d($im,$px,$py, {KERNEL=>"tanh", NOVAL=>0});

    $temp->sethdr($head);
    
    return $temp;
}

### Subroutine to calculate a spectral index image ###
sub alpha_im {
    my $im1 = $_[0];
    my $im2 = $_[1];
    my $icut1  = $_[2];   # cut level for im1
    my $icut2  = $_[3];   # cut level for im2

    my $head1 = $im1->gethdr;
    my $head2 = $im2->gethdr;

    my $freq1 = $$head1{CRVAL3};
    my $freq2 = $$head2{CRVAL3};
    
    $im1 = $im1->setbadif($im1 < $icut1 | $im2 < $icut2);
    $im2 = $im2->setbadif($im2 < $icut2 | $im1 < $icut1);

    my $si = log($im1/$im2)/log($freq1/$freq2);

    # Fix header #
    $$head2{CTYPE4} = 'ALPHA';
    $$head2{BUNIT} = '';
    $$head2{DATAMIN} = min($si);
    $$head2{DATAMAX} = max($si);
    $si->sethdr($head2);

    return $si;
}

