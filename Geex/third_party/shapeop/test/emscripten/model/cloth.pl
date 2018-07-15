#!/usr/bin/perl

$w = $ARGV[0];
$h = $ARGV[1];

# Create vertices
foreach $y (0 .. $h) {
	foreach $x (0 .. $w) {
		print "v $x $y 0\n";
	}
}

# Create faces
foreach $y (0 .. $h - 1) {
	foreach $x (0 .. $w - 1) {
		$ll = $x + $y * ($w + 1) + 1;
		$lr = $ll + 1;
		$ul = $ll + $w + 1;
		$ur = $ul + 1;
		print "f $ll $lr $ur $ul\n";
	}
}
