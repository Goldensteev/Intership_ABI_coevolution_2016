#! /usr/bin/perl

# generate histograms

# created : ?1999
# modified: 4 nov 2002 - add  negative values
# modified: 6 nov 2002 - add the option of frequency
# modified: 11 nov 2002 - small bug
# modified: may 2004    - add the cumulative '-c' option
# modified: june 2006   - correct bug for int negative values and add the reverse cumulative '-r' option
# modified: nov 2008    - add the -0 option that do not report null values

# author Guillaume Achaz <gachaz@gmail.com>

use Getopt::Std;
getopts('hu:n:fcr0');


#$# = 3;
$VERSION = "1.5 (Nov 08)";


# print help and exit
if ( $opt_h ){
	print STDERR "**  version $VERSION  **\n";
	print STDERR "$0 [-h] [-f] [-c|-r] {-u #}|{-n #} [file]\n";
	print STDERR "\t-f : give frequency instead of counts\n";
	print STDERR "\t-c : do a cumulative histogram\n";
	print STDERR "\t-r : do a reverse cumulative histogram\n";
	print STDERR "\t-0 : do NOT report null values\n";
	print STDERR "\t-u : unity for histogram construction\n";
	print STDERR "\t-n : number of bins\n";
	print STDERR "\t     -n and -u are mutually exclusive\n";
	exit (1);
}

# remember user the basic syntax
if(  (!$opt_u && !$opt_n)  ||  ($opt_u && $opt_n) || ($opt_c && $opt_r) ){
	print STDERR "$0 [-h] [-f] [-z] [-c|-r] {-u #}|{-n #} [file]\n";
	exit (1);
}



$count = 0;
$max = 0; $max_h = 0;
$min = 1e10; $min_h = 1e10;

if ($opt_u){$unite = $opt_u;}   # l'unité est de tant (le pas)
if ($opt_n){$nbe = $opt_n;}     # nbe de case


while (<>){
    chop;
    $var[$count] = $_;                              # la var est la seule valeur
    $max = ($var[$count]>$max)?$var[$count]:$max;   # get max;
    $min = ($var[$count]>$min)?$min:$var[$count];   # get min;

    $count ++;
  }


$unite = ( $opt_n ) ? ($max-$min)/$nbe : $unite; # pour avoir $nbe pts

$decalage = ($min<0)? 1-int( $min / $unite ):0;

#print STDERR "deca is $decalage\n";

for ($w=0; $w < $count; $w++){
   
    $hist = ($var[$w] >= 0 || ($var[$w] < 0 && $var[$w] == int($var[$w]) ) )? int( $var[$w] / $unite ) :  int( $var[$w] / $unite )-1;
	 
    $max_h = ($hist>$max_h)?$hist:$max_h;
    $min_h = ($hist>$min_h)?$min_h:$hist;	 
    $tablo[ $hist + $decalage ]++;          # decalage is used for nagtive values;

}


$nbe_affiche=($opt_u)?$#tablo+1:"$nbe+1";  # calcule le nombre de case


print STDERR "# histo.pl $VERSION \n";
print STDERR "# $count values from $min to $max\n";
print STDERR "# $nbe_affiche bins of $unite units\n";

if($opt_c ){
	$sum=0;
}
if( $opt_r){
	if($opt_f){
		$sum=1.0;
	}
	else{
		$sum=$count;
	}
}

for ($z=$min_h; $z<=$max_h; $z++){

	if($tablo[ $z+$decalage ] == 0){
		$tablo[$z+ $decalage]=0;   # I do not remember what it is good for but I leave it as it is -nov 08- ? :-)
	}
	
	if($opt_0 && $tablo[ $z+$decalage ] == 0){
		next;
	}

	$case = $z*$unite ;
	
	if($opt_f){
		
		$freq = $tablo[ $z+$decalage ]/$count;
		if( $opt_c ){
			$sum+=$freq;    
			printf("%g\t%.5f\n", $case,$sum); 
		}
		elsif( $opt_r ) {
			printf("%g\t%.5f\n", $case,$sum);
			$sum-=$freq; 
		}
		else { 
			printf("%g\t%.5f\n", $case,$freq);
		}
	}
	else{
		if( $opt_c ){
			$sum+=$tablo[$z+$decalage];
			printf("%g\t%d\n", $case,$sum);
		}
		elsif( $opt_r ){
			printf("%g\t%d\n", $case,$sum);
			$sum-=$tablo[$z+$decalage];
		}
		else {
			print "$case\t$tablo[ $z+$decalage ]\n";
		}
	}
}

exit(0);
