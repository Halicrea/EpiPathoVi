#!/usr/bin/perl

# Copyright INRA

# Jerome.Gouzy@inra.fr
# Erika.Sallet@inra.fr
# Sebastien.Carrere@inra.fr
# Ludovic.Cottret@inra.fr
# Ludovic.Legrand@inra.fr

# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".

# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.

# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.

# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

=head1 NAME

lipm_findmotif.pl - 

=head1 DESCRIPTION

 Find motif in DNA sequence. Output GFF3 and stats

=head1 SYNOPSIS

 lipm_findmotif.pl -fasta genome.fasta -motif GATC

 lipm_findmotif.pl -fasta genome.fasta -motif motif.csv

=head1 ARGUMENTS

=over 8

=item B<fasta> : genome fasta file

=item B<motif> : file of motifs or motif string

 format: <motif>\t<type>\t<modification position 0-based>
 
 ie: GATC	m6A	 1                         # palindromic
     TCAYNNNNNNTCC/GGANNNNNNRTGA	m6A	2  # not palindromic

=back

=head1 OPTIONS

=over 8

=item B<site_only> : print only modification position. Default site+motif.

=item B<motif_only> : print only the motif position. Default site+motif

=item B<stats> : compute motif stats

=item B<output> : output directory. Default: current directory

=item B<prefix> : prefix output file. Default: processid

=back

=head1 FUNCTIONS AND PROCEDURES

=cut

use strict;
use warnings;
use Data::Dumper;

use Pod::Usage;
use File::Basename;
use File::Spec::Functions;
use Cwd 'abs_path';
use ParamParser;

my $DEBUG = 0;

our %IUPAC_DNA = (
				   N   => '[ATGC]',
				   W   => '[AT]',
				   S   => '[CG]',
				   M   => '[AC]',
				   K   => '[GT]',
				   R   => '[AG]',
				   Y   => '[CT]',
				   B   => '[CGT]',
				   D   => '[AGT]',
				   H   => '[ACT]',
				   V   => '[ACG]',
				   '/' => '|'
);

our %DNA_PROB = (
				  N => 1,
				  W => '$DNA_PROB{A} * 2',
				  S => '$DNA_PROB{G} * 2',
				  M => '$DNA_PROB{A} + $DNA_PROB{G}',
				  K => '$DNA_PROB{A} + $DNA_PROB{G}',
				  R => '$DNA_PROB{A} + $DNA_PROB{G}',
				  Y => '$DNA_PROB{A} + $DNA_PROB{G}',
				  B => '$DNA_PROB{A} + $DNA_PROB{G} * 2',
				  D => '$DNA_PROB{A} * 2 + $DNA_PROB{G}',
				  H => '$DNA_PROB{A} * 2 + $DNA_PROB{G}',
				  V => '$DNA_PROB{A} + $DNA_PROB{G} * 2',
				  A => 0.25,
				  T => 0.25,
				  G => 0.25,
				  C => 0.25
);

my $o_param;

MAIN:
{
	$o_param = New ParamParser( 'GETOPTLONG', \&__Usage, 'fasta=s', 'motif=s', 'output=s', 'prefix=s', 'site_only', 'motif_only', 'stats', 'help', 'debug' );
	pod2usage( -exitval => '2', -verbose => 1 ) if ( $o_param->IsDefined('help') );
	$DEBUG = 1 if ( $o_param->IsDefined('debug') );
	$o_param->AssertFileExists('fasta');
	$o_param->AssertDefined('motif');
	$o_param->SetUnlessDefined( 'output', '.' );
	$o_param->SetUnlessDefined( 'prefix', $$ );

	my $prefix    = $o_param->Get('output') . '/' . $o_param->Get('prefix');
	my $fasta     = $o_param->Get('fasta');
	my $fh_outgff = General::GetStreamOut("${prefix}.motif.gff3.tmp");

	my %h_motif = LoadMotif( $o_param->Get('motif') );
	FindMotif( $fasta, \%h_motif, $fh_outgff );
	
	system("sort -k1,1 -k4,4n ${prefix}.motif.gff3.tmp > ${prefix}.motif.gff3");
	unlink("${prefix}.motif.gff3.tmp");

	exit(0) if ( !$o_param->IsDefined('stats') );

	# stats genome
	foreach my $infoseq (`infoseq --only -name -length -pgc -nohead $fasta 2>/dev/null`)
	{
		if ( $infoseq =~ /^(\S+)\s+(\d+)\s+([0-9.]+)/ )
		{
			$h_motif{data}{size}     += ( $2 + 0 );
			$h_motif{data}{gc}       += ( $3 / 100 );
			$h_motif{data}{scaffold} += 1;
		}
	}
	$h_motif{data}{gc} = $h_motif{data}{gc} / $h_motif{data}{scaffold};

	my $fh_outstat = General::GetStreamOut("${prefix}.motif.csv");
	print $fh_outstat "GC\t$h_motif{data}{gc}\n";
	print $fh_outstat "Genome\t", $h_motif{data}{size}, "\n";
	print $fh_outstat "Scaffold count\t$h_motif{data}{scaffold}\n";
	my @a_header = ( 'motif', 'original', 'type', 'modifpos', 'count', 'probability', 'expected', 'obs/exp ratio' );

	# stats motif
	print $fh_outstat join( "\t", @a_header ), "\n";
	shift @a_header;
	foreach my $key ( sort keys(%h_motif) )
	{
		next if ( $key =~ /result|data/ );

		my @a_motif = split( '', $key );
		my $expected = 0;
		$DNA_PROB{A}                = ( 1 - $h_motif{data}{gc} ) / 2;
		$DNA_PROB{T}                = $DNA_PROB{A};
		$DNA_PROB{G}                = $h_motif{data}{gc} / 2;
		$DNA_PROB{C}                = $DNA_PROB{G};
		$h_motif{$key}{probability} = 1;
		foreach (@a_motif)
		{
			last if ( $_ eq '/' );
			$h_motif{$key}{probability} = $h_motif{$key}{probability} * eval( $DNA_PROB{$_} );
		}
		$h_motif{$key}{expected} = int( $h_motif{$key}{probability} * $h_motif{data}{size} );
		$h_motif{$key}{'obs/exp ratio'} = sprintf( '%2f', ( $h_motif{$key}{count} / $h_motif{$key}{expected} ) );

		print $fh_outstat "$key\t";
		foreach my $field (@a_header)
		{
			print $fh_outstat ( ( exists $h_motif{$key}{$field} ) ? "$h_motif{$key}{$field}\t" : "\t" );
		}
		print $fh_outstat "\n";
	}

}

sub LoadMotif
{
	my $motif = shift;

	my %h_motif;

	if ( -e $motif )
	{
		my $fh_in = General::GetStreamIn($motif);
		while (<$fh_in>)
		{
			chomp;
			next if($_ =~ /^\s*$/);
			my @a_line = split( /\t|,/, $_ );
			my $motif = $a_line[0];
			$motif =~ /([a-z]+)/i;
			my $motiflen = length($1);
			my ($base) = $a_line[1] =~ /m\d([ATGC])/;
			$h_motif{$motif} = {
								 pattern     => dna2regex($motif),
								 'length'    => $motiflen,
								 methylation => $a_line[1],
								 modifpos    => $a_line[2],
								 modified    => $base,
								 motif       => $motif
			};
		}
	}
	else
	{
		my @a_motif = split( /\//, $motif );
		$h_motif{$motif} = {
							 'length' => length( $a_motif[0] ),
							 pattern  => dna2regex($motif)
		};
	}

	return %h_motif;
}

sub dna2regex
{
	my $pattern = shift;

	while ( my ( $base, $regex ) = each(%IUPAC_DNA) )
	{
		$pattern =~ s/$base/$regex/g;
	}

	return $pattern;
}

sub isPalindromic
{
	my $motif = shift;
	return ( $motif eq revcomp($motif) );
}

sub revcomp
{
	my $seq = shift;
	return reverse(comp($seq));
}

sub comp
{
	my $seq = shift;
	$seq =~ tr/ATGCNMRBDKYVH/TACGNKYVHMRBD/;
	$seq =~ tr/atgcnmrbdkyvh/tacgnkyvhmrbd/;
	return $seq;
}

sub FindMotif
{
	my ( $fasta, $rh_motif, $fh_output ) = @_;

	my $fh_fasta = General::GetStreamIn($fasta);
	my $sequence;

	my $firstline = <$fh_fasta>;
	my $seqid = ( $firstline =~ /^>(\S+)/ ) ? $1 : die 'bad fasta file';
	while (<$fh_fasta>)
	{
		my $line = $_;
		chomp $line;

		if ( $line =~ /^>(\S+)/ )
		{
			_FindMotif( $seqid, $sequence, $rh_motif, $fh_output );
			$sequence = '';
			$seqid    = $1;
			next;
		}

		# linearise
		$sequence .= $line;
	}
	_FindMotif( $seqid, $sequence, $rh_motif, $fh_output );
}

sub _FindMotif
{
	my ( $seqid, $sequence, $rh_motif, $fh_output ) = @_;
	return if ( !defined $sequence );

	print "process $seqid\n";

	my %h_data;
	foreach my $motif ( sort( keys(%$rh_motif) ) )
	{
		my $count     = 1;
		my $h_current = $rh_motif->{$motif};

		my $extra_attributs = '';
		foreach my $key ( sort( keys(%$h_current) ) )
		{
			next if ( $key =~ /count|length|pattern|original/ );
			$extra_attributs .= "$key=$h_current->{$key};";
		}
		chop($extra_attributs);
		while ( $sequence =~ m/(?=($h_current->{pattern}))/ig )
		{
			my $real_motif = $1;
			my $start      = $-[0] + 1;
			my $end        = $start + $h_current->{'length'} - 1;
			my $ID         = sprintf( "%06d", $count );

			# motif
			if ( !$o_param->IsDefined('site_only') )
			{
				print $fh_output join( "\t", ( $seqid, 'LIPM', 'DNA_motif', $start, $end, '.', '.', '.', "ID=$seqid.$motif$ID;Name=$motif;$extra_attributs" ) ), "\n";
			}

			# base
			if ( !$o_param->IsDefined('motif_only') )
			{
				if ( exists $h_current->{methylation} && exists $h_current->{modifpos} )
				{
					my $base = $start + $h_current->{modifpos};
					print $fh_output join( "\t", ( $seqid, 'LIPM', $h_current->{methylation}, $base, $base, '.', '+', '.', "ID=$seqid.$motif${ID}.1;Parent=$seqid.$motif$ID;$extra_attributs" ) ), "\n";

					if ( isPalindromic($motif) )
					{
						my $base = $end - $h_current->{modifpos};
						print $fh_output join( "\t", ( $seqid, 'LIPM', $h_current->{methylation}, $base, $base, '.', '-', '.', "ID=$seqid.$motif${ID}.2;Parent=$seqid.$motif$ID;$extra_attributs" ) ), "\n";
					}

				}
			}
			$count++;
		}
		if ( ! isPalindromic($motif) && $motif !~ /\// )
		{
			my $revcomp_motif = revcomp($h_current->{pattern});
			
			while ( $sequence =~ m/(?=($revcomp_motif))/ig )
			{
				my $real_motif = $1;
				my $start      = $-[0] + 1;
				my $end        = $start + $h_current->{'length'} - 1;
				my $ID         = sprintf( "%06d", $count );

				# motif
				if ( !$o_param->IsDefined('site_only') )
				{
					print $fh_output join( "\t", ( $seqid, 'LIPM', 'DNA_motif', $start, $end, '.', '.', '.', "ID=$seqid.$motif$ID;Name=$motif;$extra_attributs" ) ), "\n";
				}

				if ( !$o_param->IsDefined('motif_only') )
				{
					my $base = $end - $h_current->{modifpos};
					print $fh_output join( "\t", ( $seqid, 'LIPM', $h_current->{methylation}, $base, $base, '.', '-', '.', "ID=$seqid.$motif${ID}.2;Parent=$seqid.$motif$ID;$extra_attributs" ) ), "\n";
				}

				$count++;
			}
		}
		
		$rh_motif->{$motif}{count} += ( $motif =~ /\// ) ? ( $count - 1 ) / 2 : $count - 1;
	}
	return;
}

sub Debug
{
	my $msg = shift;
	print STDERR "[DEBUG] $msg\n" if ( $DEBUG == 1 );
	return;
}

sub Info
{
	my $msg = shift;
	print STDOUT "[INFO] $msg\n";
	return;
}

sub Warn
{
	my $msg = shift;
	print STDERR "[WARN] $msg\n";
	return;
}

sub __Usage
{
	pod2usage( -exitval => 'NOEXIT', -verbose => 1 );
}