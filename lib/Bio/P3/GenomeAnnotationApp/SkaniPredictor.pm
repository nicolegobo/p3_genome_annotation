package Bio::P3::GenomeAnnotationApp::SkaniPredictor;

#
# Organism prediction for the Comprehensive Genome Analysis service.
#
# Compares an assembly against a sketch database of BV-BRC Reference and
# Representative bacterial genomes using skani, and resolves the closest match
# to a taxonomy ID.
#
# STATUS: IN DEVELOPMENT. Every user-visible string emitted by this module
# carries the $STATUS marker so it can be found and removed in one pass when
# the feature graduates. See $ON_MISMATCH before changing behaviour.
#
# The policy constants in the CONFIGURATION block below are the only things
# that should need to change as evaluation continues; the logic around them is
# decision-independent.
#

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use IPC::Run;

our $VERSION = '0.1';

#
# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------
#
# Marker on every user-visible string while this feature is in development.
#
our $STATUS = 'DEV';

#
# Database and binary. The database path is Clark's deployed sketch directory;
# this is a placeholder pending a team decision on where the production copy
# lives. All three are overridable by environment variable so the module can be
# exercised outside a deployed backend.
#
our $SKANI_BIN = $ENV{BVBRC_SKANI_BIN} || 'skani';
our $SKANI_DB  = $ENV{BVBRC_SKANI_DB}
    || '/vol/bvbrc/production/application-backend/genome_annotation/skani/current/bvbrc_ref_sketches';
our $SKANI_TAXON_MAP = $ENV{BVBRC_SKANI_TAXON_MAP}
    || '/vol/bvbrc/production/application-backend/genome_annotation/skani/current/genome_taxon_map.tsv';

our $THREADS = 8;

#
# Acceptance rule: three tiers (species / genus / decline).
#
# Adopted 2026-10-14 from the 148-genome positive set and the 120-genome
# negative set (species absent from the sketch by construction); see
# organism-call-tiers.md and final_rule_stats.txt.
#
#   1. species: hits with ANI >= $SPECIES_MIN_ANI and AF_ref >= $SPECIES_MIN_AF_REF.
#      Take the highest ANI. If its confidence interval overlaps the next
#      DISTINCT species in the hit list (that species' ANI_95_percentile >=
#      this hit's ANI_5_percentile), demote to genus.
#   2. genus:   otherwise, hits with ANI >= $GENUS_MIN_ANI and
#      AF_ref >= $GENUS_MIN_AF_REF; report the best hit's genus.
#   3. decline: otherwise, naming the closest reference.
#
# There is deliberately no AF_query gate: once AF_ref is applied it rejected
# nothing further on either set and only cost positives.
#
# 92 is BELOW the ~95-96% species boundary. That is deliberate and provisional:
# the CI demotion, not the ANI floor, is what prevents congener confusion. Do
# not present 92 as a species-boundary claim, and do not loosen or remove the
# demotion without re-running both evaluation sets.
#
our $SPECIES_MIN_ANI    = 92.0;
our $SPECIES_MIN_AF_REF = 70.0;
our $GENUS_MIN_ANI      = 80.0;
our $GENUS_MIN_AF_REF   = 30.0;
our $CI_DEMOTE          = 1;

#
# Floors applied by rank_hits() when the caller passes none. Kept as the
# species-tier values so a bare rank_hits() call can never be looser than the
# rule above.
#
our $MIN_ANI      = $SPECIES_MIN_ANI;
our $MIN_AF_QUERY = 0.0;
our $MIN_AF_REF   = $SPECIES_MIN_AF_REF;

#
# How a genus-tier call is written into scientific_name.
#
# A bare genus ("Escherichia") breaks downstream steps that parse the name:
# AMRFinder's organism lookup needs two tokens (/^(\S+)\s+(\S+)/), and the
# Streptococcus workflow steps test /^Streptococcus\s/. "<Genus> sp." satisfies
# both and is the conventional way to name an organism identified to genus.
#
our $GENUS_NAME_SUFFIX = ' sp.';

#
# How to choose among the hits that pass the thresholds.
#
#   ani    - highest ANI, ties broken on reference genome id (deterministic)
#   afref  - highest reference alignment fraction
#   afmin  - highest min(AF_ref, AF_query)
#
# 'ani' reproduces the prototype's intent while removing its nondeterminism.
# 'afref' scored higher on the evaluation set (94.6% vs 91.9% species-level),
# but the gain rested on a handful of same-genus cases and is not yet trusted
# as a default. Do not change this without re-running the evaluation sweep.
#
our $RANK_BY = 'ani';

#
# What to do when a submitted taxonomy_id disagrees with the prediction.
#
#   die                  - stop the job and report both
#   warn_use_supplied    - warn, annotate with the submitted id
#   warn_use_predicted   - warn, annotate with the prediction    (current default)
#
# warn_use_predicted chosen 2026-10-26 (Nicole). In the expected deployment the
# submitted id comes from the website's minhash suggestion or the user, and
# skani is the confirming check; when they disagree at lineage level, skani's
# alignment-based call is the one we trust. Taxon AND name both switch to the
# prediction, so they cannot come from different organisms.
#
# NOTE for whoever revisits this: on the evaluation set 'die' would stop roughly
# 1 job in 20, and every one of those was a species-complex boundary case rather
# than a submitter error (Enterobacter, Burkholderia, Serratia, Enterococcus,
# E. coli/Shigella). Expect that rate in production.
#
our $ON_MISMATCH = 'warn_use_predicted';

#
# Per-job override from the organism_mismatch_policy parameter (the UI's
# policy selector). $ON_MISMATCH above is the default when the job does not
# set it, e.g. CLI submissions.
#
our %MISMATCH_POLICY = (
    auto     => 'warn_use_predicted',
    supplied => 'warn_use_supplied',
);

#
# ---------------------------------------------------------------------------
# CONSTRUCTION
# ---------------------------------------------------------------------------
#

sub new
{
    my($class, %opts) = @_;
    my $self = {
        db        => $opts{db}        // $SKANI_DB,
        taxon_map => $opts{taxon_map} // $SKANI_TAXON_MAP,
        bin       => $opts{bin}       // $SKANI_BIN,
        threads   => $opts{threads}   // $THREADS,
    };
    return bless $self, $class;
}

#
# True if the sketch database and taxon map are both present. The service
# should treat a false return as "prediction unavailable", not as an error.
#
sub enabled
{
    my($self) = @_;
    return (-d $self->{db} && -s $self->{taxon_map}) ? 1 : 0;
}

#
# ---------------------------------------------------------------------------
# PARSING  (pure functions - unit-testable without skani or the database)
# ---------------------------------------------------------------------------
#

#
# Parse skani tabular output into a list of hashrefs keyed by COLUMN NAME.
#
# Column *positions* are not stable: --detailed and --ci each insert columns,
# which silently shifts the confidence-interval fields. Anything that reads by
# index will read the wrong value the moment the invocation changes. Always go
# through the header.
#
sub parse_output
{
    my($class_or_self, $text) = @_;

    my @lines = split(/\n/, $text // '');
    shift @lines while @lines && $lines[0] !~ /\S/;
    return () unless @lines;

    my $hdr = shift @lines;
    my @cols = split(/\t/, $hdr);

    unless (@cols && $cols[0] eq 'Ref_file')
    {
        die "[$STATUS] skani output is not in the expected format " .
            "(header was '" . substr($hdr, 0, 120) . "')\n";
    }

    my @rows;
    for my $line (@lines)
    {
        next unless $line =~ /\S/;
        my @f = split(/\t/, $line, -1);
        next unless @f == @cols;
        my %r;
        @r{@cols} = @f;
        push @rows, \%r;
    }
    return @rows;
}

#
# Extract the BV-BRC genome id from a reference path: .../fasta/<genome_id>.fna
#
sub genome_id_from_ref
{
    my($class_or_self, $ref_file) = @_;
    return undef unless defined $ref_file;
    # The sketch database uses .fna, but accept the usual FASTA spellings and
    # an optional .gz so this does not quietly return undef if the database is
    # ever rebuilt with different filenames.
    return $1 if $ref_file =~ m{([^/]+?)\.(?:fna|fa|fas|fasta|fsa)(?:\.gz)?$}i;
    return undef;
}

#
# Development-status marker, exposed so callers and tests do not have to reach
# into the package variable.
#
sub status { return $STATUS }

#
# Apply thresholds and pick a single hit, deterministically.
#
# The prototype took the first row skani emitted. skani orders by ANI
# descending but its order among *tied* rows is not stable: the same query, same
# command, different thread count returns a different first row, and at -t 8 it
# varies between runs of the identical command. On the evaluation set 7 of 148
# genomes have a tie at the top ANI and for one of them the reported species
# changes. Hence the explicit sort and the genome-id tiebreak below.
#
sub rank_hits
{
    my($class_or_self, $hits, %opt) = @_;

    my $min_ani = $opt{min_ani}      // $MIN_ANI;
    my $min_afq = $opt{min_af_query} // $MIN_AF_QUERY;
    my $min_afr = $opt{min_af_ref}   // $MIN_AF_REF;
    my $rank_by = $opt{rank_by}      // $RANK_BY;

    my @elig;
    for my $h (@$hits)
    {
        my $ani = $h->{ANI};
        my $afr = $h->{Align_fraction_ref};
        my $afq = $h->{Align_fraction_query};
        next unless defined $ani && defined $afr && defined $afq;
        next unless $ani =~ /^[\d.]+$/ && $afr =~ /^[\d.]+$/ && $afq =~ /^[\d.]+$/;
        next if $ani < $min_ani || $afq < $min_afq || $afr < $min_afr;
        push @elig, $h;
    }
    return () unless @elig;

    my $key = $rank_by eq 'afref' ? sub { $_[0]->{Align_fraction_ref} }
            : $rank_by eq 'afmin' ? sub { my $h = shift;
                                          $h->{Align_fraction_ref} < $h->{Align_fraction_query}
                                              ? $h->{Align_fraction_ref}
                                              : $h->{Align_fraction_query} }
            :                       sub { $_[0]->{ANI} };

    #
    # Descending on the ranking key; ascending on genome id to break ties
    # reproducibly. The tiebreak is arbitrary but must be total and stable.
    #
    my @sorted = sort {
        $key->($b) <=> $key->($a)
            ||
        ($class_or_self->genome_id_from_ref($a->{Ref_file}) // '')
            cmp
        ($class_or_self->genome_id_from_ref($b->{Ref_file}) // '')
    } @elig;

    return @sorted;
}

#
# ---------------------------------------------------------------------------
# TAXONOMY
# ---------------------------------------------------------------------------
#

#
# Look up reference genome ids in the taxon map, in one pass over the file.
# Columns: genome_id, taxon_id, scientific_name, domain, genetic_code
#
# Returns a hashref genome_id => record. Ids absent from the map are absent
# from the result.
#
sub lookup_genome_taxa
{
    my($self, @genome_ids) = @_;

    my %want = map { $_ => 1 } grep { defined } @genome_ids;
    my %found;
    return \%found unless %want;

    open(my $fh, "<", $self->{taxon_map})
        or die "[$STATUS] Cannot open skani taxon map $self->{taxon_map}: $!\n";

    while (my $line = <$fh>)
    {
        chomp $line;
        my @f = split(/\t/, $line, -1);
        next unless @f >= 3 && $want{$f[0]};
        $found{$f[0]} = {
            genome_id       => $f[0],
            taxon_id        => $f[1],
            scientific_name => $f[2],
            domain          => $f[3],
            genetic_code    => $f[4],
        };
        last if keys %found == keys %want;
    }
    close($fh);
    return \%found;
}

#
# Resolve taxon ids to their species and genus ancestors, in one API call.
#
# Species identity MUST come from the current NCBI lineage, not from the
# reference genome's name in the taxon map. Those names are frozen at the time
# each genome was loaded: for 8,603 of the 21,365 sketch references (40%) the
# leading binomial of the stored name is no longer the taxon's species
# (Arcobacter -> Aliarcobacter, Prevotella -> Leyella, "Megasphaera sp." ->
# Megasphaera lornae, ...). Comparing names would treat renamed synonyms as
# distinct species and break the CI-overlap test.
#
# Returns a hashref taxon_id => { species_id, species_name, genus_id,
# genus_name, taxon_name, taxon_rank }, any field of which may be undef (a
# taxon above species rank has no species ancestor). Dies if the API call
# fails, because without lineages the tier rule cannot be applied.
#
sub resolve_lineages
{
    my($self, $api, @taxon_ids) = @_;

    my %uniq = map { $_ => 1 } grep { defined && /^\d+$/ } @taxon_ids;
    my %out;
    return \%out unless %uniq;

    my @t = eval { $api->query("taxonomy",
                               ["in", "taxon_id", "(" . join(",", sort keys %uniq) . ")"],
                               ["select", "taxon_id,taxon_name,taxon_rank,lineage_ids,lineage_names,lineage_ranks"]) };
    die "[$STATUS] taxonomy lineage lookup failed: $@\n" if $@;

    for my $r (@t)
    {
        my %e = (taxon_name => $r->{taxon_name}, taxon_rank => $r->{taxon_rank});
        my($ids, $names, $ranks) = @{$r}{qw(lineage_ids lineage_names lineage_ranks)};
        if (ref $ids eq 'ARRAY' && ref $names eq 'ARRAY' && ref $ranks eq 'ARRAY')
        {
            for my $i (0 .. $#$ids)
            {
                my $rk = $ranks->[$i] // '';
                if ($rk eq 'species') { @e{qw(species_id species_name)} = ($ids->[$i], $names->[$i]) }
                if ($rk eq 'genus')   { @e{qw(genus_id genus_name)}     = ($ids->[$i], $names->[$i]) }
            }
        }
        $out{$r->{taxon_id}} = \%e;
    }
    return \%out;
}

#
# Are two taxonomy ids consistent - is one an ancestor of the other?
#
# A submitter may supply a genus, a species or a strain, so this cannot be an
# equality test. Equality on the *exact* id would disagree on about 42% of jobs,
# because exact-taxon accuracy is only ~58% even when the species is right.
#
# Takes a P3DataAPI handle. Returns 1 (consistent), 0 (divergent), or undef if
# the lineage could not be resolved - callers must distinguish "divergent" from
# "unknown" and must not stop a job on undef.
#
sub taxa_consistent
{
    my($self, $api, $supplied_id, $predicted_id) = @_;

    return undef unless defined $supplied_id && defined $predicted_id;
    return 1 if $supplied_id eq $predicted_id;

    my $lineage_ids = sub {
        my($id) = @_;
        my $md = eval { $api->get_taxon_metadata($id) };
        return undef unless $md && ref $md->{taxon_lineage} eq 'ARRAY';
        return [ map { $_->[1] } @{$md->{taxon_lineage}} ];
    };

    my $sup = $lineage_ids->($supplied_id);
    my $pre = $lineage_ids->($predicted_id);
    return undef unless $sup && $pre;

    return 1 if grep { defined $_ && $_ eq $predicted_id } @$sup;
    return 1 if grep { defined $_ && $_ eq $supplied_id  } @$pre;
    return 0;
}

#
# ---------------------------------------------------------------------------
# PREDICTION
# ---------------------------------------------------------------------------
#

#
# Run skani against the sketch database and return EVERY hit it reports,
# best ANI first. No acceptance thresholds are applied here: the tier rule
# needs the whole list, because the CI-overlap test compares the chosen hit
# against the next distinct species whether or not that species passed.
#
# skani reads gzipped FASTA natively, so the caller may hand this the file it
# already downloaded; there is no need to decompress first.
#
# An empty list means skani found nothing above its own screen, which is a
# legitimate outcome, not a failure.
#
sub search
{
    my($self, $contig_file) = @_;

    #
    # skani exits 0 and writes an empty table when the input is missing or
    # unparseable, which is indistinguishable from a genuine no-hit result.
    # Check the input ourselves so a broken pipeline cannot masquerade as
    # "we could not identify your organism".
    #
    die "[$STATUS] Contig file for organism prediction does not exist: $contig_file\n"
        unless -e $contig_file;
    die "[$STATUS] Contig file for organism prediction is empty: $contig_file\n"
        unless -s $contig_file;
    die "[$STATUS] skani sketch database not found at $self->{db}\n"
        unless -d $self->{db};
    die "[$STATUS] skani taxon map not found at $self->{taxon_map}\n"
        unless -s $self->{taxon_map};

    my @cmd = ($self->{bin}, "search",
               "-d", $self->{db},
               "-q", $contig_file,
               "-t", $self->{threads},
               "--min-af", "15",
               "--detailed",
               "--ci");

    print STDERR "[$STATUS] organism prediction: @cmd\n";

    my($out, $err);
    my $ok = IPC::Run::run(\@cmd, ">", \$out, "2>", \$err);

    die "[$STATUS] skani search failed (exit status indicates error): $err\n"
        unless $ok;

    #
    # Even on exit 0, skani only WARNs when it could not read the query. Treat
    # that as a hard failure rather than as "no hits".
    #
    if (defined $err && $err =~ /is not a valid fasta\/fastq file/)
    {
        die "[$STATUS] skani could not read the contig file $contig_file. " .
            "This is a pipeline error, not an unidentifiable genome.\n";
    }

    print STDERR $err if defined $err && $err =~ /\S/;

    my @hits = $self->parse_output($out);
    return $self->rank_hits(\@hits, min_ani => 0, min_af_query => 0, min_af_ref => 0);
}

#
# Apply the three-tier rule to a list of annotated hits.
#
# Pure function: no skani, no API, no files, so it can be tested directly and
# replayed over the evaluation sets. Each hit must carry the skani columns plus
#   _genome_id, _species_id, _species_name, _genus_id, _genus_name
# as filled in by predict(). Hits may arrive in any order.
#
# Returns
#   { tier => 'species' | 'genus', hit => $chosen, demoted => 0|1,
#     next_species => $hit_or_undef, ci_overlap => 0|1|undef, reason }
#   { tier => 'decline', closest => $hit_or_undef, reason }
#
sub assign_tier
{
    my($class_or_self, $hits, %opt) = @_;

    my $sp_ani  = $opt{species_min_ani}    // $SPECIES_MIN_ANI;
    my $sp_afr  = $opt{species_min_af_ref} // $SPECIES_MIN_AF_REF;
    my $ge_ani  = $opt{genus_min_ani}      // $GENUS_MIN_ANI;
    my $ge_afr  = $opt{genus_min_af_ref}   // $GENUS_MIN_AF_REF;
    my $demote  = $opt{ci_demote}          // $CI_DEMOTE;

    #
    # Every hit, best ANI first with the deterministic tiebreak. The CI-overlap
    # test looks for the next distinct species in THIS list, not only among
    # hits that cleared the species gates: a congener that matches nearly as
    # well is evidence of ambiguity whatever its alignment fraction.
    #
    my @all = $class_or_self->rank_hits($hits, min_ani => 0, min_af_query => 0,
                                        min_af_ref => 0, rank_by => 'ani');

    my @sp = $class_or_self->rank_hits($hits, min_ani => $sp_ani, min_af_query => 0,
                                       min_af_ref => $sp_afr, rank_by => 'ani');
    if (@sp)
    {
        my $top = $sp[0];
        my $top_sp = $top->{_species_id};

        my($next) = grep { defined $_->{_species_id} && defined $top_sp
                               && $_->{_species_id} ne $top_sp } @all;

        my $overlap;
        if ($next && defined $next->{ANI_95_percentile} && defined $top->{ANI_5_percentile}
            && $next->{ANI_95_percentile} =~ /^[\d.]+$/ && $top->{ANI_5_percentile} =~ /^[\d.]+$/)
        {
            $overlap = ($next->{ANI_95_percentile} >= $top->{ANI_5_percentile}) ? 1 : 0;
        }
        elsif ($next)
        {
            #
            # A competing species exists but the CI columns are missing, so we
            # cannot show the call is separable. Treat that as overlap.
            #
            $overlap = 1;
        }
        else
        {
            $overlap = 0;
        }

        #
        # A reference whose taxon has no species ancestor cannot support a
        # species call at all, whatever its scores.
        #
        my $no_species = !defined $top_sp;

        if ((!$demote || !$overlap) && !$no_species)
        {
            return { tier => 'species', hit => $top, demoted => 0,
                     next_species => $next, ci_overlap => $overlap,
                     reason => "best hit passes species gates" .
                               ($next ? " and is separable from the next species" : "") };
        }
        return { tier => 'genus', hit => $top, demoted => 1,
                 next_species => $next, ci_overlap => $overlap,
                 reason => $no_species
                     ? "matched reference has no species-rank ancestor"
                     : "confidence interval overlaps the next species " .
                       "($next->{_species_name}); species-level assignment not supported" };
    }

    my @ge = $class_or_self->rank_hits($hits, min_ani => $ge_ani, min_af_query => 0,
                                       min_af_ref => $ge_afr, rank_by => 'ani');
    if (@ge)
    {
        return { tier => 'genus', hit => $ge[0], demoted => 0,
                 next_species => undef, ci_overlap => undef,
                 reason => "no hit passes species gates (ANI >= $sp_ani, AF_ref >= $sp_afr); " .
                           "best hit passes genus gates (ANI >= $ge_ani, AF_ref >= $ge_afr)" };
    }

    return { tier => 'decline', closest => $all[0],
             reason => "no reference genome met the genus-level thresholds " .
                       "(ANI >= $ge_ani%, reference alignment fraction >= $ge_afr%)" };
}

#
# Full prediction: search, resolve every hit to taxonomy, apply the tier rule.
#
# Takes a P3DataAPI handle, because species and genus identity come from the
# NCBI lineage (see resolve_lineages).
#
# Returns:
#   { status => 'ok', tier => 'species'|'genus', taxon_id, scientific_name,
#     rank, genome_id, reference_name, reference_taxon_id,
#     top_species_id, top_species_name, ani, af_ref,
#     af_query, ani_ci_low, ani_ci_high, demoted, ci_overlap, next_species,
#     tier_reason, n_hits }
#   { status => 'declined', reason, n_hits, closest }
#
# taxon_id and scientific_name are the species (species tier) or the genus
# (genus tier), never the matched reference's own strain taxon.
#
# Dies only on genuine failure. 'declined' is a normal outcome.
#
sub predict
{
    my($self, $contig_file, $api, %opt) = @_;

    my @hits = $self->search($contig_file);

    if (!@hits)
    {
        return {
            status  => 'declined',
            n_hits  => 0,
            reason  => "skani reported no reference genome with any usable alignment",
        };
    }

    for my $h (@hits)
    {
        $h->{_genome_id} = $self->genome_id_from_ref($h->{Ref_file});
        die "[$STATUS] Cannot parse a genome id from the skani reference path " .
            "'$h->{Ref_file}'\n" unless defined $h->{_genome_id};
    }

    my $gtax = $self->lookup_genome_taxa(map { $_->{_genome_id} } @hits);
    for my $h (@hits)
    {
        my $t = $gtax->{$h->{_genome_id}};
        die "[$STATUS] Reference genome $h->{_genome_id} is not present in the " .
            "taxon map $self->{taxon_map}. The sketch database and the taxon map " .
            "are out of step.\n" unless $t;
        $h->{_taxon_id}       = $t->{taxon_id};
        $h->{_reference_name} = $t->{scientific_name};
    }

    my $lin = $self->resolve_lineages($api, map { $_->{_taxon_id} } @hits);
    for my $h (@hits)
    {
        my $l = $lin->{$h->{_taxon_id}} || {};
        @{$h}{qw(_species_id _species_name _genus_id _genus_name)} =
            @{$l}{qw(species_id species_name genus_id genus_name)};
    }

    my $call = $self->assign_tier(\@hits, %opt);

    if ($call->{tier} eq 'decline')
    {
        my $c = $call->{closest};
        return {
            status  => 'declined',
            n_hits  => scalar(@hits),
            reason  => $call->{reason},
            closest => ($c ? $self->_hit_summary($c) : undef),
        };
    }

    my $top = $call->{hit};
    my($taxon_id, $name, $rank);
    if ($call->{tier} eq 'species')
    {
        ($taxon_id, $name, $rank) = ($top->{_species_id}, $top->{_species_name}, 'species');
    }
    else
    {
        unless (defined $top->{_genus_id})
        {
            return {
                status  => 'declined',
                n_hits  => scalar(@hits),
                reason  => "best hit reached the genus tier but its taxon " .
                           "$top->{_taxon_id} has no genus-rank ancestor",
                closest => $self->_hit_summary($top),
            };
        }
        ($taxon_id, $name, $rank) =
            ($top->{_genus_id}, $top->{_genus_name} . $GENUS_NAME_SUFFIX, 'genus');
    }

    return {
        status             => 'ok',
        tier               => $call->{tier},
        taxon_id           => $taxon_id,
        scientific_name    => $name,
        rank               => $rank,
        genome_id          => $top->{_genome_id},
        reference_name     => $top->{_reference_name},
        reference_taxon_id => $top->{_taxon_id},
        top_species_id     => $top->{_species_id},
        top_species_name   => $top->{_species_name},
        ani                => $top->{ANI},
        af_ref             => $top->{Align_fraction_ref},
        af_query           => $top->{Align_fraction_query},
        ani_ci_low         => $top->{ANI_5_percentile},
        ani_ci_high        => $top->{ANI_95_percentile},
        demoted            => $call->{demoted},
        ci_overlap         => $call->{ci_overlap},
        next_species       => ($call->{next_species}
                               ? $self->_hit_summary($call->{next_species}) : undef),
        tier_reason        => $call->{reason},
        n_hits             => scalar(@hits),
    };
}

sub _hit_summary
{
    my($self, $h) = @_;
    return {
        genome_id      => $h->{_genome_id},
        reference_name => $h->{_reference_name},
        species_id     => $h->{_species_id},
        species_name   => $h->{_species_name},
        ani            => $h->{ANI},
        af_ref         => $h->{Align_fraction_ref},
        af_query       => $h->{Align_fraction_query},
        ani_ci_low     => $h->{ANI_5_percentile},
        ani_ci_high    => $h->{ANI_95_percentile},
    };
}

#
# ---------------------------------------------------------------------------
# SUPPLIED-TAXON VALIDATION
# ---------------------------------------------------------------------------
#

#
# Validate a submitted taxonomy_id.
#
# Returns ($clean_id, $name, $reason_if_unusable).
#
# A taxonomy_id that is not an integer, or that does not resolve to an NCBI
# taxon, is NOT a fatal error: the job falls back to automatic prediction and
# the reason is recorded. The submitter is not present to correct it.
#
sub validate_supplied_taxon
{
    my($self, $api, $raw) = @_;

    return (undef, undef, undef) unless defined $raw && $raw =~ /\S/;

    unless ($raw =~ /^\s*(\d+)\s*$/)
    {
        return (undef, undef,
                "provided taxonomy id '$raw' is not an integer and did not map " .
                "to an existing NCBI taxonomy");
    }
    my $id = $1;

    my $md = eval { $api->get_taxon_metadata($id) };
    unless ($md && $md->{taxon_name})
    {
        return (undef, undef,
                "provided taxonomy id '$id' did not map to an existing NCBI taxonomy");
    }

    return ($id, $md->{taxon_name}, undef);
}

#
# ---------------------------------------------------------------------------
# REPORTING
# ---------------------------------------------------------------------------
#

#
# Build the per-job organism-prediction record.
#
# Deliberately records BOTH notions of agreement, because they answer different
# questions and must not be conflated:
#
#   exact_match  - supplied taxon_id is byte-identical to the predicted one.
#                  Disagrees on roughly 42% of jobs even when the species is
#                  right, because the sketch reference is usually a different
#                  strain. Useful as a *statistic*, useless as a *trigger*.
#
#   lineage_consistent - one taxon is an ancestor of the other. Disagrees on
#                  roughly 5%. This is what drives the $ON_MISMATCH action.
#
sub build_report
{
    my($self, %a) = @_;

    my $pred = $a{prediction};
    my %r = (
        status_marker   => $STATUS,
        mode            => $a{mode},                 # supplied | auto | auto_fallback
        supplied_taxon  => $a{supplied_taxon},
        supplied_name   => $a{supplied_name},
        supplied_reason => $a{supplied_reason},      # why a supplied id was rejected
        rule            => {
            species_min_ani    => $SPECIES_MIN_ANI,
            species_min_af_ref => $SPECIES_MIN_AF_REF,
            genus_min_ani      => $GENUS_MIN_ANI,
            genus_min_af_ref   => $GENUS_MIN_AF_REF,
            ci_demote          => $CI_DEMOTE,
        },
        prediction_status => ($pred ? $pred->{status} : 'not_run'),
    );

    if ($pred && $pred->{status} eq 'ok')
    {
        $r{tier}               = $pred->{tier};
        $r{tier_reason}        = $pred->{tier_reason};
        $r{demoted}            = $pred->{demoted};
        $r{ci_overlap}         = $pred->{ci_overlap};
        $r{predicted_taxon}    = $pred->{taxon_id};
        $r{predicted_name}     = $pred->{scientific_name};
        $r{predicted_rank}     = $pred->{rank};
        $r{reference_genome}   = $pred->{genome_id};
        $r{reference_name}     = $pred->{reference_name};
        $r{reference_taxon}    = $pred->{reference_taxon_id};
        $r{top_species_id}     = $pred->{top_species_id};      # best hit's species, even when
        $r{top_species_name}   = $pred->{top_species_name};    # the call was demoted to genus
        $r{ani}                = $pred->{ani};
        $r{af_ref}             = $pred->{af_ref};
        $r{af_query}           = $pred->{af_query};
        $r{ani_ci_low}         = $pred->{ani_ci_low};
        $r{ani_ci_high}        = $pred->{ani_ci_high};
        $r{n_hits}             = $pred->{n_hits};
        if (my $n = $pred->{next_species})
        {
            #
            # The margin to the next distinct species is the only signal that
            # separates a congener confusion from a clean call; the scores
            # themselves do not.
            #
            $r{next_species_name}    = $n->{species_name};
            $r{next_species_genome}  = $n->{genome_id};
            $r{next_species_ani}     = $n->{ani};
            $r{next_species_ci_high} = $n->{ani_ci_high};
            $r{ani_gap} = sprintf("%.2f", $pred->{ani} - $n->{ani});
        }
    }
    elsif ($pred)
    {
        $r{decline_reason} = $pred->{reason};
        $r{n_hits}         = $pred->{n_hits};
        if (my $c = $pred->{closest})
        {
            $r{closest_genome} = $c->{genome_id};
            $r{closest_name}   = $c->{reference_name};
            $r{closest_ani}    = $c->{ani};
            $r{closest_af_ref} = $c->{af_ref};
        }
    }

    #
    # What the job actually proceeds with. Differs from predicted_taxon when a
    # supplied id was kept.
    #
    $r{final_taxon} = $a{final_taxon};
    $r{final_name}  = $a{final_name};

    if (my $m = $a{minhash})
    {
        @r{keys %$m} = values %$m;
    }

    $r{exact_match}        = $a{exact_match};
    $r{lineage_consistent} = $a{lineage_consistent};
    $r{mismatch_policy}    = $a{mismatch_policy};    # auto | supplied | undef (default)
    $r{on_mismatch}        = $a{on_mismatch} // $ON_MISMATCH;
    $r{action}             = $a{action};

    return \%r;
}

#
# Summarize a Similar Genome Finder (MinHash) search recorded by the submission
# form, and compare its top species with skani's, for the per-job report.
#
# The record comes from the browser, so it is reported, never trusted: it does
# not influence which organism the job uses. Returns undef when no search was
# run.
#
# Keys (all prefixed minhash_ so they cannot collide with skani's):
#   minhash_hits            every hit shown, best first, as submitted
#   minhash_top_*           the best hit
#   minhash_selected_*      the result row the submitter selected, if any
#   minhash_taxon_source    'similar_genome_finder' if the submitted
#                           taxonomy_id is the one filled in from the finder,
#                           'manual' if the submitter changed it afterwards
#   minhash_skani_species_match       top MinHash species == skani's best-hit
#                                     species (1/0, undef if either is missing)
#   minhash_skani_lineage_consistent  top MinHash species is consistent with
#                                     skani's call (see taxa_consistent)
#
sub summarize_minhash
{
    my($self, $api, $sgf, $supplied_id, $prediction) = @_;

    return undef unless ref $sgf eq 'HASH';

    my @hits = grep { ref $_ eq 'HASH' } (ref $sgf->{hits} eq 'ARRAY' ? @{$sgf->{hits}} : ());
    my %fields = map { $_ => 1 } qw(genome_id genome_name species_taxon_id species_name
                                    distance pvalue kmer_counts est_ani);
    @hits = map { my $h = $_; +{ map { $_ => $h->{$_} } grep { $fields{$_} } keys %$h } } @hits;

    my %r = (
        minhash_ran               => 1,
        minhash_input_file        => $sgf->{input_file},
        minhash_include_bacterial => $sgf->{include_bacterial},
        minhash_include_viral     => $sgf->{include_viral},
        minhash_max_hits          => $sgf->{max_hits},
        minhash_n_hits            => scalar(@hits),
        minhash_hits              => \@hits,
        minhash_selected_genome   => $sgf->{selected_genome_id},
        minhash_selected_taxon    => $sgf->{selected_taxon_id},
    );

    if (defined $sgf->{selected_taxon_id} && defined $supplied_id)
    {
        $r{minhash_taxon_source} = $sgf->{selected_taxon_id} eq $supplied_id
            ? 'similar_genome_finder' : 'manual';
    }

    my $top = $hits[0] or return \%r;
    $r{minhash_top_genome}       = $top->{genome_id};
    $r{minhash_top_species_id}   = $top->{species_taxon_id};
    $r{minhash_top_species_name} = $top->{species_name};
    $r{minhash_top_distance}     = $top->{distance};
    $r{minhash_top_est_ani}      = $top->{est_ani};
    $r{minhash_top_kmer_counts}  = $top->{kmer_counts};

    if ($prediction && $prediction->{status} eq 'ok' && defined $top->{species_taxon_id})
    {
        $r{minhash_skani_species_match} = defined $prediction->{top_species_id}
            ? ($top->{species_taxon_id} eq $prediction->{top_species_id} ? 1 : 0)
            : undef;
        $r{minhash_skani_lineage_consistent} =
            $self->taxa_consistent($api, $top->{species_taxon_id}, $prediction->{taxon_id});
    }
    return \%r;
}

#
# Emit the record as a single greppable, parseable line so jobs can be
# aggregated into a mismatch report without scraping prose.
#
sub emit_report
{
    my($self, $report) = @_;

    my $json = eval { require JSON::XS;
                      JSON::XS->new->canonical(1)->encode($report) };
    if (!defined $json)
    {
        $json = join(" ", map { "$_=" . (defined $report->{$_} ? $report->{$_} : '') }
                     sort keys %$report);
    }
    print STDERR "[$STATUS][ORGANISM_PREDICTION] $json\n";
    return $json;
}

#
# ---------------------------------------------------------------------------
# USER-FACING MESSAGES
# ---------------------------------------------------------------------------
#

#
# Swap the organism part of a submitted genome name, keeping the submitter's
# label.
#
# The website builds scientific_name as "<taxon name> <My Label>", with
# ( ) | / : stripped from the taxon name. When that prefix is found, the
# label after it is kept: "Escherichia coli isolate-42" becomes
# "Klebsiella pneumoniae isolate-42". Any strain text in the supplied
# taxon name is dropped with it, since the prediction never names a strain.
# A name without that prefix (e.g. a CLI submission) is replaced whole.
#
sub replace_organism_name
{
    my($self, $name, $supplied_taxon_name, $new_organism) = @_;

    return $new_organism unless defined $name && $name =~ /\S/ && defined $supplied_taxon_name;

    (my $cleaned = $supplied_taxon_name) =~ s/[()|\/:]//g;
    for my $prefix ($cleaned, $supplied_taxon_name)
    {
        next unless $prefix =~ /\S/;
        return $new_organism . ($1 // '') if $name =~ /^\s*\Q$prefix\E(\s+\S.*)?$/;
    }
    return $new_organism;
}

#
# Map the job's organism_mismatch_policy to an $ON_MISMATCH-style action.
# An absent or empty policy gives the module default; an unknown one dies, so
# a typo cannot silently pick a behaviour.
#
sub resolve_mismatch_policy
{
    my($self, $policy) = @_;

    return $ON_MISMATCH if !defined $policy || $policy eq '';
    return $MISMATCH_POLICY{$policy} if exists $MISMATCH_POLICY{$policy};
    die "Invalid organism_mismatch_policy '$policy'; expected one of: " .
        join(", ", sort keys %MISMATCH_POLICY) . "\n";
}

#
# Case (b): the submitted taxonomy_id disagrees with the prediction.
#
# Behaviour is governed by the on_mismatch argument (from the job's
# organism_mismatch_policy), else $ON_MISMATCH. Returns ($taxonomy_id, $scientific_name)
# for the caller to proceed with, or dies. $scientific_name is undef when the
# caller's existing name should be kept (the supplied id was kept). When the
# prediction is used, BOTH come from it: switching the id but keeping the
# supplied name would label one organism with another's name, and AMRFinder and
# the Streptococcus workflow steps key off the name, not the id.
#
sub handle_mismatch
{
    my($self, %a) = @_;
    my($supplied, $supplied_name, $pred) = @a{qw(supplied_taxon_id supplied_name prediction)};
    my $on_mismatch = $a{on_mismatch} // $ON_MISMATCH;

    my $detail = sprintf(
        "    supplied   : %-10s %s\n" .
        "    predicted  : %-10s %s   (%s level)\n" .
        "    closest ref: %-10s %s   ANI %.2f%%, AF_ref %.1f%%\n",
        $supplied, ($supplied_name // '(name not resolved)'),
        $pred->{taxon_id}, $pred->{scientific_name}, ($pred->{tier} // '?'),
        $pred->{genome_id}, ($pred->{reference_name} // ''), $pred->{ani}, $pred->{af_ref});

    my $preamble =
        "[$STATUS - organism prediction is under active development; this check\n" .
        "is informational and is not yet authoritative.]\n\n" .
        "The supplied taxonomy_id does not match the closest genome in the\n" .
        "BV-BRC reference set:\n\n" . $detail . "\n" .
        "Closely related species frequently cannot be distinguished by this\n" .
        "method; a mismatch here does not necessarily indicate an error.\n";

    if ($on_mismatch eq 'warn_use_supplied')
    {
        print STDERR $preamble . "Proceeding with the SUPPLIED taxonomy_id.\n";
        return ($supplied, undef);
    }
    elsif ($on_mismatch eq 'warn_use_predicted')
    {
        print STDERR $preamble . "Proceeding with the PREDICTED taxonomy_id and name " .
            "($pred->{taxon_id}, $pred->{scientific_name}).\n";
        return ($pred->{taxon_id}, $pred->{scientific_name});
    }

    die $preamble .
        "\nStopping the job. To proceed, either correct the taxonomy_id or\n" .
        "resubmit without one to accept the predicted organism.\n";
}

1;
