#
# Compute genome MLST using https://github.com/tseemann/mlst
#

use Bio::KBase::AppService::AppConfig qw(application_backend_dir);
use P3DataAPI;
use strict;
use Getopt::Long::Descriptive;  
use File::Temp;
use File::Copy;
use GenomeTypeObject;
use Cwd;
use Data::Dumper;
use Time::HiRes qw(gettimeofday);
use File::Slurp;
use JSON::XS;
use Data::Dumper;
use Text::CSV_XS qw(csv);
use IPC::Run qw(run);

my ($opt, $usage) = describe_options("%c %o [< in] [> out]",
                     ["in|i=s", "Input GTO"],
                     ["out|o=s", "Output GTO"],
                     ["parallel=i", "Number of threads to use", { default => 1 }],
                     ["dry_run|n", "Dry run - print commands without executing"],
                     ["help|h", "Print this help message"]);
                     
print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if (@ARGV != 0);

chomp(my $hostname = `hostname -f`);

# helper sub to run or dry-run a command
sub run_cmd {
    my ($cmd_ref, $label) = @_;
    print STDERR "Invoke $label command: @$cmd_ref\n";
    if ($opt->dry_run) {
        print STDERR "DRY RUN - would execute: @$cmd_ref\n";
        return 1;
    }
    my($stdout, $stderr);
    my $ok = run($cmd_ref,
                 ">", \$stdout,
                 "2>", \$stderr);
    $ok or die "Error $? running $label: @$cmd_ref\n";
    return $ok;
}

my $gto_in;

if ($opt->in)
{
    $gto_in = GenomeTypeObject->new({file => $opt->in});
    $gto_in or die "Error reading GTO from " . $opt->in . "\n";
}
else
{
    $gto_in = GenomeTypeObject->new({file => \*STDIN});
    $gto_in or die "Error reading GTO from standard input\n";
}

my $in_file = File::Temp->new();
$gto_in->write_contigs_to_file($in_file);

# evaluate the mlst scheme to use based on the lineage of the genome
my %schema_map = (
    acinetobacter_baumannii => 470,
    bacillus_anthracis => 1392,
    bordetella_pertussis => 520,
    brucella_abortus => 235,
    brucella_canis => 36855,
    brucella_ceti => 120577,
    brucella_inopinata => 1218315,
    brucella_melitensis => 29459,
    brucella_microti => 444163,
    brucella_neotomae => 29460,
    brucella_ovis => 236,
    brucella_pinnipedialis => 120576,
    brucella_suis => 29461,
    burkholderia_mallei_fli => 13373,
    burkholderia_mallei_rki => 13373,
    burkholderia_pseudomallei => 111527,
    campylobacter_coli => 197,
    campylobacter_jejuni => 195,
    clostridioides_difficile => 1496,
    clostridium_perfringens => 1502,
    corynebacterium_diphtheriae => 1717,
    corynebacterium_pseudotuberculosis => 1719,
    cronobacter_malonaticus => 413503,
    cronobacter_sakazakii => 28141,
    enterococcus_faecalis => 1351,
    enterococcus_faecium => 1352,
    escherichia_coli => 562,
    francisella_tularensis => 263,
    klebsiella_grimontii => 2058152,
    klebsiella_michiganensis => 1134687,
    klebsiella_oxytoca => 571,
    klebsiella_pasteurii => 2587529,
    klebsiella_pneumoniae => 573,
    klebsiella_quasipneumoniae => 1463165,
    klebsiella_variicola => 244366,
    legionella_pneumophila => 446,
    listeria_monocytogenes => 1639,
    mycobacterium_africanum => 33894,
    mycobacterium_bovis => 1765,
    mycobacterium_canettii => 78331,
    mycobacterium_tuberculosis => 77643,
    mycobacteroides_abscessus => 36809,
    mycoplasma_gallisepticum => 2096,
    paenibacillus_larvae => 1464,
    proteus_mirabilis => 584,
    providencia_stuartii => 588,
    pseudomonas_aeruginosa => 287,
    Salmonella_enterica => 28901,
    serratia_marcescens => 615,
    staphylococcus_argenteus => 985002,
    staphylococcus_aureus => 1280,
    staphylococcus_capitis => 29388,
    streptococcus_pyogenes => 1314,
    yersinia_enterocolitica => 630
);


# grab all lineage 
my $api = P3DataAPI->new();
my $taxon_id = $gto_in->{ncbi_taxonomy_id};
my @res = $api->query("taxonomy", ['eq', 'taxon_id', $taxon_id], ['select', 'taxon_name', 'lineage_ids', 'lineage_names']);

# Create reverse hash (taxon_id => schema dir name)
my %taxon_lookup = reverse %schema_map;
my $dir_name;
foreach my $lineage_id (@{$res[0]->{lineage_ids}}) {
    if (exists $taxon_lookup{$lineage_id}) {
        $dir_name = $taxon_lookup{$lineage_id};
    }
}

print STDERR "dir_name: $dir_name\n";

# based on the schema name, determine there is a cgmlst scheme to use
if ($dir_name) {
    # my $tmp_dir = "tmp_dir_manual";
    my $tmp_dir = File::Temp->newdir(CLEANUP => 0);
    my $clean_fasta_dir = "$tmp_dir/clean_fastas";
    my $allele_call_out = "$tmp_dir/new_genomes_allele_call/";
    my $schema_path = application_backend_dir . "/CoreGenomeMLST/chewbbaca_schemas/" . $dir_name;
    
    # mkdir $tmp_dir or die "Cannot create tmp directory: $!\n";
    mkdir $clean_fasta_dir or die "Cannot create clean fastas directory: $!\n";
    
    # rename to be certain it ends with 'fasta' for chewBBACA
    run_cmd(["cp", "$in_file", "$clean_fasta_dir/input.fasta"], "Add extension");

    # new genome allele call
    my @allele_call_cmd = (
        "chewBBACA.py", "AlleleCall",
        "--input-files", $clean_fasta_dir,
        "--schema-directory", $schema_path,
        "--output-directory", $allele_call_out,
        "--cpu", "4",
        "--output-unclassified",
        "--output-missing",
        "--output-novel",
        "--no-inferred"
    );
    run_cmd(\@allele_call_cmd, "Allele Call");

    # check 95% of loci have an exact allele match
    my $new_allele_call = "$tmp_dir/new_genomes_allele_call/results_alleles.tsv";
    if (!$opt->dry_run) {
        open(my $check_fh, '<', $new_allele_call) or die "Cannot open $new_allele_call for validation: $!";
        my $check_header = <$check_fh>;  # skip header
        my $check_data = <$check_fh>;    # grab second line
        close($check_fh);
        chomp $check_data;
        my @check_fields = split(/\t/, $check_data);
        shift @check_fields;  # remove first column (filename)

        my $total_loci = scalar @check_fields;
        my $exact_matches = 0;
        for my $val (@check_fields) {
            $exact_matches++ if $val =~ /^\d+$/;
        }

        my $exact_pct = ($total_loci > 0) ? ($exact_matches / $total_loci) * 100 : 0;
        print STDERR sprintf("Allele call check: %d / %d loci (%.1f%%) have exact allele matches\n",
                             $exact_matches, $total_loci, $exact_pct);

        if ($exact_pct < 95) {
            die sprintf("Only %.1f%% of loci have exact allele matches (threshold: 95%%). Aborting.\n", $exact_pct);
        }
    }

    # Join new allele call with master allele call using chewBBACA JoinProfiles
    my $master_table_dir = "/home/nbowers/bvbrc-dev/dev_container/cgmlst_for_all/rerunning_joins_11_25_2025/raw_genome_ids/";
    my $master_table = $master_table_dir . $dir_name . "_11_25_2025_joined.tsv";
    my $master_joined = "$tmp_dir/master_joined.tsv";
    my $copy_of_master = "$tmp_dir/master_copy.tsv";
    my $clean_master_joined = "$tmp_dir/master_joined_clean.tsv";

    run_cmd([
        "cp",  $master_table, $copy_of_master 
    ], "Copy Master Table");

    run_cmd([
        "chewBBACA.py", "JoinProfiles",
        "--profiles", $copy_of_master, $new_allele_call,
        "--output-file", $master_joined
    ], "Join Profiles");

    # cleaning for pHierCC
    run_cmd([
        "python3", "/home/nbowers/bvbrc-dev/dev_container/modules/bvbrc_CoreGenomeMLST/service-scripts/core-genome-mlst-utils.py",
        "clean-allelic-profile", $master_joined
    ], "Clean Allele Call");

    my $precomputed_clusters_dir = application_backend_dir . "/CoreGenomeMLST/precomputed_clusters/";
    my $precomputed_clusters_path = $precomputed_clusters_dir . lc($dir_name) . ".cgMLSTv1.npz";
    my $local_clusters_path = "$tmp_dir/precomputed_clusters.npz";
    my $heircc_out = "$tmp_dir/cluster";

    # copy over npz file
        
    run_cmd([
        "cp",  $precomputed_clusters_path, $local_clusters_path 
    ], "Copy Master Table");

    # clustering with pHierCC
    run_cmd([
        "python3", "/home/nbowers/bvbrc-dev/dev_container/cgmlst_for_all/dev_clustering/testing_heirCC/testing_pHierCC/git_repo/pHierCC/pHierCC.py",
        "--profile", $clean_master_joined,
        "--output", $heircc_out,
        "--append", $precomputed_clusters_path
    ], "Cluster");

    my $heircc_tsv = "$tmp_dir/cluster.HierCC.gz";
    run_cmd(["gunzip", $heircc_tsv], "gunzip");

    # get genome ID
    my $genome_id = $gto_in->{id};

    # Step 1: Parse 2nd line of new_allele_call into a string
    my $allele_call_string;
    if ($opt->dry_run) {
        print STDERR "DRY RUN - would parse allele calls from $new_allele_call\n";
        $allele_call_string = "DRY_RUN_ALLELE_CALLS";
    } else {
        open(my $allele_fh, '<', $new_allele_call) or die "Cannot open $new_allele_call: $!";
        my $header_line = <$allele_fh>;  # skip header
        my $data_line   = <$allele_fh>;  # grab second line
        close($allele_fh);
        chomp $data_line;
        $data_line =~ s/\t/,/g;  # replace tabs with commas
        $data_line =~ s/^[^,]*,//;  # remove the first column (filename)
        $allele_call_string = $data_line;
        print STDERR "Allele call string: $allele_call_string\n";
    }

    # Step 2: Parse heircc_tsv (unzipped) into key value pairs
    # find the row matching the genome_id and zip with header
    my $heircc_unzipped = "$tmp_dir/cluster.HierCC";
    my %cluster_kv;
    if ($opt->dry_run) {
        print STDERR "DRY RUN $heircc_unzipped\n";
        %cluster_kv = (HC1 => "DRY_RUN", HC5 => "DRY_RUN");
    } else {
        open(my $heircc_fh, '<', $heircc_unzipped) or die "Cannot open $heircc_unzipped: $!";
        my $heircc_header = <$heircc_fh>;
        chomp $heircc_header;
        my @heircc_keys = split(/\t/, $heircc_header);

        # find the row where the first column matches the genome_id
        my $heircc_data;
        my $line_num = 0;
        while (my $line = <$heircc_fh>) {
            chomp $line;
            $line_num++;
            my @cols = split(/\t/, $line);
            if ($cols[0] eq 'input') {  # use eq for exact match instead of regex
                print STDERR "Found input row at line $line_num\n";
                $heircc_data = $line;
                last;
            }
        }
        close($heircc_fh);

        die "Could not find $genome_id row in $heircc_unzipped\n" unless $heircc_data;

        my @heircc_values = split(/\t/, $heircc_data);

        shift @heircc_values;  # remove the first column (genome ID)
        shift @heircc_keys;    # remove the first column (genome ID)
        @cluster_kv{@heircc_keys} = @heircc_values;

        print STDERR "Cluster key-value pairs:\n";
        foreach my $key (sort keys %cluster_kv) {
            print STDERR "  $key => $cluster_kv{$key}\n";
        }
    }

    # Build the event and typing objects
    my $event = {
        tool_name      => "p3x-compute-cgmlst",
        parameters     => [map { "$_" } @allele_call_cmd],
        execution_time => scalar gettimeofday,
        hostname       => $hostname,
    };
    my $event_id = $gto_in->add_analysis_event($event);

    my $typing = {
        typing_method           => "sequence_typing",
        allele_calls            => $allele_call_string,
        schema_name             => $dir_name,
        cluster_key_value_pairs => \%cluster_kv,
        event_id                => $event_id,
    };
    push(@{$gto_in->{typing}}, $typing);

} else {
    print STDERR "Schema does not exist for this speices\n";
}

if ($opt->out)
{
    $gto_in->destroy_to_file($opt->out);
}
else
{
    $gto_in->destroy_to_file(\*STDOUT);
}