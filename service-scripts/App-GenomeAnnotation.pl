#
# The Genome Annotation application.
# temp copy that rips out the checkm/eval code
#

use Bio::KBase::AppService::AppScript;
use Bio::P3::GenomeAnnotationApp::GenomeAnnotationCore;
use Bio::P3::GenomeAnnotationApp::SkaniPredictor;
use Bio::KBase::AppService::AppConfig qw(data_api_url db_host db_user db_pass db_name seedtk);
use Bio::KBase::AppService::FastaParser 'parse_fasta';
use Bio::KBase::AppService::LongestCommonSubstring qw(BuildString BuildTree LongestCommonSubstring);
use IPC::Run;
use SolrAPI;
use DBI;

use strict;
use Data::Dumper;
use gjoseqlib;
use File::Basename;
use File::Temp;
use LWP::UserAgent;
use JSON::XS;
use IPC::Run;
use IO::File;
use Module::Metadata;
use GenomeTypeObject;

my $script = Bio::KBase::AppService::AppScript->new(\&process_genome, \&preflight);

my $rc = $script->run(\@ARGV);

exit $rc;

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;

    #
    # Ensure the contigs are valid, and look up their size.
    #

    my $ctg = $params->{contigs};
    $ctg or die "Contigs must be specified\n";

    my $res = $app->workspace->stat($ctg);
    $res->size > 0 or die "Contigs not found\n";

    #
    # Size estimate based on conservative 500 bytes/second aggregate
    # compute rate for contig size, with a minimum allocated
    # time of 60 minutes (to account for non-annotation portions).
    #
    my $time = $res->size / 500;
    $time = 7200 if $time < 7200;

    my $ram = "16G";
    if ($res->size > 10_000_000)
    {
	$ram = "128G";
    }

    #
    # Request 8 cpus for some of the fatter bits of the compute.
    #
    return {
	cpu => 8,
	memory => $ram,
	runtime => int($time),
	storage => 10 * $res->size,
    };
}

sub process_genome
{
    my($app, $app_def, $raw_params, $params) = @_;

    print "Proc genome ", Dumper($app_def, $raw_params, $params);

    my $json = JSON::XS->new->pretty(1)->canonical(1);

    #
    # Do some sanity checking on params.
    #
    # Both recipe and workflow may not be specified.
    #
    if ($params->{workflow} && $params->{recipe})
    {
	die "Both a workflow document and a recipe may not be supplied to an annotation request";
    }

    my $core = Bio::P3::GenomeAnnotationApp::GenomeAnnotationCore->new(app => $app,
								       app_def => $app_def,
								       params => $params);

    if (exists($raw_params->{tax_id}) && !exists($params->{taxonomy_id}))
    {
	print STDERR "Fixup incorrect taxid in parameters\n";
	$params->{taxonomy_id} = $raw_params->{tax_id};
    }

    #
    # Validate any submitted taxonomy_id, but do NOT fail on a bad one.
    #
    # A non-integer value, or an integer that does not resolve to an NCBI
    # taxon, falls back to automatic prediction with the reason recorded. The
    # submitter is not present to correct it, so failing the job would just
    # throw the work away.
    #
    # (The previous code assigned $1 unconditionally, which leaked whatever the
    # last successful match in scope had captured whenever taxonomy_id was
    # absent -- precisely the case prediction now has to handle.)
    #
    my $user_id = $core->user_id;
    my $api = P3DataAPI->new();

    my $predictor = Bio::P3::GenomeAnnotationApp::SkaniPredictor->new();

    #
    # What to do if the supplied taxonomy_id disagrees with the prediction:
    # "auto" uses the prediction, "supplied" keeps the submitted id. Checked
    # up front so a bad value fails before any work is done.
    #
    my $mismatch_policy = $params->{organism_mismatch_policy};
    my $on_mismatch = $predictor->resolve_mismatch_policy($mismatch_policy);

    my $supplied_raw = $params->{taxonomy_id};
    my($supplied_id, $supplied_name, $supplied_reason) =
        $predictor->validate_supplied_taxon($api, $supplied_raw);

    if (defined $supplied_reason && $on_mismatch eq 'warn_use_supplied')
    {
	#
	# The submitter asked us to keep their taxon, so substituting a predicted
	# one would override that choice. Stop here, before the contigs are
	# downloaded, with the report line written first.
	#
	my $fatal = "[" . $predictor->status . "] The job was stopped because the " .
	    "$supplied_reason. organism_mismatch_policy is 'supplied', so a predicted " .
	    "organism will not be substituted. Resubmit with a valid NCBI taxonomy_id, " .
	    "or set organism_mismatch_policy to 'auto' to use the predicted organism.\n";
	$predictor->emit_report($predictor->build_report(
	    mode            => 'supplied_invalid',
	    minhash         => scalar(eval { $predictor->summarize_minhash(
				   $api, $params->{similar_genome_finder}) }),
	    supplied_taxon  => $supplied_raw,
	    supplied_reason => $supplied_reason,
	    mismatch_policy => $mismatch_policy,
	    on_mismatch     => $on_mismatch,
	    action          => 'failed_invalid_supplied_taxon'));
	die $fatal;
    }
    elsif (defined $supplied_reason)
    {
	print STDERR "[" . $predictor->status . "] $supplied_reason; " .
	    "falling back to automatic organism prediction\n";
	delete $params->{taxonomy_id};
    }
    else
    {
	$params->{taxonomy_id} = $supplied_id if defined $supplied_id;

	#
	# scientific_name is optional (CLI jobs may send only a taxonomy_id).
	# Fill it from the taxon so downstream steps that key off the name --
	# AMRFinder's organism flag, the Streptococcus workflow steps -- see
	# the same organism as the id.
	#
	if (defined $supplied_id && !($params->{scientific_name} && $params->{scientific_name} =~ /\S/))
	{
	    $params->{scientific_name} = $supplied_name;
	}
    }

    my $prediction_mode = defined $supplied_reason ? 'auto_fallback'
                        : defined $supplied_id     ? 'supplied'
                        :                            'auto';

    #
    # Stage the contigs once, here. Organism prediction needs them before the
    # genome object can be created, and the annotation pipeline below reuses
    # this same copy rather than downloading a second one.
    #
    my $ws = $app->workspace();
    my($input_path) = $params->{contigs};

    my $temp = File::Temp->new();
    $ws->copy_files_to_handles(1, $core->token, [[$input_path, $temp]]);
    close($temp);

    #
    # Organism prediction (IN DEVELOPMENT).
    #
    #   a. taxonomy_id supplied and consistent with the prediction -> confirm
    #   b. taxonomy_id supplied and divergent                      -> organism_mismatch_policy
    #   c. no usable taxonomy_id                                   -> use prediction
    #
    # Runs on every job, including those that supply a correct taxonomy_id,
    # because cases (a) and (b) cannot exist otherwise. Cost is ~15s wall and
    # ~1.6GB RSS per job, of which roughly 95% is container startup and the
    # marker-index load; the comparison itself is ~0.5s.
    #
    # The prediction is reported at species or genus rank (see the tier rule in
    # SkaniPredictor), never at the matched reference's strain.
    #
    # Any fatal outcome is held in $fatal until the organism report below has
    # been written, so every job emits exactly one report line, including the
    # ones that stop here.
    #
    my $prediction;
    my($exact_match, $lineage_consistent, $action, $fatal, $taxon_overridden);

    if (!$predictor->enabled)
    {
	$action = 'prediction_unavailable';
	print STDERR "Organism prediction unavailable (sketch database not found)\n";
	if (!$params->{taxonomy_id})
	{
	    $fatal = "Cannot determine the organism for these contigs: no usable " .
		"taxonomy_id was supplied and organism prediction is not configured " .
		"on this deployment.\n";
	}
    }
    else
    {
	$prediction = $predictor->predict("$temp", $api);

	if ($prediction->{status} eq 'declined')
	{
	    if (!$params->{taxonomy_id})
	    {
		#
		# Nothing to fall back to: prediction found no reference above the
		# genus-level thresholds and no usable taxonomy_id was supplied.
		#
		$action = 'failed_unidentifiable';
		my $c = $prediction->{closest};
		$fatal = "[" . $predictor->status . "] Could not determine the organism " .
		    "from these contigs: $prediction->{reason}. The assembly did not " .
		    "match any genome in the BV-BRC reference set closely enough to " .
		    "identify it, even at genus level." .
		    ($c ? sprintf(" Closest reference: %s (%s), ANI %.2f%%, reference " .
				  "alignment fraction %.1f%%.",
				  $c->{reference_name}, $c->{genome_id}, $c->{ani}, $c->{af_ref})
			: "") .
		    " Resubmit with a taxonomy_id to annotate it anyway.\n";
	    }
	    else
	    {
		$action = 'declined_kept_supplied';
		print STDERR "[" . $predictor->status . "] organism prediction declined " .
		    "($prediction->{reason}); proceeding with taxonomy_id " .
		    "$params->{taxonomy_id}\n";
	    }
	}
	else
	{
	    #
	    # Record BOTH notions of agreement. exact_match is the statistic for
	    # the mismatch report; lineage consistency is what drives the action.
	    #
	    if (defined $params->{taxonomy_id})
	    {
		$exact_match = ($params->{taxonomy_id} eq $prediction->{taxon_id}) ? 1 : 0;
		$lineage_consistent = $predictor->taxa_consistent(
		    $api, $params->{taxonomy_id}, $prediction->{taxon_id});
	    }

	    if (!defined $params->{taxonomy_id})
	    {
		#
		# Case (c). Taxon and name both come from the prediction, which
		# resolved them from the same NCBI lineage, so they cannot disagree.
		# The submitter's label follows the organism, as the website does
		# when a taxonomy is chosen.
		#
		$action = "used_prediction_$prediction->{tier}";
		$params->{taxonomy_id}     = $prediction->{taxon_id};
		$params->{scientific_name} = $prediction->{scientific_name};
		$params->{scientific_name} .= " $1"
		    if defined $params->{my_label} && $params->{my_label} =~ /^\s*(\S.*?)\s*$/;

		printf STDERR "[%s] organism predicted at %s level: %s (taxon %s) " .
		    "from reference %s [%s], ANI %.2f%%, AF_ref %.1f%%%s\n",
		    $predictor->status, $prediction->{tier},
		    $params->{scientific_name}, $params->{taxonomy_id},
		    $prediction->{genome_id}, $prediction->{reference_name},
		    $prediction->{ani}, $prediction->{af_ref},
		    ($prediction->{tier} eq 'genus' ? "; $prediction->{tier_reason}" : "");
	    }
	    elsif (!defined $lineage_consistent)
	    {
		$action = 'lineage_unresolved_kept_supplied';
		print STDERR "[" . $predictor->status . "] could not resolve lineages to " .
		    "compare supplied taxonomy_id $params->{taxonomy_id} against " .
		    "predicted $prediction->{taxon_id}; proceeding with the supplied " .
		    "value\n";
	    }
	    elsif ($lineage_consistent)
	    {
		#
		# Case (a). A genus-level prediction is consistent with any supplied
		# species inside that genus.
		#
		$action = 'confirmed';
		printf STDERR "[%s] supplied taxonomy_id %s is consistent with the " .
		    "predicted organism %s (%s level, ANI %.2f%%, exact id match: %s)\n",
		    $predictor->status, $params->{taxonomy_id},
		    $prediction->{scientific_name}, $prediction->{tier}, $prediction->{ani},
		    ($exact_match ? "yes" : "no");
	    }
	    else
	    {
		#
		# Case (b).
		#
		$action = "mismatch_$on_mismatch";
		my($use_id, $use_name) = eval {
		    $predictor->handle_mismatch(
			supplied_taxon_id => $params->{taxonomy_id},
			supplied_name     => $supplied_name,
			prediction        => $prediction,
			on_mismatch       => $on_mismatch);
		};
		if ($@)
		{
		    $fatal = $@;
		}
		else
		{
		    #
		    # Taxon and name move together, or not at all. The submitter's
		    # label (the part of the name after the organism) is kept.
		    #
		    $taxon_overridden = 1 if $use_id ne $params->{taxonomy_id};
		    $params->{taxonomy_id}     = $use_id;
		    $params->{scientific_name} = $predictor->replace_organism_name(
			$params->{scientific_name}, $supplied_name, $use_name) if defined $use_name;
		}
	    }
	}
    }

    #
    # One structured line per job, for the mismatch report. Greppable on
    # [ORGANISM_PREDICTION] and parseable as JSON.
    #
    my $minhash = eval {
	$predictor->summarize_minhash($api, $params->{similar_genome_finder}, $supplied_id, $prediction);
    };
    print STDERR "[" . $predictor->status . "] could not summarize similar_genome_finder: $@" if $@;

    my $organism_report = $predictor->build_report(
	mode               => $prediction_mode,
	minhash            => $minhash,
	supplied_taxon     => $supplied_raw,
	supplied_name      => $supplied_name,
	supplied_reason    => $supplied_reason,
	prediction         => $prediction,
	exact_match        => $exact_match,
	lineage_consistent => $lineage_consistent,
	mismatch_policy    => $mismatch_policy,
	on_mismatch        => $on_mismatch,
	final_taxon       => ($fatal ? undef : $params->{taxonomy_id}),
	final_name         => ($fatal ? undef : $params->{scientific_name}),
	action             => $action);
    $predictor->emit_report($organism_report);
    $params->{_organism_prediction} = $organism_report;

    die $fatal if $fatal;

    #
    # If we are missing domain and/or genetic code, look up from taxonomy.
    # When the prediction replaced the supplied taxon, any submitted domain
    # and code belonged to the replaced organism, so look them up regardless.
    #
    my $def = $api->get_taxon_metadata($params->{taxonomy_id});

    print Dumper($params);
    if ($taxon_overridden ||
	!$params->{domain} ||
	$params->{domain} eq 'auto' ||
	!$params->{code})
    {
	print Dumper($def);
	$params->{domain} = $def->{domain} if $def->{domain};
	$params->{code} = $def->{genetic_code} if $def->{genetic_code};
    }
	
    print Dumper($params);
    #
    # Construct genome object metadata and create a new genome object.
    #

    my $meta = {
	scientific_name => $params->{scientific_name},
	genetic_code => $params->{code},
	domain => $params->{domain},
	($def->{taxon_lineage} ? (ncbi_lineage => $def->{taxon_lineage}) : ()),
	($params->{taxonomy_id} ? (ncbi_taxonomy_id => $params->{taxonomy_id}) : ()),
	($user_id ? (owner => $user_id) : ()),

    };
    my $genome = $core->impl->create_genome($meta);

    #
    # Determine workspace paths for our input and output.
    # $ws and $input_path were established above, with the contig download.
    #

    my $output_folder = $app->result_folder();

    my $output_base = $params->{output_file};

    if (!$output_base)
    {
	$output_base = basename($input_path);
    }

    #
    # Read contig data. The file was downloaded once, above, so that organism
    # prediction could run before the genome object was created.
    #

    my $contig_data_fh = open_contigs("$temp");


    #
    # Use the state-machine based parser ported from RAST.
    #
    # We do an initial read of the first 64 ids of the file. If any of these are too long,
    # we compute the longest common substring to use to replace clip out and replace with
    # "contigs_" to shorten the ids. We keep a mapping from new name to original name.
    #
    # This should hit the common cases where we have IDs like
    # SA-B-8-4-Neg_un-mapped_reads_[SA-B-8-4-Neg_S4_L001_R1_001]_(paired)_contig_22
    # or
    # Ahmed6_GGACTCC_L005_R1_001_(paired)_trimmed_(paired)_merged_contig_1
    #

    my $n = 0;
    my @ids;
    my $too_long = 0;
    my $max_id_len = 60;
    parse_fasta($contig_data_fh, undef, sub {
	my($id, $seq) = @_;
	push(@ids, $id);
	$too_long++ if length($id) > $max_id_len;
	$n++;
	return ($n < 64);
    });
    close($contig_data_fh);

    if ($n == 0)
    {
	die "No contigs loaded from $temp $input_path\n";
    }

    my %name_map;
    my %name_rev_map;
    my $lcs;
    my $qlcs;
    if ($too_long)
    {
	if ($n == 1)
	{
	    $name_map{'contig'} = $ids[0];
	    $name_rev_map{$ids[0]} = 'contig';
	}
	else
	{
	    BuildString(@ids);
	    my $tree = BuildTree();
	    $lcs = LongestCommonSubstring($tree);
	    $qlcs = quotemeta($lcs);
	    print STDERR "Shortening contig names using longest substring '$lcs'\n";
	}
    }

    #
    # Reopen the file and load data, remapping names if necessary.
    #

    $contig_data_fh = open_contigs("$temp");

    my $n = 0;
    my @contigs;
    parse_fasta($contig_data_fh, undef, sub {
	my($id, $seq) = @_;

	my $orig_id = $id;
	my @orig;
	if ($name_rev_map{$id})
	{
	    $id = $name_rev_map{$id};
	    @orig = (original_id => $orig_id);
	}
	elsif ($lcs)
	{
	    $id =~ s/$qlcs/contig_/;
	    @orig = (original_id => $orig_id);
	}
	if (length($id) > $max_id_len)
	{
	    die "Contig id $orig_id too long even after shortening to $id via longest substring $lcs\n";
	}

	push(@contigs, { id => $id, dna => $seq, @orig });
	$n++;
	return 1;
    });
    close($contig_data_fh);
    $core->impl->add_contigs($genome, \@contigs);

    local $Bio::KBase::GenomeAnnotation::Service::CallContext = $core->ctx;

    my $result;

    #
    # Set overrides. This was the method used before we stashed
    # params in the context; unclear if this is better or worse from
    # an engineering perspective.
    #
    my $override;
    if (my $ref = $params->{reference_genome_id})
    {
	$override->{evaluate_genome} =  {
	    evaluate_genome_parameters => { reference_genome_id => $ref },
	};
    }
    if (my $ref = $params->{reference_virus_name})
    {
	$override->{call_features_vigor4} =  {
	    vigor4_parameters => { reference_name => $ref },
	};
    }
    
    $result = $core->run_pipeline($genome, $params->{workflow}, $params->{recipe}, $override);

    my($gto_path, $index_queue_id) = $core->write_output($genome, $result, {}, undef,
							 $params->{public} ? 1 : 0,
							 $params->{queue_nowait} ? 1 : 0,
							 $params->{skip_indexing} ? 1 : 0);

    #
    # Write the genome quality data as a standalone JSON file to aid in downstream
    # quality summarization.
    #

    if (ref($result->{quality}) && !$app->donot_create_result_folder())
    {
	$ws->save_data_to_file($json->encode($result->{quality}),
		           {}, "$output_folder/quality.json", "json", 1, 1, $core->token);
    }

    #
    # Determine if we are one of a peer group of jobs that was started
    # on behalf of a parent job. If we are, and if we are the last job running,
    # invoke post-parent-job processing.
    #

    my $parent_output_folder;
    my $run_last;
    
    if (my $parent = $params->{_parent_job})
    {
	my $dsn = "DBI:mysql:database=" . db_name . ";host=" . db_host;
	my $dbh = DBI->connect($dsn, db_user, db_pass, { RaiseError => 1, AutoCommit => 0 });
	
	#
	# Save information about this genome.
	#
	$dbh->do(qq(INSERT INTO GenomeAnnotation_JobDetails (job_id, parent_job, genome_id, genome_name, gto_path)
		    VALUES (?, ?, ?, ?, ?)), undef,
		 $app->task_id, $parent, $genome->{id}, $genome->{scientific_name}, $gto_path);
	$dbh->commit();

	my $sth = $dbh->prepare(qq(SELECT children_created, children_completed, parent_app, app_spec, app_params
				   FROM JobGroup
				   WHERE parent_job = ?
				   FOR UPDATE));
	my $res = $sth->execute($parent);
	my $last_job = 0;
	my($created, $completed, $app, $spec, $params);
	if ($res != 1)
	{
	    warn "Missing parent job $parent in database\n";
	}
	else
	{
	    ($created, $completed, $app, $spec, $params) = $sth->fetchrow_array();
	    print "Created=$created completed=$completed\n";

	    eval {
		my $params_dat = decode_json($params);
		$parent_output_folder = $params_dat->{output_path} . "/." . $params_dat->{output_file};
		print STDERR "Found parent output folder $parent_output_folder\n";
	    };
	    if ($@)
	    {
		warn "Error parsing parent job params data : $@\n$params\n";
	    }

	    if ($completed == $created - 1)
	    {
		print "We are the last one out!\n";
		$last_job = 1;
	    }
	    elsif ($completed < $created - 1)
	    {
		print "Not so many gone\n";
	    }
	    else
	    {
		warn "completed=$completed created=$created - should not happen here\n";
	    }
	    my $n = $dbh->do(qq(UPDATE JobGroup
				SET children_completed = children_completed + 1
				WHERE parent_job = ? AND children_completed = ?), undef,
			     $parent, $completed);
	    if ($n == 0)
	    {
		print "Failed on update - not us!\n";
	    }
	    else
	    {
		print "Completed with n=$n\n";
	    }
	}

	$dbh->commit();
	if ($last_job)
	{
	    #
	    # We defer this until we are outside this block so we can ensure
	    # the genome report is written first.
	    #
	    $run_last = sub { $core->run_last_job_processing($parent, $app, $spec, $params); };
	}
    }

    
    #
    # Do last-job processing if needed.
    #
    &$run_last if $run_last;

    $core->ctx->stderr(undef);

    return {
	gto_path => $gto_path,
	index_queue_id => $index_queue_id,
	genome_id => $genome->{id},
    };
}

sub open_contigs
{
    my($file) = @_;
    my $contig_data_fh;
    open($contig_data_fh, "<", $file) or die "Cannot open contig file $file: $!";

    #
    # Read first block to see if this is a gzipped file.
    #
    my $block;
    $contig_data_fh->read($block, 256);
    if ($block =~ /^\037\213/)
    {
	#
	# Gzipped. Close and reopen from gunzip.
	#
	
	close($contig_data_fh);
	undef $contig_data_fh;
	open($contig_data_fh, "-|", "gzip", "-d", "-c", $file) or die "Cannot open gzip from $file: $!";
    }
    else
    {
	$contig_data_fh->seek(0, 0);
    }

    return $contig_data_fh;
}

sub run_seedtk_cmd
{
    my(@cmd) = @_;
    local $ENV{PATH} = seedtk . "/bin:$ENV{PATH}";
    my $ok = IPC::Run::run(@cmd);
    $ok or die "Failure $? running seedtk cmd: " . Dumper(\@cmd);
}
