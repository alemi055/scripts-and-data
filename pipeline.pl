#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Config::General;
use Cwd;

#
# NGS PIPELINE IN PERL
#
# credit: Ilya Plyusnin, University of Helsinki, Ilja.Pljusnin@helsinki.fi
#

my $install_dir 	= "/home/alemi055/lazypipe_ete2021";
my $perl_scripts 	= "/home/alemi055/lazypipe_ete2021/lazypipe/perl";
my $cpp_tools 		= "/home/alemi055/lazypipe_ete2021/lazypipe/cpp";
my $R_scripts 		= "/home/alemi055/lazypipe_ete2021/lazypipe/R";
my $config_default_file	= "/home/alemi055/lazypipe_ete2021/lazypipe/pipeline.default.config";
my $config_my_file	= "/home/alemi055/lazypipe_ete2021/lazypipe/pipeline.my.config";
my $ktImportText	= "~/lazypipe_ete2021/Krona/KronaTools/scripts/ImportText.pl";	

my $usage= 	"USAGE: $0 -1 file [-2 file --hostgen file --hgtaxid int --inlen int] --t|numth num -p|pipe str --pre trimm|fastp --ass velvet|megahit --gen mga|mgm --ann blast|sans|centrifuge -w|weights str [-r|res dir -l|label str]\n".
		"\n".
		"-1|read1 file    : paired reads, forward reads\n".
		"-2|read2 file    : paired reads, backward reads. If ommited tries to guess file name from --read1\n".
		"-r|res   dir     : results will be printed to res/\$opt{'label'}. If ommited print to RES/\$opt{'label'}\n".
		"-l|label str     : name result files using this label. If ommited use prefix of the read1 filename\n".
		"-t|numth num     : number of threads\n".
		"--pre str        : use fastp|trimmomatic to preprocess reads. Default Trimmomatic\n".
		"--inlen int      : insert length. Default 300\n".
		"--ass str        : assembler: velvet|megahit|spades. Default megahit\n".
		"--gen str        : an initio gene identification: mga|mgm (MetaGeneAnnotator|MetaGeneMark). Default mga\n".
		"--ann str        : gene annotation: blast|sans|centrifuge. Default sans\n".
		"--th_bits int    : threshold for SANSparallel bitscore. Default 120\n".
		"--hostgen file   : *.fna file containing host genome. Filtering is turned on by --hostgen file -p 2\n".
		"--hgtaxid int    : NCBI taxon id for the host genome taxon. When given, hostgen filtered reads will be assigned to this taxon\n".
		"-w|weights str   : model for abundance estimation: taxacount|bitscore|bitscore2. Default: bitscore2\n".
		"-p|pipe str      : pipeline steps to process. These can be either in int:int or int,int[,int] format\n".
		"        1     : filter reads with Trimmomatic|fastp\n".
		"        2     : filter host genome using --hostgen file\n".
		"        3     : assemble\n".
		"        4     : run MetaGeneAnnotator|MetaGeneMark\n".
		"        5     : run BLASTP/SANSparallel/centrifuge\n".
		"        6     : realign reads to contigs\n".
		"        7     : generate abundance + seq annotation tables + Krona graph + sort contigs by taxa\n".
		"        8     : IGV report aligning contigs to refgen (requires local database of viral refgenomes)\n".
		"        9     : add blastn annotations for viral contigs\n".
		"        10    : add blastn annotations for contigs with no db hits in step 5\n".
		"        11    : assembly stats + QC plots\n".		
		"        12    : collect files for sharing\n".		
		"        13    : cleanup (delete intermediate file)\n\n".
		"NOTE: default options will be read from $config_my_file, or, if this is not available, from $config_default_file\n\n";

# Default options defined in configuration file
my $config_file	= ( -e($config_my_file)) ? $config_my_file: $config_default_file;
my $conf = Config::General->new(
   -ConfigFile      => "$config_file",
   -InterPolateVars => 1,
   -LowerCaseNames  => 1,
   -AutoTrue        => 1);
my %opt = $conf->getall();

# Command line options overwrite config file options
GetOptions(\%opt,
        'read1|1=s','read2|2=s','hostgen=s','hgtaxid=i',
	'pipe|p=s','res|r=s','label|l=s','pre=s','ass=s','th_bits=i',
	'centrifuge_th=i','gen=s','ann=s','inlen=i','numth|t=i', 'weights|w=s') or die $usage;
$opt{'ann'} = lc($opt{'ann'});
$opt{'ass'} = lc($opt{'ass'});
$opt{'gen'} = lc($opt{'gen'});
$opt{'weights'} = lc($opt{'weights'});
if(!$opt{'read1'} || !$opt{'pipe'}){
	print STDERR "\nMissing arguments\n\n";
        die $usage;
}
if( $opt{'pre'} =~ /trimm/gi ){
	$opt{'pre'} = 'trimm';
}
elsif( $opt{'pre'} =~ /fastp/gi ){
	$opt{'pre'} = 'fastp';
}
else{
	print STDERR "\ninvalid option --pre $opt{'pre'}. Running with --pre trimmomatic\n";
	$opt{'pre'} = "trimm";
}
	
if( !(($opt{'ass'} eq 'velvet') || ($opt{'ass'} eq 'megahit') || ($opt{'ass'} eq 'spades')) ){
	print STDERR "\ninvalid option --ass $opt{'ass'}. Running with --ass megahit\n";
	$opt{'ass'} = 'megahit';
}
if( !(($opt{'gen'} eq 'mga') || ($opt{'gen'} eq 'mgm') )){
	print STDERR "\ninvalid option --gen $opt{'gen'}. Running with --gen mga\n\n";
	$opt{'gen'} = 'mga';
}
if( !(($opt{'ann'} eq 'blastp') || ($opt{'ann'} eq 'sans') || ($opt{'ann'} eq 'centrifuge'))){
	print STDERR "\ninvalid option --ann $opt{'ann'}. Running with --ann sans\n\n";
	$opt{'ann'} = 'sans';
}

if( !(($opt{'weights'} eq 'taxacount') || ($opt{'weights'} eq 'bitscore') || ($opt{'weights'} eq 'bitscore2')) ){
	print STDERR "\ninvalid option --weights $opt{'weights'}. Running with --weights bitscore2\n\n";
	$opt{'weights'} = 'bitscore2';
}



# GUESSING MISSING PARAMETERS
if( !defined($opt{'read2'}) ){
	$opt{'read2'}= $opt{'read1'};
	if( !($opt{'read2'} =~ s/_R1/_R2/i) ){
		die "ERROR: failed to guess read2 file\n";
	}
}

# CHECKING INPUT FILES
if( !(-e $opt{'read1'}) ){
	die "ERROR: check read1 file: $opt{'read1'} does not exist\n";}
if( !(-e $opt{'read2'}) ){
	die "ERROR: check read2 file: $opt{'read2'} does not exist\n";}
if( defined($opt{'hostgen'}) && !(-e $opt{'hostgen'})){
	die "ERROR: check hostgen file: $opt{'hostgen'} does not exist\n";}

if(!$opt{'label'}){
	$opt{'label'}= $opt{'read1'};
	my @tmp= split(/\//,$opt{'label'});	# remove /*/* prefix
	$opt{'label'} = $tmp[$#tmp];
	@tmp= split(/\./,$opt{'label'});	# remove .*[.*] suffix
	$opt{'label'}= $tmp[0];
	@tmp= split(/_/,$opt{'label'},-1);	# remove _.* suffix
	$opt{'label'}= $tmp[0];	
}	
if($opt{'res'}){
	$opt{'res'}= "$opt{'res'}/$opt{'label'}";
}
else{
	$opt{'res'}= "RES/$opt{'label'}";
}

# DEBUG: PRINT OPT
#foreach my $k(sort keys %opt){ print "$k\t:$opt{$k}\n";} exit(1);

# Setting TMPDIR and R environment
$ENV{TMPDIR} = $opt{'wrkdir'};
if( !(-e "./.Renviron") ){
	system_call("echo \"TMPDIR=$opt{'wrkdir'}\" > .Renviron");
}


# # UPDATING TAXONOMY FILES: this will also load taxonomy on the very first usage
# if($opt{'update'}){
	# system_call("perl $perl_scripts/update_database.pl");
# }



# PIPELINE STEPS
my %pipeh;
for(my $i=1; $i<=12; $i++){  $pipeh{$i}= !1; }
my @tmp= split(/,/,$opt{'pipe'});
foreach my $t(@tmp){
	my @pair= split(/:/,$t);
	if(scalar @pair > 1){
		if($pair[0] < 1  ||  $pair[1] > 13){
			die "invalid pipeline steps: $t: we only understand steps from 1 to 13\n";
		}
		for(my $i=$pair[0]; $i<=$pair[1]; $i++){
			$pipeh{$i} = 1;
		}
	}
	else{
		if($t < 1 || $t > 13){
			die "invalid pipeline steps: $t: we only understand steps from 1 to 13\n";
		}
		$pipeh{$t} = 1;
	}
}


if( !(-e $opt{'res'})){		system_call("mkdir -p $opt{'res'}"); }
if( !(-e $opt{'wrkdir'})){	system_call("mkdir -p $opt{'wrkdir'}"); }
my %timer_hash;
if( (-e "$opt{'res'}/timer.txt") ){
	read_timer("$opt{'res'}/timer.txt",\%timer_hash);
}


if($pipeh{1}){
	my @time1= times;

	print STDERR "\n# PREPROCESS READS\n\n";
	
	if( $opt{'pre'} eq 'trimm' ){
		system_call("java -jar $opt{'trimm_jar'} PE -threads $opt{'numth'} $opt{'read1'} $opt{'read2'} -baseout $opt{'res'}/trimmed $opt{'trimm_par'} &> $opt{'res'}/trimmomatic.log");
	}
	elsif( $opt{'pre'} eq 'fastp' ){
		system_call("fastp --thread $opt{'numth'} -i $opt{'read1'} -I $opt{'read2'} -o $opt{'res'}/trimmed_1P -O $opt{'res'}/trimmed_2P --unpaired1 $opt{'res'}/trimmed_1U --unpaired2 $opt{'res'}/trimmed_2U $opt{'fastp_par'} 2> $opt{'res'}/fastp.log");
	}

	my @time2= times; $timer_hash{1}= timediff(\@time1,\@time2);
}

if($pipeh{2}){
	my @time1= times;
	
	print STDERR "\n# FILTER HOST GENOME\n\n";
	if(!defined($opt{'hostgen'})  ||  !$opt{'hostgen'} ){
		print STDERR "No hostgen specified: no filtering\n";
	}
	else{
		if( !((-e "$opt{'hostgen'}.amb") && (-e "$opt{'hostgen'}.sa"))){
		  system_call("bwa index -a bwtsw $opt{'hostgen'}");
		}
		system_call("bwa mem -t $opt{'numth'} -T $opt{'bwa_t1'} $opt{'hostgen'} $opt{'res'}/trimmed_1P $opt{'res'}/trimmed_2P 1> $opt{'res'}/rg.sam 2> $opt{'res'}/rg.log");
		system_call("samtools view -@ $opt{'numth'} -F4 -Shb $opt{'res'}/rg.sam > $opt{'res'}/rg.bam");
		system_call("samtools sort -@ $opt{'numth'} -n $opt{'res'}/rg.bam  -T $opt{'wrkdir'} -o $opt{'res'}/rg.sorted.bam");
		system_call("samtools view -F4 $opt{'res'}/rg.sorted.bam | perl $perl_scripts/filter_tophits.pl > $opt{'res'}/rg.top.T$opt{'bwa_t1'}.sam");
		system_call("cat $opt{'res'}/rg.top.T$opt{'bwa_t1'}.sam | cut -f1 | sort -T $opt{'wrkdir'} | uniq > $opt{'res'}/rg.reads");
		
		#system_call("perl $perl_scripts/filter_reads.pl $opt{'res'}/trimmed_1P $opt{'res'}/trimmed_2P $opt{'res'}/rg.reads rgflt");
		system_call("$cpp_tools/fqfilt $opt{'res'}/trimmed_1P $opt{'res'}/trimmed_2P $opt{'res'}/rg.reads rgflt");
	}
	
	my @time2= times; $timer_hash{2}= timediff(\@time1,\@time2);	
}

if($pipeh{3}){
	my @time1= times;
	
	print STDERR "\n# ASSEMBLE\n\n";	
	
	my $rgflt = (-e "$opt{'res'}/trimmed_1P_rgflt") && (-e "$opt{'res'}/trimmed_2P_rgflt");
	my $r1 = ($rgflt) ? "$opt{'res'}/trimmed_1P_rgflt" : "$opt{'res'}/trimmed_1P";
	my $r2 = ($rgflt) ? "$opt{'res'}/trimmed_2P_rgflt" : "$opt{'res'}/trimmed_2P";
	
	if($opt{'ass'} eq 'velvet'){
		system_call("velveth $opt{'res'} $opt{'kmer'} -fastq -shortPaired $r1 $r2 &> $opt{'res'}/velveth.log");
		system_call("velvetg $opt{'res'} -exp_cov auto -cov_cutoff $opt{'cov_cutoff'} -ins_length $opt{'inlen'} &> $opt{'res'}/velvetg.log");
		system_call("mv $opt{'res'}/contigs.fa $opt{'res'}/contigs.velvet.fa");
		system_call("cat $opt{'res'}/contigs.velvet.fa | sed 's/>NODE_\\([0-9]\\+\\)_length_\\([0-9]\\+\\)_cov_\\([\\.0-9]\\+\\)/>contig=\\1_length=\\2_coverage=\\3/' 1> $opt{'res'}/contigs.fa");
	}
	elsif($opt{'ass'} eq 'megahit'){
		system_call("rm -fR $opt{'res'}/assembler_out"); # Megahit will complain if that dir exists
		system_call("$opt{'megahit'} -t $opt{'numth'} -1 $r1 -2 $r2 --out-dir $opt{'res'}/assembler_out &> $opt{'res'}/assembler.log");
		system_call("mv $opt{'res'}/assembler_out/final.contigs.fa $opt{'res'}/contigs.megahit.fa");
		system_call("cat $opt{'res'}/contigs.megahit.fa | sed 's/>\\([[:alnum:]]\\+\\)_\\([0-9]\\+\\)/>contig=\\1.\\2/' | sed 's/multi/coverage/' | sed 's/len/length/' | sed 's/\\s\\+/_/g' 1> $opt{'res'}/contigs.fa");
		#system_call("rm -fR $opt{'res'}/assembler_out");
	}
	elsif($opt{'ass'} eq 'spades'){
		system_call("rm -fR $opt{'res'}/assembler_out");
		system_call("ln -f $r1 $r1.fq");
		system_call("ln -f $r2 $r2.fq");
		system_call("$opt{'spades'} -t $opt{'numth'} -1 $r1.fq -2 $r2.fq -o $opt{'res'}/assembler_out &> $opt{'res'}/assembler.log");
		system_call("cat $opt{'res'}/assembler_out/contigs.fasta | sed 's/>NODE_/>contig=/' | sed 's/length_/length=/' | sed 's/cov_/coverage=/' 1> $opt{'res'}/contigs.fa");
		system_call("mv $opt{'res'}/assembler_out/scaffolds.fasta $opt{'res'}/scaffolds.fa");
		#system_call("rm -fR $opt{'res'}/assembler_out");
	}
	else{
		die "ERROR: invalid --ass $opt{'ass'}";
	}
	
	my @time2= times; $timer_hash{3}= timediff(\@time1,\@time2);
}

if($pipeh{4}){
	my @time1= times;
	
	if($opt{'gen'} eq 'mga'){
	print STDERR "\n# RUN METAGENEANNOTATOR\n\n";
	system_call("$opt{'mga'} $opt{'res'}/contigs.fa -m 1> $opt{'res'}/genes.mga");
	system_call("perl $perl_scripts/metagene_parser.pl $opt{'res'}/contigs.fa $opt{'res'}/genes.mga --format dna --minlen $opt{'th_glen'} > $opt{'res'}/genes.dna.fa");
	system_call("perl $perl_scripts/metagene_parser.pl $opt{'res'}/contigs.fa $opt{'res'}/genes.mga --format aa --minlen $opt{'th_glen'}  > $opt{'res'}/genes.aa.fa");
	}
	elsif($opt{'gen'} eq 'mgm'){
	print STDERR "\n# RUN METAGENEMARK\n\n";
	system_call("$opt{'mgm_root'}/gmhmmp -m $opt{'mgm_root'}/MetaGeneMark_v1.mod -f 3 -a -A $opt{'res'}/genes.aa.fa.tmp -d -D $opt{'res'}/genes.dna.fa.tmp -o $opt{'res'}/genes.gff3 $opt{'res'}/contigs.fa");
	# converting seqid lines to expected format
	system_call("cat $opt{'res'}/genes.aa.fa.tmp | sed -r 's/^>(\\w+)\\|([A-Za-z0-9.]+)\\|(\\w+)\\|([+-])\\|([0-9]+)\\|([0-9]+)\\t(\\S+)\$/\\7_start=\\5_end=\\6_strand=\\4_model=\\2/' > $opt{'res'}/genes.aa.fa");
	system_call("cat $opt{'res'}/genes.dna.fa.tmp | sed -r 's/^>(\\w+)\\|([A-Za-z0-9.]+)\\|(\\w+)\\|([+-])\\|([0-9]+)\\|([0-9]+)\\t(\\S+)\$/\\7_start=\\5_end=\\6_strand=\\4_model=\\2/' > $opt{'res'}/genes.dna.fa");
	}
	else{
		die "ERROR: invalid --gen $opt{'gen'}";
	}	
	
	my @time2= times; $timer_hash{4}= timediff(\@time1,\@time2);
}

if($pipeh{5}){
	my @time1= times;

	if(lc($opt{'ann'}) eq 'blastp'){
	print STDERR "\n# RUN BLASTP\n\n";
	system_call("blastp -num_threads $opt{'numth'} -evalue 10 -use_sw_tback -max_target_seqs 1 -outfmt '6 qseqid sacc bitscore pident length stitle sscinames staxids ' -query $opt{'res'}/genes.aa.fa -db $opt{'data'}/blastdb/nr 1> $opt{'res'}/db_hits1.txt 2> $opt{'res'}/blast.log");
	system_call("cat $opt{'res'}/db_hits1.txt | grep -v \"^#\" | sort -T $opt{'wrkdir'} -t\$'\\t' -k1,1 -uV > $opt{'res'}/db_hits2.txt");
	system_call("perl $perl_scripts/get_seq_ids.pl $opt{'res'}/genes.aa.fa 1> $opt{'res'}/genes.ids");
	system_call("join -t \$'\\t' -j1 -a1 <(sort -T $opt{'wrkdir'} $opt{'res'}/genes.ids) <(sort -T $opt{'wrkdir'} -t\$'\\t' -k1,1 $opt{'res'}/db_hits2.txt) | perl $perl_scripts/tab.pl > $opt{'res'}/db_hits3.txt");
	
	system_call("cat $opt{'res'}/db_hits2.txt | awk -F'\\t' '{print \$1\"\\t\"\$7\"\\t\"\$3}' | sed 's/^contig=\\([A-Za-z0-9\\.]\\+\\)_[^\\t]*/\\1/' | sort -T $opt{'wrkdir'} | uniq > $opt{'res'}/contid_taxid_score");
	}
	
	elsif( lc($opt{'ann'}) eq 'centrifuge'){
	print STDERR "\n# RUN CENTRIFUGE\n\n";
	system_call("centrifuge --threads $opt{'numth'} -f -x $opt{'centrifuge_db'} -U $opt{'res'}/contigs.fa --report-file $opt{'res'}/centrifuge.report -S $opt{'res'}/db_hits1.txt &> $opt{'res'}/centrifuge.log");
	system_call("cat $opt{'res'}/db_hits1.txt | awk -F'\\t' '\$6>=$opt{'centrifuge_th'}' 1> $opt{'res'}/tmp");	# First filter hits by score!
	system_call("cat $opt{'res'}/tmp | awk -F'\\t' '{print \$1\"\\t\"\$2\"\\t\"\$4\"\\t\"\$6\"\\t\"\$7\"\\t\"\$3}' > $opt{'res'}/db_hits2.txt");
	system_call("perl $perl_scripts/get_seq_ids.pl $opt{'res'}/contigs.fa 1> $opt{'res'}/contig.ids");
	system_call("echo -e \"qid\\tsid\\tscore\\thitLength\\tqueryLength\\ttaxid\" > $opt{'res'}/db_hits3.txt");
	system_call("join -t \$'\\t' -j1 -a1 <(sort -T $opt{'wrkdir'} $opt{'res'}/contig.ids) <(sort -T $opt{'wrkdir'} -t\$'\\t' -k1,1 $opt{'res'}/db_hits2.txt) | perl $perl_scripts/tab.pl >> $opt{'res'}/db_hits3.txt");
	
	system_call("cat $opt{'res'}/db_hits2.txt | awk -F'\\t' '{print \$1\"\\t\"\$6\"\\t\"\$3}'  | sed 's/^contig=\\([A-Za-z0-9\\.]\\+\\)_[^\\t]*/\\1/' | sort -T $opt{'wrkdir'} | uniq > $opt{'res'}/contid_taxid_score");
	}
		
	elsif( lc($opt{'ann'}) eq 'sans'){
	print STDERR "\n# RUN SANS\n\n";
	
#	system_call("python $opt{'sanspanz'} -i $opt{'res'}/genes.aa.fa -o $opt{'res'}/db_hits1.txt &> $opt{'res'}/sans.log"); # newer sans returns hits sorted by descending bits
	system_call("awk -F'\\t' '\$1>$opt{'th_bits'}' $opt{'res'}/db_hits1.txt > $opt{'res'}/tmp0");	# filter hits by sans bit-score
	system_call("awk -F'\\t' '\$3==0'  $opt{'res'}/tmp0 > $opt{'res'}/tmp1");			# filter self hits marked as "1" in col3
	#system_call("cat $opt{'res'}/tmp1 | sort -T $opt{'wrkdir'} -t\$'\\t' -suV -k4,4 > $opt{'res'}/tmp2");	# sort col4:qpid , report only the top hit
	system_call("cat $opt{'res'}/tmp1 | sort -T $opt{'wrkdir'} -t\$'\\t' -sV -k4,4 | perl $perl_scripts/filter_tophits2.pl --qcol 4 --bitscol 1 --ties 1> $opt{'res'}/tmp2"); # select tophits including ties

	
	system_call("cat $opt{'res'}/tmp2 | awk -F'\\t' '{print \$4\"\\t\"\$5\"\\t\"\$1\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10}' > $opt{'res'}/tmp3");	# rearrange cols
	system_call("$cpp_tools/linkids2 -v -d $opt{'res'}/tmp3 -c 2 -m $opt{'taxonomy_idmap'} 1> $opt{'res'}/db_hits2.txt");
	
	system_call("perl $perl_scripts/get_seq_ids.pl $opt{'res'}/genes.aa.fa 1> $opt{'res'}/genes.ids");
	system_call("echo -e \"qpid\\tspid\\tbits\\tpide\\tlali\\tdesc\\ttaxid\" > $opt{'res'}/db_hits3.txt");
	system_call("join -t \$'\\t' -j1 -a1 <(sort -T $opt{'wrkdir'} $opt{'res'}/genes.ids) <(sort -T $opt{'wrkdir'} -t\$'\\t' -k1,1 $opt{'res'}/db_hits2.txt) | perl $perl_scripts/tab.pl >> $opt{'res'}/db_hits3.txt");
	
	system_call("cat $opt{'res'}/db_hits2.txt | awk -F'\\t' '{print \$1\"\\t\"\$7\"\\t\"\$3}' | sed 's/^contig=\\([A-Za-z0-9\\.]\\+\\)_[^\\t]*/\\1/' | sort -T $opt{'wrkdir'} | uniq > $opt{'res'}/contid_taxid_score"); # fields: contid + taxid + bitscore
	#system_call("rm -f $opt{'res'}/tmp*");
	}
	else{
		die "ERROR: invalid --ann $opt{'ann'}";
	}
	
	my @time2= times; $timer_hash{5}= timediff(\@time1,\@time2);
}

if($pipeh{6}){
	my @time1= times;
	
	print STDERR "\n# REALIGN READS TO CONTIGS\n\n";
	system_call("bwa index $opt{'res'}/contigs.fa &> $opt{'res'}/bwa.log");
	if( -e "$opt{'res'}/trimmed_1P_rgflt"){
	   system_call("bwa mem -t $opt{'numth'} -T $opt{'bwa_t2'} $opt{'res'}/contigs.fa $opt{'res'}/trimmed_1P_rgflt $opt{'res'}/trimmed_2P_rgflt 1> $opt{'res'}/contigs.sam 2> $opt{'res'}/bwa.log");
	}
	else{
	   system_call("bwa mem -t $opt{'numth'} -T $opt{'bwa_t2'} $opt{'res'}/contigs.fa $opt{'res'}/trimmed_1P $opt{'res'}/trimmed_2P 1> $opt{'res'}/contigs.sam 2> $opt{'res'}/bwa.log");
	}
	system_call("samtools view -@ $opt{'numth'} -F4 -bSh $opt{'res'}/contigs.sam -o $opt{'res'}/contigs.bam");
	system_call("samtools sort -@ $opt{'numth'} -n $opt{'res'}/contigs.bam -T $opt{'wrkdir'} -o $opt{'res'}/contigs.sorted.bam");	# sort by read name for filter_tophits.pl
	system_call("samtools view -F4 -h $opt{'res'}/contigs.sorted.bam | perl $perl_scripts/filter_tophits.pl > $opt{'res'}/contigs.flt.sam");
	system_call("samtools view -@ $opt{'numth'} -bS $opt{'res'}/contigs.flt.sam -o $opt{'res'}/contigs.flt.bam");
	system_call("samtools sort -@ $opt{'numth'} $opt{'res'}/contigs.flt.bam -T $opt{'wrkdir'} -o $opt{'res'}/contigs.flt.sorted.bam");		# sort by chromosome posit
	system_call("samtools index -@ $opt{'numth'} $opt{'res'}/contigs.flt.sorted.bam");
	system_call("samtools idxstats $opt{'res'}/contigs.flt.sorted.bam > $opt{'res'}/contigs.bam.stats");
	
	my @time2= times; $timer_hash{6}= timediff(\@time1,\@time2);
}

if($pipeh{7}){
	my @time1= times;
	
	print STDERR "\n# GENERATE SUMMARY TABLES + TAXONOMIC PROFILE + KRONA GRAPH + SORT CONTIGS BY TAXA\n\n";
	
	# ABUNDANCE TABLES
	my $hostgen_pars1 = "";
	my $hostgen_pars2 = "";
	my $host_taxid = "";
	if($opt{'hostgen'} && $opt{'hgtaxid'}){ 
		my $readn_flt 	= (`cat $opt{'res'}/trimmed_1P | wc -l`)/4;
		my $readn_rgflt	= $readn_flt;
		if(-e "$opt{'res'}/trimmed_1P_rgflt"){
			$readn_rgflt = (`cat $opt{'res'}/trimmed_1P_rgflt | wc -l`)/4;
		}
		my $rgabund 	= $readn_flt - $readn_rgflt;
		$hostgen_pars1	= "--rgabund $rgabund --rgtaxid $opt{'hgtaxid'}";
		$hostgen_pars2	= "--rgtaxid $opt{'hgtaxid'}";
		$host_taxid	= $opt{'hgtaxid'};
	}
	system_call("perl $perl_scripts/get_abund_table.pl $hostgen_pars1 --weights $opt{'weights'} --cont_score_tail 5 $opt{'res'}/contid_taxid_score $opt{'res'}/contigs.bam.stats 1> $opt{'res'}/readn_taxid");
	system_call("perl $perl_scripts/link_taxinfo.pl -h -v -i $opt{'res'}/readn_taxid $opt{'taxonomy_nodes'} $opt{'taxonomy_names'} species,genus,family,order,class,phylum,superkingdom 1> $opt{'res'}/abund_table.txt");
	
	# TAXONOMIC PROFILE
	system_call("perl $perl_scripts/abundtable2taxprofile.pl --sample $opt{'label'} --tail 2 $hostgen_pars2 $opt{'res'}/abund_table.txt 1> $opt{'res'}/taxprof.txt");


	# CREATE KRONA GRAPH	
	system_call("perl $perl_scripts/abundtable2krona.pl --tail 2 $hostgen_pars2 $opt{'res'}/abund_table.txt 1> $opt{'res'}/krona_data.txt");
#	system_call("ktImportText $opt{'res'}/krona_data.txt -o $opt{'res'}/krona_graph.html -u \"http://krona.sourceforge.net\"");
  system_call("$ktImportText $opt{'res'}/krona_data.txt -o $opt{'res'}/krona_graph.html -u \"http://krona.sourceforge.net\"");

	# CONTIG ANNOTATION TABLES
	system_call("perl $perl_scripts/link_taxinfo.pl -h -i $opt{'res'}/db_hits3.txt $opt{'taxonomy_nodes'} $opt{'taxonomy_names'} species,genus,family,order,superkingdom | perl $perl_scripts/split_qid.pl --header > $opt{'res'}/annot_table.txt");

	# EXCEL TABLES: abund and annot
	system_call("Rscript $R_scripts/print_abund_table.R $opt{'res'}/abund_table.txt $opt{'res'}/annot_table.txt $opt{'res'}/summary1-$opt{'label'}.xlsx 2 $host_taxid"); # NOTE: contign field is calc from annot_table.txt
	system_call("Rscript $R_scripts/print_annot_table.R $opt{'res'}/annot_table.txt $opt{'res'}/summary2-$opt{'label'}.xlsx");
	
	# NCA ANALYSIS: PREDICT CELLULAR HOST
	if($opt{'nca_predict_host'} && defined($opt{'nca_db'}) ){
		system_call("perl $perl_scripts/nca_link_hostclass.pl $opt{'res'}/annot_table.txt $opt{'res'}/contigs.fa $opt{'nca_db'} $opt{'res'} 1> $opt{'res'}/annot_table.tmp");
		system_call("mv $opt{'res'}/annot_table.tmp $opt{'res'}/annot_table.txt");
	}	
	
	# SORT CONTIGS TO DIRS USING TAXONOMY CLASSIFICATION
	system_call("perl $perl_scripts/sort_contigs_bytaxa.pl --res $opt{'res'}/contigs --blast $opt{'res'}/annot_table.txt --cont $opt{'res'}/contigs.fa");

	
	my @time2= times; $timer_hash{7}= timediff(\@time1,\@time2);
}


if($pipeh{8}){
	print STDERR "\n# IGV REPORT ALIGNING CONTIGS TO REFGEN\n\n";
	system_call("perl $perl_scripts/refgen_graph.pl --numth $opt{'numth'} --wrk $opt{'wrkdir'} --th_csumq 2  -v $opt{'res'}");	
}

if($pipeh{9}){
	print STDERR "\n# ADD BLASTN ANNOTATION FOR VIRAL CONTIGS\n\n";
	system_call("perl $perl_scripts/get_blastn_vi_ann.pl --numth $opt{'numth'} --wrk $opt{'wrkdir'}  -v $opt{'res'}");
}

if($pipeh{10}){
	print STDERR "\n# ADD BLASTN ANNOTATION FOR UNANNOTATED CONTIGS\n\n";
	system_call("perl $perl_scripts/get_blastn_un_ann.pl --numth $opt{'numth'} --wrk $opt{'wrkdir'}  -v $opt{'res'}");
}

if($pipeh{11}){
	my @time1= times;
	
	print STDERR "\n# ASSEMBLY STATS + QUALITY CONTROL PLOTS\n\n";
	
	system_call("perl $perl_scripts/assembly_stats_B.pl --format col,names --fastq $opt{'read1'} $opt{'res'} 1> $opt{'res'}/assembly.stats.txt");
	system_call("perl $perl_scripts/assembly_stats.pl $opt{'res'}/contigs.fa $opt{'res'}/genes.dna.fa mean,sum,N50,LN500,Lbp500,LN1000,Lbp1000 col,names 1>> $opt{'res'}/assembly.stats.txt");

	system_call("cat $opt{'read1'} | awk '{if(NR%4==2) print length(\$1)}' | sort -T $opt{'wrkdir'} -n > $opt{'res'}/r1.readlen");
	system_call("cat $opt{'read2'} | awk '{if(NR%4==2) print length(\$1)}' | sort -T $opt{'wrkdir'} -n > $opt{'res'}/r2.readlen");
	system_call("cat $opt{'res'}/contigs.fa | awk '{if(NR%2==0) print length(\$1)}' | sort -T $opt{'wrkdir'} -n > $opt{'res'}/contig.len");
	
	system_call("Rscript $R_scripts/qc_plots.R seqlen $opt{'res'}/r1.readlen $opt{'res'}/qc.r1hist.jpeg");
	system_call("Rscript $R_scripts/qc_plots.R seqlen $opt{'res'}/r2.readlen $opt{'res'}/qc.r2hist.jpeg");
	system_call("Rscript $R_scripts/qc_plots.R seqlen $opt{'res'}/contig.len $opt{'res'}/qc.conhist.jpeg");
	system_call("Rscript $R_scripts/qc_plots.R readn $opt{'res'}/assembly.stats.txt $opt{'res'}/qc.rsurv.jpeg");
	
	my @time2= times; $timer_hash{9}= timediff(\@time1,\@time2);
}


if($pipeh{12}){
	my @time1= times;
	
	print STDERR "\n# PACK FILES FOR SHARING\n\n";

	my $wd= getcwd();
	chdir($opt{'res'}) or die "$!";	# chdir needed to keep dir structure made by tar simple
	system_call("rm -fR $opt{'label'}");
	system_call("mkdir -p $opt{'label'}");
	my @files_share= ("abund_table.txt",
			"annot_table.txt",
			"taxprof.txt",
			"contigs",
			"contigs.fa",
#			"assembly.stats.txt",
			"genes.aa.fa",
			"genes.dna.fa",
#			"*.jpeg",
			"*.html",
			"*.xlsx");
	if( -e "blastn_vi_annot_table.txt" ){ 		push(@files_share, "blastn_vi_annot_table.txt"); }
	if( -e "blastn_un_annot_table.txt" ){ 		push(@files_share, "blastn_un_annot_table.txt"); }
	my $files_share_str = join(" ",@files_share);
	system_call("cp -r $files_share_str $opt{'label'}/");
	system_call("tar -czf $opt{'label'}.tar.gz $opt{'label'}");
	chdir($wd) or die "$!";
	my @time2= times; $timer_hash{10}= timediff(\@time1,\@time2);
}


if($pipeh{13}){
	# CLEANUP
	system_call("rm -fR $opt{'res'}/*.bam  $opt{'res'}/*.bai");
	system_call("rm -fR $opt{'res'}/*Graph*");
	system_call("rm -fR $opt{'res'}/Roadmaps");
	system_call("rm -fR $opt{'res'}/*tmp*");
	system_call("rm -fR $opt{'res'}/assembler_out");
}		
	

#system_call("rm -fR $opt{'wrkdir'}/* ");	# interfiers with multiple runs
print_timer("$opt{'res'}/timer.txt",\%timer_hash);


sub system_call{
	my $call= shift;
	print STDERR "$call\n";
	my @args= ("bash","-c",$call);
	system(@args) == 0 or die "";
}

sub read_timer{
	my $file= shift;
	my $timer_hashp= shift;
	open(IN,"<$file") or die;
	while(my $l=<IN>){
		chomp($l);
		my($var,$val) = split(/\t/,$l);
		$timer_hashp->{$var}= $val;
	}close(IN);
}

sub print_timer{
	my $file= shift;
	my %timer_hash= %{shift()};
	open(OUT,">$file") or die;
	foreach my $k(sort {$a <=> $b} keys %timer_hash){
		print OUT "$k\t",$timer_hash{$k},"\n";
	}
	close(OUT);
}
sub timediff{
	my @t1= @{shift() };
	my @t2= @{shift() };
	my $td= ($t2[0]-$t1[0]) + ($t2[1]-$t1[1]) + ($t2[2]-$t1[2]) + ($t2[3]-$t1[3]);
	return $td;
}

