#!/usr/bin/env nextflow
/*
========================================================================================
                         PhilPalmer/exomedepth
========================================================================================
 PhilPalmer/exomedepth Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/PhilPalmer/exomedepth
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     PhilPalmer/exomedepth //v
    =========================================
    Description:
    Exome Depth Per-Chromosome Analysis

    Runs ExomeDepth on a single chromosome for a set of BAM files.
    NB: each BAM file is assumed to contain a single sample. Multisample BAM files are not supported.

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run PhilPalmer/exomedepth --target_bed bedfile --chr 20 --ref reference.fasta -profile standard,docker

    Mandatory arguments:
      --target_bed                  BED file containing target regions to analyse (eg: exome design covered regions)
      --ref                         Genome reference in FASTA format, indexed using samtools faindex
      --bam_folder                  Path to folder containing bam files

    Options:
      -profile                      Available: standard, conda, docker, singularity, awsbatch, test
      --transition_probability      0.0001
      --msg                         Using Exome Depth transition probability //transition_probability
      --help                        Display help message

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.target_bed = 'data/genes_hg19.bed'
params.ref = 'data/test.fasta'
target_bed = file(params.target_bed)
ref = file(params.ref)
chrs = [1..22,'X']
//chrs = [ 'chr1', 'chr2', 'chrM' ]
params.transition_probability = 0.0001
transition_probability = params.transition_probability


// Validate inputs
if ( params.target_bed ){
    target_bed = file(params.target_bed)
    if( !target_bed.exists() ) exit 1, "Target bed file not found: ${params.target_bed}"
}



/*--------------------------------------------------
  Bam input files
  The input must be a path to a folder containing multiple bam files
---------------------------------------------------*/
params.bam_folder="data"
params.bam_file_prefix="*"


Channel.fromPath("${params.bam_folder}/${params.bam_file_prefix}*.bam").map{ file -> tuple(file.name, file) }.set{bamChannel}


// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/exomedepth
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/exomedepth'
summary['Target bed']    = params.target_bed
summary['Fasta Ref']    = params.ref
summary['Bams']    = "${params.bam_folder}/${params.bam_file_prefix}*.bam"
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="
Channel
/*
 * STEP 1 - Produce ExomeDepth tsv file for each chromosome
 */
process exomedepth{
    publishDir "${params.outdir}/ExomeDepth", mode: 'copy'

    input:
    file target_bed from target_bed
    set val(name), file(bam) from bamChannel
    each chr from chrs
    file ref

    output:
    //file "exome_depth.${chr}.tsv" into records
    set val(chr), file("exome_depth${chr}.tsv") into records

    script:
    """
    #!/usr/bin/env Rscript


    source("https://bioconductor.org/biocLite.R")
    biocLite("Rsamtools")

    library(ExomeDepth)
    library(Rsamtools)
    read.bed = function(f) {
      read.table(f, col.names=c("chr","start","end","id"), fill=1)
    }
    # Reference sequence
    hg19.fasta = "$ref"
    # Read the target / covered region
    print("Reading target regions for $chr from $target_bed")
    chr.covered = read.bed(pipe("grep '$chr[^0-9]' $target_bed"))
    # ExomeDepth wants the columns named differently
    chr.covered = data.frame(
        chromosome=chr.covered\$chr,
        start=chr.covered\$start,
        end=chr.covered\$end,
        name=paste(chr.covered\$chr,chr.covered\$start,chr.covered\$end,sep="-")
    )
    # BAM files are the primary input - convert to R vector
    bam.files = c('${bam}')
    all.samples = sapply(bam.files, function(file.name) {
        # Scan the bam header and parse out the sample name from the first read group
        read.group.info = strsplit(scanBamHeader(file.name)[[1]]\$text[["@RG"]],":")
        names(read.group.info) = sapply(read.group.info, function(field) field[[1]])
        return(read.group.info\$SM[[2]])
    })
    print(sprintf("Processing %d samples",length(all.samples)))
    # Finally we can call ExomeDepth
    bam.counts <- getBamCounts(bed.frame = chr.covered,
                              bam.files = bam.files,
                              include.chr = F,
                              referenceFasta = hg19.fasta)
    print("Successfully counted reads in BAM files")
    # Note: at this point bam.counts has column names reflecting the
    # file names => convert to actual sample names which is more convenient
    colnames(bam.counts) = c("GC", all.samples)
    # Problem: sample names starting with numbers get mangled. So convert them back,
    # but ignore the first column which is actually the GC percentage
    all.samples = colnames(bam.counts)[-1]
    for(test.sample in all.samples) {
        print(sprintf("Processing sample %s", test.sample))
        reference.samples = all.samples[-match(test.sample, all.samples)]
        bam.counts.df = as.data.frame(bam.counts[,reference.samples])[,-1:-6]
        test.sample.counts = bam.counts[,test.sample][[1]]
        print(sprintf("Selecting reference set for %s ...", test.sample ))
        reference = select.reference.set(
                                 test.counts = bam.counts[,test.sample][[1]],
                                 reference.counts = as.matrix(as.data.frame(bam.counts[,reference.samples])[,-1:-6]),
                                 bin.length = chr.covered\$end - chr.covered\$start
                                )
        # Get counts just for the reference set
        reference.counts = apply(bam.counts.df[,reference\$reference.choice,drop=F],1,sum)
        print(sprintf("Creating ExomeDepth object ..."))
        ed = new("ExomeDepth",
                      test = test.sample.counts,
                      reference = reference.counts,
                      formula = "cbind(test, reference) ~ 1")
        print(sprintf("Calling CNVs ..."))
        found.cnvs = CallCNVs(x = ed,
                                transition.probability = $transition_probability,
                                chromosome = chr.covered\$chromosome,
                                start = chr.covered\$start,
                                end = chr.covered\$end,
                                name = chr.covered\$name)
        results = found.cnvs@CNV.calls
        results\$sample = rep(test.sample, nrow(results))
        print(sprintf("Writing %d results to exome_depth.${chr}.tsv ...", nrow(results)))
        write.table(file="exome_depth.${chr}.tsv",
                    x=results,
                    row.names=F,
                    append=T)
    }
    print("Finished")
    """
}

/*
 * STEP 2 - Merge results from multiple ExomeDepth analyses together
 */
process merge {
    publishDir "${params.outdir}/Merge", mode: 'copy'

    input:
    //file (chrdepth) from records
    set val(chr), file(rec) from records

    output:
    file ('exome_depth.cnvs.tsv')

    script:
    """
    cat $rec | grep -v '^"sample"' | awk '{ if(NR==1 || \$1 != "\\"start.p\\"") print \$0 }' > $exome_depth.cnvs.tsv
    """
}
