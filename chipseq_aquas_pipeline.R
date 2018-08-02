suppressMessages(library(GenomicFeatures))
suppressMessages(library(GenomicAlignments))

# TODO: remove absolute paths and turn them into inputs
scripts.dir <- "/storage/cylin/bin/aquas/chip"
fn.dir <- "/storage/cylin/anaconda3/envs/aquas_chipseq/bin/"

##This functions runs aquas up until alignment and cross-correlation
pushAquasAlignQC <- function(project.dir, dt.location, chip.type = "histone", 
	no.pseudo.rep = T, ncores = 32){

	df.sample.info <- read.delim(dt.location, sep = "\t", header = T)
	sample.names <- as.character(df.sample.info$NAME)
	sapply(sample.names, function(sample.name){

		##Get the genome
		genome.name <- as.character(subset(df.sample.info, NAME == sample.name)$GENOME)

		##location of the fastq file
		fastq.file.location <- as.character(subset(df.sample.info, NAME == sample.name)$FASTQ_FILE)
		fastq.files <- unlist(strsplit(split = "::", fastq.file.location))
		##Find the endedness of the fastqs
		if(length(fastq.files) == 2){
			endedness <- "-pe"
		} else if (length(fastq.files) == 1){
			endedness <- "-se"
		}
		##Check if the files exists
		if(!all(file.exists(fastq.files)))
		stop("Error: Fastq files missing")

		##location of the bam file
		file.id <- as.character(subset(df.sample.info, NAME == sample.name)$UNIQUE_ID)
		bam.file.dir <- as.character(subset(df.sample.info, UNIQUE_ID == file.id)$FILE_PATH)
  		all.bam.files <- list.files(path = bam.file.dir, pattern = "*.bam$")
		bam.file.logical <- grepl(pattern = file.id, x = all.bam.files)
		##Check if the bam file already exists
		if(all(!bam.file.logical)){
			aligned <- FALSE
		} else {
			bam.file.location <- file.path(bam.file.dir, all.bam.files[bam.file.logical])
			if(file.exists(bam.file.location)){
				print(sprintf("BAM already exists for %s", sample.name))
			}
			aligned <- TRUE
		}

		aquas.dir <- file.path(project.dir, "aquas_output")
		if(!dir.exists(aquas.dir)){
 			dir.create(aquas.dir)
 		}

		##setting parameters
		output.dir <- file.path(aquas.dir, sample.name)
		if(aligned){
				input.file.cmd <- sprintf("-bam1 %s", bam.file.location)
		} else if(endedness == "-se"){
			input.file.cmd <- sprintf("-fastq1 %s", fastq.files[1])
			} else {
				input.file.cmd <- sprintf("-fastq1_1 %s -fastq1_2 %s", fastq.files[1], fastq.files[2])
			}

		##Check if everything needs to be done for true replicates
		if(no.pseudo.rep){
			input.file.cmd <- sprintf("%s -true_rep", input.file.cmd)
		}

		##Prepare the final aquas command
		screen.name <- sprintf("align%d%d", sample(1:1000000, size = 1), sample(1:1000000, size = 1))
		aquas.cmd <- sprintf("python %s/chipseq.py -screen %s", scripts.dir, screen.name)
		aquas.cmd <- sprintf("%s -type %s -final_stage xcor", aquas.cmd, chip.type)
		if(endedness == "-pe"){
			aquas.cmd <- sprintf("%s -pe_no_trim_fastq", aquas.cmd)
		}
		aquas.cmd <- sprintf("%s -out_dir %s %s -species %s", aquas.cmd, output.dir, endedness, genome.name)
		aquas.cmd <- sprintf("%s -system slurm -nth %d %s", aquas.cmd, ncores, input.file.cmd)

		#open the bashfile to write to
    	bash.file.location <- file.path(aquas.dir, sprintf("%s_align_aquas.sh", sample.name))
		file.conn <- file(bash.file.location)
    	writeLines(c('#!/usr/bin/bash', aquas.cmd), file.conn)
		close(file.conn)

		##Execute
		system(sprintf("chmod u+x %s", bash.file.location))
		system(bash.file.location)
		})
}

##This functions uses aquas to call peaks
callAquasPeaks <- function(project.dir, dt.location, 
	no.pseudo.rep = T, peak.call.threshold = 0.00001, ncores = 32){

	aquas.dir <- file.path(project.dir, "aquas_output")
	if(!dir.exists(aquas.dir)){
 		dir.create(aquas.dir)
 	}

	df.sample.info <- read.delim(dt.location, sep = "\t", header = T)
	sample.names <- as.character(subset(df.sample.info, ENRICHED_MACS != "NONE")$NAME)
	sapply(sample.names, function(sample.name){

		##Get the genome
		genome.name <- as.character(subset(df.sample.info, NAME == sample.name)$GENOME)

		##location of the fastq file
		chip.type <- as.character(subset(df.sample.info, NAME == sample.name)$TYPE)
		if(!(chip.type %in% c("TF", "histone"))){
			stop("ChIP type not included")
		}

		##location of the fastq file
		fastq.file.location <- as.character(subset(df.sample.info, NAME == sample.name)$FASTQ_FILE)
		fastq.files <- unlist(strsplit(split = "::", fastq.file.location))
		##Find the endedness of the fastqs
		if(length(fastq.files) == 2){
			endedness <- "-pe"
		} else if (length(fastq.files) == 1){
			endedness <- "-se"
		}

		##location of the tag align file for experiment
		expt.file.id <- as.character(subset(df.sample.info, NAME == sample.name)$UNIQUE_ID)
		expt.bam.file.dir <- as.character(subset(df.sample.info, UNIQUE_ID == expt.file.id)$FILE_PATH)
  		expt.all.bam.files <- list.files(path = expt.bam.file.dir, pattern = "*.nodup.tagAlign.gz$", full.names = T)
		expt.bam.file.logical <- grepl(pattern = expt.file.id, x = expt.all.bam.files)
		##Check if the tag align file already exists
		if(all(!expt.bam.file.logical)){
			expt.aligned <- FALSE
		} else {
			expt.bam.file.location <- expt.all.bam.files[expt.bam.file.logical]
			expt.aligned <- TRUE
		}

		if(!expt.aligned){
			stop("Experiment bam file not found")
		}

		##setting parameters
		output.dir <- file.path(aquas.dir, sample.name)
		input.file.cmd <- sprintf("-tag1 %s", expt.bam.file.location)

		bkg.sample.name <- as.character(subset(df.sample.info, NAME == sample.name)$BACKGROUND)
		if(bkg.sample.name != "NONE"){
			bkg.file.id <- as.character(subset(df.sample.info, NAME == bkg.sample.name)$UNIQUE_ID)
			bkg.bam.file.dir <- as.character(subset(df.sample.info, UNIQUE_ID == bkg.file.id)$FILE_PATH)
  			bkg.all.bam.files <- list.files(path = bkg.bam.file.dir, pattern = "*.nodup.tagAlign.gz$", full.names = T)
			bkg.bam.file.logical <- grepl(pattern = bkg.file.id, x = bkg.all.bam.files)
			##Check if the bam file already exists
			if(all(!bkg.bam.file.logical)){
				bkg.aligned <- FALSE
			} else {
				bkg.bam.file.location <- bkg.all.bam.files[bkg.bam.file.logical]
				bkg.aligned <- TRUE
			}
			if(!bkg.aligned){
				stop("background bam file not found")
			}
			input.file.cmd <- sprintf("%s -ctl_tag %s", input.file.cmd, bkg.bam.file.location)
		}

		##Check if everything needs to be done for true replicates
		if(no.pseudo.rep){
			input.file.cmd <- sprintf("%s -true_rep", input.file.cmd)
		}

		##Prepare the final aquas command
		screen.name <- sprintf("peakcall%d%d", sample(1:1000000, size = 1), sample(1:1000000, size = 1))
		aquas.cmd <- sprintf("python %s/chipseq.py -screen %s", scripts.dir, screen.name)
		aquas.cmd <- sprintf("%s -type %s", aquas.cmd, chip.type)
		aquas.cmd <- sprintf("%s -out_dir %s %s -species %s", aquas.cmd, output.dir, endedness, genome.name)
		aquas.cmd <- sprintf("%s -system slurm -nth %d %s", aquas.cmd, ncores, input.file.cmd)
		aquas.cmd <- sprintf("%s -peak_caller macs2 -pval_thresh_macs2 %f", aquas.cmd, peak.call.threshold)

		print(aquas.cmd)

		#open the bashfile to write to
    	bash.file.location <- file.path(aquas.dir, sprintf("%s_peaks_aquas.sh", sample.name))
		file.conn <- file(bash.file.location)
    	writeLines(c('#!/usr/bin/bash', aquas.cmd), file.conn)
		close(file.conn)

		##Execute
		system(sprintf("chmod u+x %s", bash.file.location))
		system(bash.file.location)
		})
}

##This function moves all the qc data to qc_summary directory
##and writes all the qc information pre_peak_call_qc.txt
collectAquasQC <- function(project.dir){

	aquas.dir <- file.path(project.dir, "aquas_output")
	if(!dir.exists(aquas.dir)){
		stop("No aquas output")
	}

	qc.dir <- file.path(project.dir, "qc_summary")
	if(!dir.exists(qc.dir)){
 		dir.create(qc.dir)
 	}

 	##Get all the qc summary json files
 	summary.json.files <- list.files(path = aquas.dir, pattern = "ENCODE_summary.json", full.names = T,
 		recursive = T)
 	new.summary.json.files <- file.path(qc.dir, paste(basename(dirname(summary.json.files)), "_aquas_summary.json", sep = ""))

 	##copy them over
 	for(idx in 1:length(summary.json.files)){
 		cp.cmd <- sprintf("cp %s %s", summary.json.files[idx], new.summary.json.files[idx])
 		system(cp.cmd)
 	}

 	combined.output.file <- file.path(qc.dir, "pre_peak_call_qc.txt")
	combine.cmd <- sprintf("%s/utils/parse_summary_qc_recursively.py --search-dir %s ", scripts.dir, qc.dir)
	combine.cmd <- sprintf("%s --json-file *_aquas_summary.json > %s", combine.cmd, combined.output.file)
	system(combine.cmd)

}

##Moves all the aquas generated bams and tag aligns to appropriate location
mvBams <- function(project.dir, dt.location){

	##Check for presence of necessary directories
	aquas.dir <- file.path(project.dir, "aquas_output")
	if(!dir.exists(aquas.dir)){
		stop("No aquas output")
	}

	bam.dir <- file.path(project.dir, "bams")
	if(!dir.exists(bam.dir)){
 		dir.create(bam.dir)
 	}

 	df.sample.info <- read.delim(dt.location, sep = "\t", header = T)
	sample.names <- as.character(df.sample.info$NAME)

	##Move the bams to appropriate location
	sapply(sample.names, function(sample.name) {
		
		file.id <- as.character(subset(df.sample.info, NAME == sample.name)$UNIQUE_ID)

		##location of the fastq file
		fastq.file.location <- as.character(subset(df.sample.info, NAME == sample.name)$FASTQ_FILE)
		fastq.files <- unlist(strsplit(split = "::", fastq.file.location))

		##Move the bam files
		if(length(fastq.files) == 2){
			all.bam.files <- list.files(path = aquas.dir, pattern = "*.PE2SE.nodup.bam$", 
  				full.names = T, recursive = T)
		} else if (length(fastq.files) == 1){
			all.bam.files <- list.files(path = aquas.dir, pattern = "*.nodup.bam$", 
  				full.names = T, recursive = T)
		}
  		
		bam.file.logical <- grepl(pattern = file.id, x = all.bam.files)
		bam.file.location <- file.path(all.bam.files[bam.file.logical])

		if(length(bam.file.location) == 1){

			new.bam.file.location <- file.path(bam.dir, basename(bam.file.location))
 			##Move the bam files
 			mv.bam.cmd <- sprintf("mv %s %s", bam.file.location, new.bam.file.location)
 			system(mv.bam.cmd)

 			##Move the bai files
 			mv.bai.cmd <- sprintf("mv %s %s", sprintf("%s.bai", bam.file.location), sprintf("%s.bai", new.bam.file.location))
 			system(mv.bai.cmd)
		} else {
			return(sample.name)
		}
	})

	##Move the tagalign files to appropriate location and find the samples
	##whoes tag align files are missing
	missing.samples <- sapply(sample.names, function(sample.name) {
		
		file.id <- as.character(subset(df.sample.info, NAME == sample.name)$UNIQUE_ID)

		##location of the fastq file
		fastq.file.location <- as.character(subset(df.sample.info, NAME == sample.name)$FASTQ_FILE)
		fastq.files <- unlist(strsplit(split = "::", fastq.file.location))

		##Move the tagAlign files
		if(length(fastq.files) == 2){
			all.tagAlign.files  <- list.files(path = aquas.dir, pattern = "*.PE2SE.nodup.tagAlign.gz$",
  				full.names = T, recursive = T)
		} else if (length(fastq.files) == 1){
			all.tagAlign.files <- list.files(path = aquas.dir, pattern = "*.nodup.tagAlign.gz$",
  				full.names = T, recursive = T)
		}

		tagAlign.file.logical <- grepl(pattern = file.id, x = all.tagAlign.files)
		tagAlign.file.location <- file.path(all.tagAlign.files[tagAlign.file.logical])

		if(length(tagAlign.file.location) == 1){

			new.tagAlign.file.location <- file.path(bam.dir, basename(tagAlign.file.location))
 			##Move the tagAlign files
 			mv.tagAlign.cmd <- sprintf("mv %s %s", tagAlign.file.location, new.tagAlign.file.location)
 			system(mv.tagAlign.cmd)
 			return(NULL)
		} else {
			return(sample.name)
		}
	})

	return(unlist(missing.samples))
}

##Moves all the aquas generated peak files to appropriate location
mvPeakFiles <- function(project.dir, dt.location){

	aquas.dir <- file.path(project.dir, "aquas_output")
	if(!dir.exists(aquas.dir)){
		stop("No aquas output")
	}

	peak.dir <- file.path(project.dir, "macsEnriched")
	if(!dir.exists(peak.dir)){
 		dir.create(peak.dir)
 	}

	df.sample.info <- read.delim(dt.location, sep = "\t", header = T)
	sample.names <- as.character(subset(df.sample.info, ENRICHED_MACS != "NONE")$NAME)

	missing.samples <- sapply(sample.names, function(sample.name) {
		print(sample.name)
		file.id <- as.character(subset(df.sample.info, NAME == sample.name)$UNIQUE_ID)

		##Move the peak file
		all.peak.files <- list.files(path = aquas.dir, pattern = "*.narrowPeak.gz$", full.names = T, recursive = T)
		peak.file.logical <- grepl(pattern = file.id, x = all.peak.files)
		peak.file.location <- file.path(all.peak.files[peak.file.logical])
		if(length(peak.file.location) > 1){
			stop("Error: More than 1 peak file")
		}
		if(length(peak.file.location) == 1){
			new.peak.file.name <- as.character(subset(df.sample.info, NAME == sample.name)$ENRICHED_MACS)
			new.peak.file.location <- file.path(peak.dir, sprintf("%s.gz", new.peak.file.name))
			mv.peak.cmd <- sprintf("mv %s %s", peak.file.location, new.peak.file.location)
 			unzip.peak.cmd <- sprintf("gunzip %s", new.peak.file.location)
 			system(mv.peak.cmd)
 			system(unzip.peak.cmd)
 			return(NULL)
 			} else {
 				return(sample.name)
 			}
	})
	return(unlist(missing.samples))	
}

createSignalTrack <- function(project.dir, dt.location){

	aquas.dir <- file.path(project.dir, "aquas_output")
	if(!dir.exists(aquas.dir)){
		stop("No aquas output")
	}

	bw.dir <- file.path(project.dir, "wiggles")
	if(!dir.exists(bw.dir)){
 		dir.create(bw.dir)
 	}

	df.sample.info <- read.delim(dt.location, sep = "\t", header = T)
	sample.names <- as.character(subset(df.sample.info, ENRICHED_MACS != "NONE")$NAME)
	sapply(sample.names, function(sample.name){

		##Get the genome
		genome.name <- as.character(subset(df.sample.info, NAME == sample.name)$GENOME)

		##location of the bam file for experiment
		expt.file.id <- as.character(subset(df.sample.info, NAME == sample.name)$UNIQUE_ID)
		expt.bam.file.dir <- as.character(subset(df.sample.info, UNIQUE_ID == expt.file.id)$FILE_PATH)
  		expt.all.bam.files <- list.files(path = expt.bam.file.dir, pattern = "*.nodup.bam$", full.names = T)
		expt.bam.file.logical <- grepl(pattern = expt.file.id, x = expt.all.bam.files)
		##Check if the bam file already exists
		if(all(!expt.bam.file.logical)){
			expt.aligned <- FALSE
		} else {
			expt.bam.file.location <- expt.all.bam.files[expt.bam.file.logical]
			expt.aligned <- TRUE
		}

		if(!expt.aligned){
			stop("Experiment bam file not found")
		}

		bkg.sample.name <- as.character(subset(df.sample.info, NAME == sample.name)$BACKGROUND)
		if(bkg.sample.name != "NONE"){
			bkg.file.id <- as.character(subset(df.sample.info, NAME == bkg.sample.name)$UNIQUE_ID)
			bkg.bam.file.dir <- as.character(subset(df.sample.info, UNIQUE_ID == bkg.file.id)$FILE_PATH)
  			bkg.all.bam.files <- list.files(path = bkg.bam.file.dir, pattern = "*.nodup.bam$", full.names = T)
			bkg.bam.file.logical <- grepl(pattern = bkg.file.id, x = bkg.all.bam.files)
			##Check if the bam file already exists
			if(all(!bkg.bam.file.logical)){
				bkg.aligned <- FALSE
			} else {
				bkg.bam.file.location <- bkg.all.bam.files[bkg.bam.file.logical]
				bkg.aligned <- TRUE
			}
			if(!bkg.aligned){
				stop("Background bam file not found")
			}
		}

		##get estimated fragment length from cross-corr. analysis log
		cc.qc.log.files <- list.files(path = aquas.dir, pattern = "*.cc.qc$", full.names = T,
			recursive = T)
		expt.cc.file.logical <- grepl(pattern = expt.file.id, x = cc.qc.log.files)
		expt.cc.file.location <- cc.qc.log.files[expt.cc.file.logical]
		fraglen <- read.delim(expt.cc.file.location, header = F)[1, 3]

		macs2.location <- file.path(fn.dir, "macs2")
		bdgTobw.location <- file.path(fn.dir, "bedGraphToBigWig")
		bwTowig.location <- "/storage/cylin/home/irtishas/tools/ucsc/bigWigToWig"

		if(genome.name %in% c("hg19", "hg38")){
			genome.abbrv <- "hs"
		} else if (genome.name == c("mm10")){
			genome.abbrv <- "mm"
		}
		chrom.size.file <- "/storage/cylin/grail/genomes/chrom_sizes/"
		chrom.size.file <- sprintf("%s%s.chrom.sizes", chrom.size.file, genome.name)
		
		bdg.output.file <- file.path(bw.dir, sample.name)
		#bw.output.file <- sprintf("%s%s", bdg.output.file, "_treat_pileup.bw")
		#wig.output.file <- sprintf("%s%s", bdg.output.file, "_treat_pileup.wig")

		signal.cmd <- sprintf("%s callpeak -t %s -c %s -B --nomodel --extsize %d --trackline --SPMR -g %s -n %s",
			macs2.location, expt.bam.file.location, bkg.bam.file.location, fraglen, genome.abbrv, bdg.output.file)
		# TODO: make it work if there is no background
		track.name.cmd <- sprintf("sed -i '/name/ s/treatment pileup/%s/' %s_treat_pileup.bdg", 
			sample.name ,bdg.output.file)
		gzip.cmd <- sprintf("gzip %s_treat_pileup.bdg", bdg.output.file)
		rm.cmd <- sprintf("rm %s_peaks.xls %s_peaks.narrowPeak %s_summits.bed %s_control_lambda.bdg", 
			bdg.output.file, bdg.output.file, bdg.output.file, bdg.output.file)	
		#sort.cmd1 <- "export LC_COLLATE=C"
		#sort.cmd2 <- sprintf("sort -k 1,1 -k 2,2n %s_treat_pileup.bdg > %s_sorted_treat_pileup.bdg", 
		#	bdg.output.file, bdg.output.file)
		#mv.cmd <- sprintf("mv %s_sorted_treat_pileup.bdg %s_treat_pileup.bdg", bdg.output.file, bdg.output.file)
		#bdgTobw.cmd <- sprintf("%s %s_treat_pileup.bdg %s %s", 
		#	bdgTobw.location, bdg.output.file, chrom.size.file, bw.output.file)
		#bwTowig.cmd <- sprintf("%s %s %s", 
		#	bwTowig.location, bw.output.file, wig.output.file)
		#gzip.cmd <- sprintf("gzip %s", wig.output.file)
	
		#open the bashfile to write to
    	bash.file.location <- file.path(aquas.dir, sprintf("%s_signal_tracks.sh", sample.name))
		file.conn <- file(bash.file.location)
    	writeLines(c('#!/usr/bin/bash', signal.cmd, track.name.cmd, gzip.cmd, rm.cmd), file.conn)
		close(file.conn)

		system(sprintf("chmod u+x %s", bash.file.location))
		system(sprintf("sbatch %s", bash.file.location))

		})

}

##This function will read in the bam file for eveyr sample
##It will shift the reads by aqurate fragment length
##Then it will count the number of reads 
getCustomQC <- function(project.dir, dt.location, peaks.dirname){
	aquas.dir <- file.path(project.dir, "aquas_output")
	if(!dir.exists(aquas.dir)){
		stop("No aquas output")
	}
	df.sample.info <- read.delim(dt.location, sep = "\t", header = T)
	sample.names <- as.character(subset(df.sample.info, ENRICHED_MACS != "NONE")$NAME)
	df <- do.call(rbind, mclapply(sample.names, function(sample.name){
		print(sample.name)
		##Get the genome
		genome.name <- as.character(subset(df.sample.info, NAME == sample.name)$GENOME)

		##location of the bam file for experiment
		expt.file.id <- as.character(subset(df.sample.info, NAME == sample.name)$UNIQUE_ID)
		expt.bam.file.dir <- as.character(subset(df.sample.info, UNIQUE_ID == expt.file.id)$FILE_PATH)
  		expt.all.bam.files <- list.files(path = expt.bam.file.dir, pattern = "*.nodup.bam$", full.names = T)
		expt.bam.file.logical <- grepl(pattern = expt.file.id, x = expt.all.bam.files)
		##Check if the bam file already exists
		if(all(!expt.bam.file.logical)){
			expt.aligned <- FALSE
		} else {
			expt.bam.file.location <- expt.all.bam.files[expt.bam.file.logical]
			expt.aligned <- TRUE
		}

		if(!expt.aligned){
			stop("Experiment bam file not found")
		}

		##location of the fastq file
		fastq.file.location <- as.character(subset(df.sample.info, NAME == sample.name)$FASTQ_FILE)
		fastq.files <- unlist(strsplit(split = "::", fastq.file.location))

		##get estimated fragment length from cross-corr. analysis log
		if(length(fastq.files) == 2){
			cc.qc.log.files <- list.files(path = aquas.dir, pattern = "*.R1.cc.qc$", full.names = T,
				recursive = T)
		} else if (length(fastq.files) == 1){
			cc.qc.log.files <- list.files(path = aquas.dir, pattern = "*.cc.qc$", full.names = T,
				recursive = T)
		}
		expt.cc.file.logical <- grepl(pattern = expt.file.id, x = cc.qc.log.files)
		expt.cc.file.location <- cc.qc.log.files[expt.cc.file.logical]
		fraglen <- read.delim(expt.cc.file.location, header = F)[1, 3]
		half.fraglen <- round((fraglen + 1)/2, digit = 0)

		##Read the bam
		g1 <- readGAlignments(expt.bam.file.location)
  		gr.reads <- as(g1, "GRanges")
  		rm(g1)
  		Nreads <- length(gr.reads)

		##Shift the reads by the fragment lengh
		gr.pos.reads <- GenomicRanges::shift(gr.reads[which(strand(gr.reads) == "+")], half.fraglen)
		gr.neg.reads <- GenomicRanges::shift(gr.reads[which(strand(gr.reads) == "-")], -half.fraglen)
		gr.shifted.reads <- c(gr.pos.reads, gr.neg.reads)

		##Find how many reads overlap a promoter
		gr.shifted.reads <- annotatePromoterPeaks(gr.shifted.reads, genome = genome.name, 
			tss.flank = 1000)
  		ans <- table(gr.shifted.reads$annotation)
  		proportion <- round(ans/length(gr.shifted.reads), digit = 2) * 100
  		promoter.n <- ans["promoter"]
  		promoter <- proportion["promoter"]

		Npeaks <- 0
		frip <- 0
		promoter.peaks.n <- 0
		promoter.peaks <- 0

		##Find how many reads overlap peaks and how many peaks are withing 1KB of a TSS
		bed.file.name <- as.character(subset(df.sample.info, UNIQUE_ID == expt.file.id)$ENRICHED_MACS)
		bed.file.path <- file.path(project.dir, peaks.dirname, bed.file.name)
		nlines.str <- system(sprintf("wc -l %s", bed.file.path), intern = TRUE)
		nlines <- as.numeric(strsplit(split = " ", nlines.str)[[1]][1])
		if(nlines > 0){
			df.ranges <- read.delim(bed.file.path, sep = "\t", header = F)
				df.ranges <- df.ranges[ ,1:3]
				colnames(df.ranges) <- c("seqnames", "start", "end")
				gr.peaks <- as(df.ranges, "GRanges")
				Npeaks <- length(gr.peaks)
				cnt <- sum(countOverlaps(gr.peaks, gr.shifted.reads))
				frip <- round(cnt/Nreads, digit = 2)
				gr.peaks <- annotatePromoterPeaks(gr.peaks, genome = genome.name, 
				tss.flank = 1000)
  	  		ans.peaks <- table(gr.peaks$annotation)
  	  		proportion.peaks <- round(ans.peaks/Npeaks, digit = 2) * 100
  	  		promoter.peaks <- proportion.peaks["promoter"]
  	  		promoter.peaks.n <- ans.peaks["promoter"]
		}


		df.stats <- data.frame(sample = sample.name, frip.score = frip,
		number.of.peaks = Npeaks,
		number.reads.mapped.to.promoter = promoter.n,
		fraction.reads.mapped.to.promoter = promoter, 
		number.peaks.mapped.to.promoter = promoter.peaks.n,
		fraction.peaks.mapped.to.promoter = promoter.peaks)
		rownames(df.stats) <- NULL
		print(df.stats)
		return(df.stats)
	}, mc.cores = 10))

	return(df)

}

##Takes in peaks and checks if any of the peaks overlaps within 1KB of a a know TSS
annotatePromoterPeaks <- function(gr.peaks, genome = NULL, tss.flank = 1000){

	gr.peaks$.idx <- 1:length(gr.peaks)
	gr.peaks$annotation <- rep(NA_character_, length(gr.peaks))

	if(is.null(genome)){
	  stop("genome name absent")
	}

	annotation.dir <- ("/storage/cylin/grail/annotations/txdb_sqlite")
	genome.db.name <- sprintf("%s.refGene.sqlite", genome)
	genome.db.location <- file.path(annotation.dir, genome.db.name)
	
	if(!file.exists(genome.db.location)){
	  stop("annotation db does not exist")
	}

	txdb <- loadDb(genome.db.location)
	##Get the transcripts by gene
	ls.gr.tx <- transcriptsBy(txdb, by = "gene")
	gids <- rep(names(ls.gr.tx), elementNROWS(ls.gr.tx))
	df.tx <- transform(as.data.frame(unname(unlist(ls.gr.tx))), gene_id=gids)
	gr.tx <- as(df.tx, "GRanges")

	##Annotate the nearest TSS
	gr.tss <- resize(gr.tx, width = 1L, fix = "start")
	gr.tss <- resize(gr.tss, width = (2*tss.flank), fix = "center")
	qt <- queryHits(findOverlaps(gr.peaks, gr.tss, ignore.strand = T))
	pro.annotation <- rep(NA_character_, length(gr.peaks))
	pro.annotation[qt] <- "promoter"
	gr.peaks$annotation <- pro.annotation
	return(gr.peaks)

}
