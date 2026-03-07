### Boas Pucker ###
### pucker@uni-bonn.de ###
__version__ = "v0.1"

__usage__ = """
						PBB plasmid sequencing workflow (""" + __version__ +""")
						python pbb_plasmid_seq.py
						--reads <READS_FASTQ>
						--ref <REFERENCE_FASTA>
						--out <OUTPUT_DIRECTORY>
						--tmp <TMP_FOLDER>
						
						optional:
						--minimap <PATH_TO_MINIMAP2>
						--samtools <PATH_TO_SAMTOOLS>
						--threads <NUMBER_THREADS>[8]
						
						bug reports and feature requests: pucker@uni-bonn.de
					"""

import os, sys, subprocess, gzip

# --- end of imports --- #

def clean_header( header ):
	"""! @brief clean given header """
	
	illegal_character = [ " ", ":", ".", ";", ",", "(", ")", "|" ]
	for each in illegal_character:
		header = header.replace( each, "_" )
	return header


def clean_input_fasta_file( reference_file, clean_fasta_file ):
	"""! @brief clean FASTA file """
	
	sequences = {}
	with open( reference_file, "r" ) as f:
		header = clean_header( f.readline().strip()[1:] )
		line = f.readline()
		seq = []
		while line:
			if line[0] == ">":
				sequences.update( { header: "".join( seq ) } )
				seq = []
				header = clean_header( line.strip()[1:] )
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )
	
	with open( clean_fasta_file, "w" ) as out:
		for key in list( sequences.keys() ):
			out.write( '>' + key + "\n" + sequences[ key ] + "\n" )
	return sequences


def count_bases_in_FASTQ( fastq_file, qual_status=False ):
	"""! @brief count number of bases in FASTQ file """
	
	with gzip.open( fastq_file, "rb" ) as f:
		total_length = []
		average_quality = []
		line = f.readline().decode("utf-8")	#header
		while line:
			seq = f.readline().decode("utf-8").strip().upper()
			total_length.append( len( seq ) )
			f.readline()	#useless line
			qual = f.readline().decode("utf-8")	#quality line
			if qual_status:
				average_quality.append( calc_avg_qual( qual ) )
			line = f.readline().decode("utf-8")
		total_len = sum( total_length )
	return total_len


def gfa_to_fasta( gfa_file, assembly_file ):
	"""! @brief convert GFA to FASTA """
	
	with open( gfa_file ) as gfa, open( assembly_file, "w") as fasta:
		for line in gfa:
			if line.startswith("S"):
				parts = line.strip().split("\t")
				name = parts[1]
				seq = parts[2]
				fasta.write(f">{name}\n{seq}\n")


def main( arguments ):
	"""! @brief runs everything """
	
	reads_file = arguments[ arguments.index('--reads')+1 ]
	reference_file = arguments[ arguments.index('--ref')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	tmp_folder = arguments[ arguments.index('--tmp')+1 ]
	
	if '--minimap' in arguments:
		minimap = arguments[ arguments.index('--minimap')+1 ]
	else:
		minimap = "minimap2"
	
	if '--samtools' in arguments:
		samtools = arguments[ arguments.index('--samtools')+1 ]
	else:
		samtools = "samtools"
	
	if '--flye' in arguments:
		flye = arguments[ arguments.index('--flye')+1 ]
	else:
		flye = "flye"
	
	seqkit = "seqkit"
	seqtk = "seqtk"
	miniasm = "miniasm"
	racon = "racon"
	
	if '--threads' in arguments:
		threads = arguments[ arguments.index('--threads')+1 ]
	else:
		threads = "8"
	
	assembly_status = True
	if '--assembly' in arguments:
		if arguments[ arguments.index('--assembly')+1 ] == "off":
			assembly_status = False
	
	if output_folder[-1] != "/":
		output_folder += "/"
	
	if tmp_folder[-1] != "/":
		tmp_folder += "/"
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	if not os.path.exists( tmp_folder ):
		os.makedirs( tmp_folder )
	
	sys.stdout.write( "PBB plasmid sequencing: "+ __version__ + "\n" )
	sys.stdout.write( "reads: "+ reads_file + "\n" )
	sys.stdout.write( "reference: "+ reference_file + "\n" )
	sys.stdout.flush()
	
	# --- clean input FASTA file --- #
	clean_fasta_file = tmp_folder + "clean_reference_sequences.fasta"
	ref_seqs = clean_input_fasta_file( reference_file, clean_fasta_file )
	sys.stdout.write( "reference cleaning done\n" )
	sys.stdout.flush()
	
	# --- run minimap --- #
	sam_file = tmp_folder + "mapping.sam"
	doc_file =  tmp_folder + "mapping.errors.txt"
	if not os.path.isfile( sam_file ):	#skip if SAM already exists
		p = subprocess.Popen( args= " ".join( [ 	minimap,
																		"-ax map-ont -t", threads,
																		"--secondary=no",
																		clean_fasta_file,
																		reads_file,
																		">", sam_file,
																		"2>", doc_file+".err" 
																		] ), shell=True )
		p.communicate()
	sys.stdout.write( "read mapping done\n" )
	sys.stdout.flush()
	
	# --- convert SAM to BAM --- #
	bam_file = tmp_folder + "mapping.bam"
	doc_file2 =  tmp_folder + "samtools.conversion_errors.txt"
	if not os.path.isfile( bam_file ):	#skip if BAM already exists
		p = subprocess.Popen( args= " ".join( [ 	samtools, "view -bS", sam_file, ">", bam_file, "2>", doc_file2 ] ), shell=True )
		p.communicate()
	sys.stdout.write( "SAM to BAM conversion done\n" )
	sys.stdout.flush()
	
	# --- sort BAM --- #
	sorted_bam_file = tmp_folder + "mapping.sorted.bam"
	doc_file3 =  tmp_folder + "samtools.sorting_errors.txt"
	if not os.path.isfile( sorted_bam_file ):	#skip if sorted BAM already exists
		p = subprocess.Popen( args= " ".join( [ 	samtools, "sort", bam_file, "-o", sorted_bam_file, "2>", doc_file3 ] ), shell=True )
		p.communicate()
	sys.stdout.write( "BAM sorting done\n" )
	sys.stdout.flush()
	
	# --- index BAM --- #
	doc_file4 =  tmp_folder + "samtools.indexing_errors.txt"
	p = subprocess.Popen( args= " ".join( [ 	samtools, "index", sorted_bam_file, "2>", doc_file4 ] ), shell=True )
	p.communicate()
	sys.stdout.write( "BAM indexing done\n" )
	sys.stdout.flush()
	
	# --- split BAM by reference sequence into multiple BAMs (one per plasmid) --- #
	doc_file5 =  tmp_folder + "bam_file_splitting_errors.txt"
	for idx, seq in enumerate( list( ref_seqs.keys() ) ):
		
		sys.stdout.write( "Starting to process individual plasmids: " + seq + "\t(" + str( idx+1 ) +"/"+ str( len( list( ref_seqs.keys() ) ) ) + ")\n" )
		sys.stdout.flush()
		
		# generate BAM file
		bam_per_ref = output_folder + seq + ".bam"
		if not os.path.isfile( bam_per_ref ):
			p = subprocess.Popen( args= " ".join( [ 	samtools, "view -b", sorted_bam_file, seq, ">", bam_per_ref, "2>>", doc_file5 ] ), shell=True )
			p.communicate()
			sys.stdout.write( "BAM done\n" )
			sys.stdout.flush()
		
			# index file
			p = subprocess.Popen( args= " ".join( [ 	samtools, "index", bam_per_ref, "2>>", doc_file5 ] ), shell=True )
			p.communicate()
		
		#generate FASTA file
		fasta_file = output_folder + seq + ".fasta"
		with open( fasta_file, "w" ) as out:
			out.write( '>' + seq + "\n" + ref_seqs[ seq ] + "\n" )
		p = subprocess.Popen( args= " ".join( [ 	samtools, "faidx", fasta_file, "2>>", doc_file5 ] ), shell=True )
		p.communicate()
		sys.stdout.write( "FASTA done\n" )
		sys.stdout.flush()
		
		#extract mapped reads into FASTQ file
		fastq_file = tmp_folder + seq + ".fastq.gz"
		p = subprocess.Popen( args= " ".join( [ samtools, "fastq -F 2308", bam_per_ref, "|", seqkit, "rmdup -n | gzip >", fastq_file, "2>>", doc_file5 ] ), shell=True )
		#-F 2308: unmapped (4); secondary (256); supplementary (2048)
		p.communicate()
		sys.stdout.write( "FASTQ done\n" )
		sys.stdout.flush()
		
		#calculate coverage and calculate read reduction factor
		genome_size = len( ref_seqs[ seq ] )
		target_coverage = 200
		total_bases = count_bases_in_FASTQ( fastq_file )
		if total_bases == 0:	#avoid ZeroDivisionError
			total_bases += 1
		current_coverage = total_bases / genome_size
		ratio_to_keep = target_coverage / float( current_coverage )
		
		sys.stdout.write( "Genome size:" + str( genome_size ) + "\n" )
		sys.stdout.write( "Total bases:" + str( total_bases ) + "\n" )
		sys.stdout.write( "Current coverage:" + str( current_coverage ) + "\n" )
		sys.stdout.write( "Ratio:" + str( ratio_to_keep ) + "\n" )
		sys.stdout.flush()
		
		
		#randomly reduce number of reads
		subset_fastq_file = output_folder + seq + ".subset.fastq.gz"
		p = subprocess.Popen( args= " ".join( [ 	seqtk, "sample -s42", fastq_file, str( ratio_to_keep ), "| gzip >", subset_fastq_file, "2>>", doc_file5 ] ), shell=True )
		#randomly take X% of reads
		p.communicate()
		sys.stdout.write( "FASTQ filtering done\n" )
		sys.stdout.flush()		
		
		#assemble plasmid sequence with flye
		if assembly_status:	#would allow to switch this off
			
			miniasm_result_folder = tmp_folder + seq + ".miniasm/"
			if not os.path.exists( miniasm_result_folder ):
				os.makedirs( miniasm_result_folder )
				paf_file = miniasm_result_folder + "overlaps.paf"
				p = subprocess.Popen( args= " ".join( [ 	minimap,
																				"-x ava-ont",
																				subset_fastq_file,
																				subset_fastq_file,
																				">",
																				paf_file,
																				"2>>", doc_file5 
																			] ), shell=True )
				p.communicate()
				sys.stdout.write( "Minimap for miniasm done\n" )
				sys.stdout.flush()
				
				gfa_file = miniasm_result_folder + "assembly.gfa"
				p = subprocess.Popen( args= " ".join( [ 	miniasm, "-f", subset_fastq_file, paf_file, ">", gfa_file, "2>>", doc_file5 ] ), shell=True )
				p.communicate()
				sys.stdout.write( "Initial miniasm assembly done\n" )
				sys.stdout.flush()
				
				assembly_file = miniasm_result_folder + "assembly.fasta"
				gfa_to_fasta( gfa_file, assembly_file )
				sys.stdout.write( "Contig extraction done\n" )
				sys.stdout.flush()			
				
				assembly_paf_file = miniasm_result_folder + "assembly.paf"
				p = subprocess.Popen( args= " ".join( [ 	minimap,
																				"-x map-ont",
																				assembly_file,
																				subset_fastq_file,
																				">",
																				assembly_paf_file,
																				"2>>", doc_file5
																			] ), shell=True )
				p.communicate()
				sys.stdout.write( "Mapping to assembly done\n" )
				sys.stdout.flush()
				
				polished_assembly_file = miniasm_result_folder + "polished_assembly.fasta"
				p = subprocess.Popen( args= " ".join( [ racon, subset_fastq_file, assembly_paf_file, assembly_file, ">", polished_assembly_file, "2>>", doc_file5 ] ), shell=True )
				p.communicate()
				sys.stdout.write( "Assembly polishing done\n" )
				sys.stdout.flush()
				
				#copy assembly from that output folder				
				p = subprocess.Popen( args= " ".join( [ "cp", polished_assembly_file, output_folder + seq + ".assembly.fasta", "2>>", doc_file5 ] ), shell=True )
				p.communicate()
				
				sys.stdout.write( "Miniasm assembly process and correction done\n" )
				sys.stdout.flush()


if '--reads' in sys.argv and '--ref' in sys.argv and '--out' in sys.argv and '--tmp' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
	
