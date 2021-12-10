process download_fastq {
    output:
        file '*.sra' into get_sample_channel

    script:
    """
    #!/usr/bin/env bash

    wget -O SRR628582.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR628582/SRR628582.1 ;
    """
}

process split_fastq {
    input:
	file '*' from get_sample_channel

    output:
        file '*.gz' into reads_channel

    script:
    """
    #!/usr/bin/env bash

    fastq-dump --gzip --split-files ./SRR628582.sra
    """
}

process download_human_genome {
    output:
        file '*.fa' into human_genome_channel

    script:
    """
    #!/usr/bin/env bash

    wget ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ;
    gunzip *.fa.gz
    """
}

process download_human_annotation {
    output:
        file '*.gtf' into human_annotation_channel


    script:
    """
    #!/usr/bin/env bash

    wget ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz
    gunzip Homo_sapiens.GRCh38.104.chr.gtf.gz
    """
}


process create_index {
    input:
	file fasta from human_genome_channel
	file annotation from human_annotation_channel

    output:
	path 'GenomeDir' into index_path

    script:
    """
    #!/usr/bin/env bash

    STAR --runThreadN 12 \
      --runMode genomeGenerate \
      --sjdbGTFfile $annotation \
      --genomeFastaFiles $fasta
    """
}

process mapping {
    input:
	path GenomeDir from index_path
	file fastq from reads_channel

    output:
    	file "*.bam" into bam_files

    script:
    """
    #!/usr/bin/env bash

    STAR --outSAMstrandField intronMotif \
      --outFilterMismatchNmax 4 \
      --outFilterMultimapNmax 10 \
      --readFilesIn <(gunzip ${fastq[0]}) <(gunzip ${fastq[1]}) \
      --runThreadN 12 \
      --genomeDir GenomeDir \
      --outSAMunmapped None \
      --outSAMtype BAM SortedByCoordinate \
      --outStd BAM_SortedByCoordinate \
      --genomeLoad NoSharedMemory \
      --limitBAMsortRAM '50GB' \
      > ${fastq[0]}.bam
    """
}

bam_files.view()
