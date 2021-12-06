process download_fastq {
    output:
        file *.gz into reads_channel

    script:
    """
    #!/usr/bin/env bash

    SRAID='SRR628582' ; #other samples : SRR628583,SRR628584,SRR628585,SRR628586,SRR628587,SRR628588,SRR628589,SRR628589
    wget -O $SRAID.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/$SRAID/$SRAID.1 ;
    fastq-dump --gzip --split-files ./$SRAID.sra
    """
}

process download_human_genome {
    output:
        *fa file into human_genome_channel

    script:
    """
    #!/usr/bin/env bash

    wget ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ;
    gunzip *fa.gz > *fa
    """
}

process download_human_annotation {
    output:
        file Homo_sapiens.GRCh38.104.chr.gtf file into human_annotation_channel


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
    path "GenomeDir" into index_path

    script:
    """
    #!/usr/bin/env bash

    STAR --runThreadN 4 \
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
    file *.bam

    script:
    """
    #!/usr/bin/env bash

    STAR --outSAMstrandField intronMotif \
      --outFilterMismatchNmax 4 \
      --outFilterMultimapNmax 10 \
      --readFilesIn <(gunzip ${fastq[1][0]}) <(gunzip ${fastq[1][1]}) \
      --runThreadN 4 \
      --genomeDir GenomeDir \
      --outSAMunmapped None \
      --outSAMtype BAM SortedByCoordinate \
      --outStd BAM_SortedByCoordinate \
      --genomeLoad NoSharedMemory \
      --limitBAMsortRAM '10GB' \
      > ${fastq[0]}.bam
    """
}
mapping.view()
