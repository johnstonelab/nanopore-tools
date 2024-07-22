version 1.0

workflow bamtobigwig {

    input {
        File genome
        File chromsizes
        File fastq
        String modmotif
        String modkitoptions
        String mapqmin
    }

    call minimapalign {
        input:
            bam = bam,
            genome = genome
    }
    call filter {
        input:
            sortedbam = minimapalign.sortedbam,
            mapqmin = mapqmin
    }
    call tobedgraph {
        input:
            genome = genome,
            filteredbai = filter.filteredbai,
            filteredbam = filter.filteredbam,
            modmotif = modmotif,
            modkitoptions = modkitoptions
    }
    call tobigwig {
        input:
            filteredbedgraph = tobedgraph.filteredbedgraph,
            chromsizes = chromsizes
    }
    output {
        File filteredcovbedgraph = tobedgraph.filteredcovbedgraph
        File filteredbw = tobigwig.filteredbw
    }

    meta {
        author: "Martin Aryee"
        email:"martin.aryee@gmail.com"
    }
}
task minimapalign {
    input {
        File bam
        File genome
    }
    command <<<
    samtools fastq -T "*" ~{bam} | minimap2 -ax map-ont -y ~{genome} - | samtools sort -T tmp -o sorted.bam
    samtools index sorted.bam
    >>>
    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/minimap2:latest"
		memory: "64G"
		disks: "local-disk 500 SSD"
		cpu: 8
    }
    output {
        File sortedbam = "sorted.bam"
    }
}
task filter {
    input {
        File sortedbam
        String mapqmin
    }
    command <<<
    samtools view -bh -q ~{mapqmin} ~{sortedbam} > filtered.bam
    samtools index filtered.bam
    >>>
    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/samtools:latest"
		memory: "64G"
		disks: "local-disk 500 SSD"
		cpu: 8
    }
    output {
        File filteredbam = "filtered.bam"
        File filteredbai = "filtered.bam.bai"
    }
}
task tobedgraph {
    input {
        File genome
        File filteredbam
        File filteredbai
        String modmotif
        #specify motif with modmotif - i.e. CG 0 for cpg methylation or GC 1 for gpc
        String modkitoptions
        #any other options that are needed
    }
    command <<<
    modkit pileup ~{filteredbam} filtered.bed --motif ~{modmotif} ~{modkitoptions} --ref ~{genome} --ignore h -t 12 --combine-strands
    awk '$10 > 0 {printf "%s\t%d\t%d\t%d\t%2.3f\n" , $1,$2,$3,$10,$11}' filtered.bed | sort -k1,1 -k2,2n > filtered_coverage.bedgraph
    awk '$10 > 0 {printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$11}' filtered.bed | sort -k1,1 -k2,2n > filtered.bedgraph
    >>>
    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/modkit:latest"
		memory: "64G"
		disks: "local-disk 500 SSD"
		cpu: 8
    }
    output {
        File filteredcovbedgraph = "filtered_coverage.bedgraph"
        File filteredbedgraph = "filtered.bedgraph"
    }
}
task tobigwig {
    input {
        File filteredbedgraph
        File chromsizes
    }
    command <<<
    bedGraphToBigWig ~{filteredbedgraph} ~{chromsizes} filtered.bw
    >>>
    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/bedtools:latest"
		memory: "64G"
		disks: "local-disk 500 SSD"
		cpu: 8
    }
    output {
        File filteredbw = "filtered.bw"
    }
}