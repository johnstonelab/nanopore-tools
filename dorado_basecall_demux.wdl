version 1.0

workflow dorado_basecall {

    input {
        String sample_id
        File fast5_archive
        String basecall_model
        String kitname
        Int disk_gb
    }

    call basecall {
        input:
            sample_id = sample_id,
            fast5_archive = fast5_archive,
            basecall_model = basecall_model,
            kitname = kitname,
            disk_gb = disk_gb
    }
    call demux {
        input:
            unmapped_bam = basecall.unmapped_bam
    }
    output {
        File unmapped_bam = basecall.unmapped_bam
        #File? duplex_unmapped_bam = basecall.duplex_unmapped_bam
        #File? pairs = basecall.pairs
        Array[File] bams = demux.bams
    }

    meta {
        author: "Martin Aryee"
        email:"martin.aryee@gmail.com"
    }
}
task basecall  {
    input {
        String sample_id
        File fast5_archive
        String basecall_model
        String kitname
        Int disk_gb #= ceil(size(fast5_archive, "GB")*3) + 5
    }
    command <<<
        filetype=$(file ~{fast5_archive})

        if [[ ~{fast5_archive} == *.pod5 ]]; then
            mkdir pod5s
            ln -s ~{fast5_archive} pod5s/reads.pod5
        fi

        if [[ "$filetype" == *"gzip compressed data"* ]]; then
          echo "FAST5s appear to be compressed with gzip. Decompressing..."
          mkdir fast5s
          tar zxvf ~{fast5_archive} -C fast5s
          # Convert FAST5 into POD5
          pod5 convert fast5 -r fast5s -o pod5s/reads.pod5 --threads 12 
        fi

        if [[ "$filetype" == *"Zip archive data"* ]]; then
          echo "FAST5s appear to be compressed with zip. Decompressing..."
          mkdir fast5s 
          unzip ~{fast5_archive} -d fast5s
          # Convert FAST5 into POD5
          pod5 convert fast5 -r fast5s -o pod5s/reads.pod5 --threads 12 
        fi
        
        # Simplex call with --emit-moves and --modified-bases
        dorado basecaller ~{kitname} --no-trim ~{basecall_model} pod5s --emit-moves > ~{sample_id}.unmapped.bam
        #/dorado_models/~{basecall_model}
        # Identify potential pairs
        # duplex_tools pair --output_dir ./pairs ~{sample_id}.unmapped.bam
    
        # Stereo duplex basecall
        # if [ -f "pairs/pair_ids_filtered.txt" ]; then
        #     NUM_PAIRS=$(wc -l pairs/pair_ids_filtered.txt)
        #     echo "${NUM_PAIRS} pairs found. Duplex calling..."
        #     dorado duplex /dorado_models/~{basecall_model} pod5s --pairs pairs/pair_ids_filtered.txt | samtools view -Sh > ~{sample_id}.duplex.unmapped.bam    
        # else 
        #     echo "No pairs found."
        # fi
    >>>
    runtime {
        gpuType: "nvidia-tesla-v100"
        gpuCount: 1
        cpu: 12
        disks: "local-disk " + disk_gb + " SSD" 
        memory: "32GB"  
        nvidiaDriverVersion: "470.161.03"
        zones: ["us-central1-a"] 
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/dorado"
    }
    output {
        File unmapped_bam = "~{sample_id}.unmapped.bam"
        # File? duplex_unmapped_bam = "~{sample_id}.duplex.unmapped.bam"
        # File? pairs = "pairs/pair_ids_filtered.txt"
    }
}
task demux {
    input {
        File unmapped_bam
    }
    command <<<
        mkdir out
        dorado demux -o ./out --no-classify ~{unmapped_bam} 
    >>>
    runtime {
        gpuType: "nvidia-tesla-v100"
        gpuCount: 1
        cpu: 12
        disks: "local-disk 500 SSD"
        memory: "32GB"  
        nvidiaDriverVersion: "470.161.03"
        zones: ["us-central1-a"] 
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/dorado"
    }
    output {
        Array[File] bams = glob("./out/*.bam")
    }
}