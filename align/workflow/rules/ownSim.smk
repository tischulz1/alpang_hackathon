#Pipeline to simulate reads using own script

rule simulateReads:
    input:
        pjoin(OWN_SCRIPT_DIR, "extracted_genome_sequences", "{haplotype}.fa")
    output:
        pjoin(OWN_SCRIPT_DIR, "simulated_read_sequences", "{haplotype}.testReads.fa")
    shell:
        """
        python3 scripts/simReads.py -dp 10 -lmn 250 -lmx 250 -lavg 250 -ls 1 -sr 0.01 -dr 0.005 -ir 0.005 -r {input} -o {output}
        """
