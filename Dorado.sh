 #!/bin/bash
dorado basecaller $MODELS_DIR/dna_r10.4.1_e8.2_400bps_sup@v4.3.0 pod5/ --min-qscore 10 --kit-name SQK-RBK114-24 -x "cuda:$CUDA_VISIBLE_DEVICES" > calls.bam

dorado demux --output-dir demux/ --no-classify --emit-fastq calls.bam
