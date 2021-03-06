args = commandArgs(trailingOnly=TRUE)
library("seqinr")
dna.files <- dir(pattern =".dna.fasta")
protein.files <- dir(pattern =".protein.fasta")
distance <- c()
for(i in 1:length(dna.files))
{
  outbase="aln.fa"
  dna.out=paste (dna.files[i],outbase, sep=".")
  reverse.align(nucl.file = dna.files[i], out.file = dna.out, align.prot=TRUE, clustal.path= args[1])
}

