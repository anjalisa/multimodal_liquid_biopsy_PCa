### Script for Downsampling 

# Downsampling; script: https://github.com/sdchandra/tMAD
# since total reads are counted, the amount of paired reads have to be doupled

### !!!! samtools is used for this downsampling; samtools has to be in the path or rather has to be loaded as module: 
# module load samtools/1.9

###########################################

# Define variables
args <- commandArgs()

bam_path <- args[1]
total_reads_needed <- args[2]

###########################################

command <- paste0("ls ", bam_path, "/*.bam")
bamfiles <- as.character(system(command, intern=TRUE, wait=TRUE))

bamfiles <- unlist(lapply(bamfiles, function(x){
  y <- strsplit(x, "/")[[1]]
  y[length(y)]
  
}))

random_seed <- 8
total_reads_label <- paste0(round(total_reads_needed/1000000, 2), "Million")


for(t in 1:length(total_reads_needed)){
  
  total_reads <- total_reads_needed[t]
  
  for(i in 1:length(bamfiles)){
    
    bamfile <- paste0(bam_path, "/", bamfiles[i])
    
    if(file.exists(bamfile)){
      
      # Get the read count total from BAM
      command <- paste0("samtools view -c -q 20 -F 4 -F 2048 -F 256 -F 1024 ", bamfile)
      
      BAM_total <- as.numeric(system(command, intern=TRUE, wait=TRUE))
      
      if(total_reads >= BAM_total){
        
        
        stop(paste0("BAM has lower read count than ", 
                    total_reads_label[t], 
                    bamfiles[i]))         
        
      }else{
        
        sampling_prop <- 8 + round(total_reads/BAM_total, 3)
        
        command <- 
          paste0("samtools view -s ", 
                 sampling_prop, 
                 " -q 20 -F 4 -F 2048 -F 256 -F 1024 -b ", bamfile, " > ", 
                 "/samp_", total_reads_label[t], "_", bamfiles[i])
        
        system(command, intern=FALSE, wait=TRUE)
      }        
      
    }else{
      
      stop(paste0("File does not exists : ", bamfile))
      
    }
  }
}