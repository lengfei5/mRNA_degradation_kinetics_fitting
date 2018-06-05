###
### real data
###
do.step = TRUE;
Real.Data = FALSE
TEST = FALSE

if(Real.Data)
{
  #data.version = '_total_v2'
  #res.version = '_total_all_v2';
  data.version = '_total_counts_v2'
  res.version = '_total_counts_all_genes_norm_params_v6';
  
}else{
  #fake.version = '_total_counts_s1'
  #fake.version = '_total_counts_s4_nonident_model3';
  #res.version = '_total_counts_s4_nonident_model3_v1';
  fake.version = '_total_counts_s5_all'
  res.version = '_total_counts_s5_all_v1';
}

if(do.step)
{
  #source('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/model_RNA_seq_total/functions.R')
  if(Real.Data)
  {
    cat('FIT THE REAL DATA ON VITAL-IT\n')
    
    load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_analysis_sel_alphas', data.version, '.Rdata', sep=''))
    #source("/Users/jiwang/Degradation_Liver/Main_Code_Degradation/f24_modified_1.0.r")
    T = R;
    gene_list = T$gene;
    #load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_analysis_sel', data.version, '.Rdata', sep=''))
    #T = Tt;
    
    ### Prepare the data for the fits
    if(TEST)
    {
      #examples = c('Arntl','Clock', 'Per1', 'Per2', 'Per3', 'Cry1', 'Cry2','Npas2', 'Nr1d1', 'Nr1d2', 'Rora', 'Rorc', 'Dbp', 'Hlf', 'Tef', 'Nfil3', 
      #'Bhlhe40', 'Bhlhe41', 'Klf10', 'Klf3', 'Klf13', 'Hsf1', 'Crebbp',
      #             'Abcb11', 'Tfrc', 'Nedd4l', 'Loxl4', 'Tubb2a', 'Cbs', 'Rcan1', 'Mfsd2', 'Gsta3', 'Cdkn1a', 'Ppara', 'Ppard', 'Nr3c1', 'Nr2c6')
      examples = c('Per2', 'Cry1', 'Per3','Per1',  'Cry2')   
                   
      #mm = match(unique(examples), T[,1])
      #mm = mm[order(-T[mm, 5])]
      #examples = examples[which(!is.na(mm)==TRUE)]
      #examples = c(examples, as.character(T$gene[order(T$qv.rpkm.mRNA)][1:3000]))
      examples = unique(examples)
      
      mm = match(examples, gene_list)
      mm = mm[which(!is.na(mm)==TRUE)]
      #mm = mm[c(1:20)]
      gene_list = gene_list[mm]
      T = T[mm, ]
    }
  }else{
    cat('FIT THE SIMULATED DATA ON VITAL-IT\n')
    #load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_simulated_data', fake.version, '.Rdata', sep=''));
    #load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_simulated_data', fake.version, '.Rdata', sep=''))
    load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_simulated_data', fake.version, '.Rdata', sep='')) 
    #T = F0;
    T = F;
    gene_list = T$gene
    #gene_list = T$gene[c(1:10)];
    #gene_list = 'fake_m4_157';
    #T = T[match(gene_list, T$gene), ]
    #gene_list = T$gene[c(1:10)];
  }
  
  save(T, file = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/stuff_for_vital_IT/genes_ready_for_fits.Rdata')
  write.table(gene_list, file = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/stuff_for_vital_IT/gene_list.txt', quote = FALSE, sep = '\n', row.names = FALSE, col.names = FALSE)
  
  ## ON VITAL-IT 
	
	#clear the directory
	cmd = "ssh jwang@dev.vital-it.ch 'rm -r /scratch/cluster/monthly/jingkui/microarray/*' " 
	system(cmd)
	
	#copy the data and the code on vital-IT
	cmd = 'scp -r /Users/jiwang/Degradation_Liver/Main_Code_Degradation/stuff_for_vital_IT/* jwang@dev.vital-it.ch:/scratch/cluster/monthly/jingkui/microarray/.'
	system(cmd)
	cmd = 'scp -r /Users/jiwang/Degradation_Liver/Main_Code_Degradation/model_RNA_seq_total/functions.R jwang@dev.vital-it.ch:/scratch/cluster/monthly/jingkui/microarray/.'
	system(cmd)
  cmd = 'scp -r /Users/jiwang/Degradation_Liver/Main_Code_Degradation/model_RNA_seq_total/Scaling_factors_48_samples.Rdata jwang@dev.vital-it.ch:/scratch/cluster/monthly/jingkui/microarray/.'
  system(cmd)
  #cmd = 'scp -r /Users/jiwang/Degradation_Liver/Main_Code_Degradation/model_RNA_seq_total/Libary_size_48_samples.Rdata jwang@dev.vital-it.ch:/scratch/cluster/monthly/jingkui/microarray/.'
  #system(cmd)
  
	#########################
  # run the fits : bash -s <to_run_on_vital_IT.sh
	cmd = "ssh jwang@dev.vital-it.ch 'bash -s' < /Users/jiwang/Degradation_Liver/Main_Code_Degradation/stuff_for_vital_IT/to_run_on_vital_IT.sh"
	system(cmd)
	##########################
	
	# concatenate all fits results together when fits are all done
	cmd = "ssh jwang@dev.vital-it.ch 'perl /scratch/cluster/monthly/jingkui/microarray/cat_all_fits_results.pl'"
	system(cmd)
	
	# wait for everything to be done and copy the data on the local computer
	cmd = "ssh jwang@dev.vital-it.ch 'ls -l /scratch/cluster/monthly/jingkui/microarray/fits/all_genes_fits_results_with_header.txt'"
	res = system(cmd, intern = TRUE)
	
	while(length(res) == 0)
	{
		Sys.sleep(60*5) # wait for 5 minutes
		cat('... wait for the fits to finish\n')
		cmd = "ssh jwang@dev.vital-it.ch 'ls -l /scratch/cluster/monthly/jingkui/microarray/fits/all_genes_fits_results_with_header.txt'"
		res = system(cmd, intern = TRUE)
	}
	
	cat('........................................ fits finished on vital-IT\n')	
	
	cmd = paste("scp jwang@dev.vital-it.ch:/scratch/cluster/monthly/jingkui/microarray/fits/all_genes_fits_results_with_header.txt /Users/jiwang/Degradation_Liver/Main_Code_Degradation/Vital-IT_fits_optim/all_genes_fits_results_with_header_RNA_seq_beta", res.version, ".txt", sep = '')
	system(cmd)
	
	cat('........................................ fits results imported on local computer\n')
	
	file =  paste("/Users/jiwang/Degradation_Liver/Main_Code_Degradation/Vital-IT_fits_optim/all_genes_fits_results_with_header_RNA_seq_beta", res.version, ".txt", sep = '')
	fits.results = read.table(file = file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
	#fits.results = fits.results[-8528, ]
    #colnames(fits.results) = c('gene', names)
    
	if(nrow(fits.results)<0.9*length(gene_list)){cat('........................................ fits results missing for more than 10% of genes\n')}
	m = match(gene_list, fits.results$gene)
	if(sum(is.na(m))>0){cat('........................................ fits results missing for',sum(is.na(m)),'genes\n')}
	fits.results = fits.results[m,]
	
	#kk = match(gene_list, T$gene)
	#Tt = Tt[kk,]
	T = data.frame(T, fits.results[,-1],stringsAsFactors = FALSE)	
	save(T, file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
	
}
