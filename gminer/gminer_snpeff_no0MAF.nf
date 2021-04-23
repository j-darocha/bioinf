#!/usr/bin/env nextflow

vcf_in_ch = Channel.fromPath("${params.VCFDIR}/*.vcf.gz")
outputdir = params.OUTFOLDER
dirplugins = params.dirplugins
dircache = params.cache
gene_list = params.gene_list

process check_nb_chromosomes {
        input:
                file vcf from vcf_in_ch
        output:
                file "*.vcf" into decompressed_vcf
                file "$vcf" into input_to_down_prcesses

        """
        output_filename=\$(basename $vcf .gz)
        zcat $vcf >\$output_filename 
        """
        }

// Check if vcf files contains variants from one chromosome 
// Workflow will stop at exit state 1 
decompressed_vcf.buffer{ myFile = file(it)
                  allLines  = myFile.readLines()
                  def myList = []
                  for( line in  allLines ) {
                        if (  !(line  =~ /\#+.+/ )  )
                            { mychromosome =  line.split()
                              chrom =  mychromosome[0].replaceAll("[^a-zA-Z0-9 ]+","")
                                if ( !(chrom in myList)  )
                                        { myList = myList + chrom  }
                                                }
                                                        }
                if ( myList.size() == 1  ) { println " file $it.name contains one chromosome  " }
                else { println "WARNING: file ${it.name}.gz contains more than one chrmosome, workflow is stopping"
                        System.exit(1)  }
                 myFile.delete()
                         }

input_to_down_prcesses.into { vcf_to_check; vcf_for_process ; vcf_for_freq ; vcfforstat;  vcf_TOtrim; origvcf_forVEP }

process split_samples {   
    input:
         file(vcf) from vcf_for_process
    output: 
         set val(name), file("${name}_split.vcf") into splitted_samples_vcf
        script:
                name = vcf.baseName.replaceFirst(".vcf","")
   """
   bcftools norm -m - $vcf > ${name}_split.vcf
   """
}

// will add the ID to the empty column AND remove non zero
process add_IDS {
        input:
                set val(name), file(vcf) from splitted_samples_vcf
        output:
                set val(name), file("${name}_IDs.vcf") into (vcf_for_stat, vcf_for_allel_freq, VEP_ready_vcfs, SnpEff_ready_vcfs, vcf_for_phasing )

        """
        bcftools annotate --set-id "%CHROM\\_%POS\\_%REF\\_%FIRST_ALT"  $vcf -O z -o output.vcf.gz
	bcftools view output.vcf.gz -i 'MAF[0]!=0' > ${name}_IDs.vcf
	
        """
}

// generates stats file
process getStats {
        input:
                set val(name), file(vcf) from vcf_for_stat
        output:
                set val(name), file("${name}_stats.txt") into stat_output

        """
        bcftools stats $vcf  > ${name}_stats.txt
        """
}

// generates stats file
process parseStats {
        input: 
                set val(name), file(stats) from  stat_output
                
        output: 
                set val(name), file("${name}_statFile") into statFile
        publishDir "${outputdir}/${name}", mode:'copy'

        script:
        """
        #!/usr/bin/env python
        from re import match 
        with open( "$stats" ) as file:
                myfile = file.readlines()
        stats = []
        for index, line in enumerate(myfile) : 
                is_SN = bool( match(r'^SN', line) )
                if is_SN : 
                        stats.append(line.split('\\t')[-1].replace('\\n', '') ) 
                if '# QUAL,' in line: 
                        stats.append( myfile[index+2].split()[2] ) 
        col_names = ['nb_samples', 'mb_records', 'no_ALTs', 'nb_SNPs', 'mb_MNPs', 'nb_indels', 'nb_others', 'nb_multAllelic', 'nb_multAllelicSNPs' , 'lowest_qaulity_score']
        
        with open( "${name}_statFile" , 'w') as output:
                output.writelines( ','.join( col_names  )+'\\n' )
                output.writelines( ','.join( stats  )+'\\n' )            
        """
}

// calculate allele frequencies 
process allelFreq {
        input: 
                set val(name), file(vcf) from vcf_for_allel_freq 
        output: 
                set val(name), file("${name}_freqFilename.frq") into allelFreqFile

        publishDir "${outputdir}/${name}", mode:'copy'
                
        """
        plink --vcf $vcf --freq --keep-allele-order --out ${name}_freqFilename
        """
}


// annotation with SnpEff 
process snpeff_annotation {
	maxForks 50
        input:
            set val(name), file(vcf) from SnpEff_ready_vcfs
        output:
              file("${name}_sneff.vcf.gz") into vcf_snpeff
              file("*.txt") into summary_text
              file("*.html") into html_report
              file("*.csv") into csv_snpeff

                    
        publishDir "${outputdir}/${name}", mode:'copy'
        """
        java -Xmx4g  -jar "${params.snpeff}"/snpEff.jar  hg19 -canon -csvStats mystats.csv -stats ${name}  $vcf   > ${name}_sneff.vcf.gz 
        mv ${name} ${name}.html
	mv mystats.csv ${name}.csv  
        """
}

//run vep
process vep {
        module 'perl526'
  	maxForks 8
  	input:
        	set val(name), file(vcf) from VEP_ready_vcfs            
  	output:
   		set val(name), file("${name}_vepFilename.vep") into vep_ch  

  	publishDir "${outputdir}/${name}", mode: 'copy'
  	script:
                
    """ 
    vep --ASSEMBLY GRCh37  -i $vcf -o ${name}_vepFilename.vep --dir_cache $dircache --cache --port 3337 --symbol --canonical --ccds --tab --check_existing --sift s --dir_plugins $dirplugins --plugin CADD,/dataB/aux/37/cadd/whole_genome_SNVs.tsv.gz --sift b --polyphen b --plugin Condel,/dataB/aux/37/cadd/config/Condel/,b,2 --plugin dbNSFP,/dataB/aux/37/cadd/dbNSFP_hg19.gz,MetaSVM_pred,MetaSVM_score,MetaSVM_rankscore,MetaLR_pred,MetaLR_score,MetaLR_rankscore,VEST3_score,VEST3_rankscore,FATHMM_pred,FATHMM_score,FATHMM_converted_rankscore,CADD_phred,CADD_raw,CADD_raw_rankscore,DANN_rankscore,DANN_score,LRT_converted_rankscore,LRT_pred,LRT_score,fathmm-MKL_coding_group,fathmm-MKL_coding_pred,fathmm-MKL_coding_rankscore,fathmm-MKL_coding_score,PROVEAN_converted_rankscore,PROVEAN_pred,PROVEAN_score,MutationAssessor_UniprotID,MutationAssessor_pred,MutationAssessor_score,MutationAssessor_score_rankscore,MutationAssessor_variant --plugin LoF,loftee_path:/opt/exp_soft/bioinf/ensembl-vep/plugins/loftee,human_ancestor_fa:/opt/exp_soft/bioinf/ensembl-vep/plugins/loftee/human_ancestor.fa.gz

    """
}



joined_vepfrqs = vep_ch.join(allelFreqFile)

//vep_freq merging - you need the python script from bin
process merging {
	input:
		set val(name), file(vep), file(frq) from joined_vepfrqs
	output:	
		set val(name), file("${name}_vepfreqmerged.txt") into merged_ch 
		publishDir "${outputdir}/${name}", mode: 'copy'
	script:
	"""
	vep-frq_merge_args.py --vepfile $vep --frqfile $frq --output ${name}
	"""
}


//variant counting and filtering - vepbased - you need the python script from bin
process countANDfilter {
        input:
                set val(name), file(mvepfrq) from merged_ch
        output:
		set val(name), file("${name}_canonical_adme_only_1p.txt") into counted_ch
		set val(name), file("${name}_canonical_adme_only_5p.txt")
		set val(name), file("${name}_canonical_adme_only.txt")
		set val(name), file("${name}_canonical_counts_adme_only_1p.txt")
		set val(name), file("${name}_canonical_counts_adme_only_5p.txt")
		set val(name), file("${name}_canonical_counts_adme_only.txt")
		set val(name), file("${name}_missense_can_1of4_only_1p.txt")
		set val(name), file("${name}_missense_can_1of4_only_5p.txt")
		set val(name), file("${name}_missense_can_1of4_only.txt")
		set val(name), file("${name}_missense_can_adme_only_1p.txt")
		set val(name), file("${name}_missense_can_adme_only_5p.txt")
		set val(name), file("${name}_missense_can_adme_only.txt")
		set val(name), file("${name}_missense_can_all4tools_only_1p.txt")
		set val(name), file("${name}_missense_can_all4tools_only_5p.txt")
		set val(name), file("${name}_missense_can_all4tools_only.txt")
		set val(name), file("${name}_ADME_model_Pred_1o5.txt")
                set val(name), file("${name}_ADME_model_Pred_1o5_1p.txt")
                set val(name), file("${name}_ADME_model_Pred_1o5_5p.txt")
                set val(name), file("${name}_ADME_model_Pred_ALL.txt")
                set val(name), file("${name}_ADME_model_Pred_ALL_1p.txt")
                set val(name), file("${name}_ADME_model_Pred_ALL_5p.txt")
                set val(name), file("${name}_HC_LOF_noncan.txt")
                set val(name), file("${name}_LC_LOF_noncan.txt")
                set val(name), file("${name}_HC_LOF_can.txt")
                set val(name), file("${name}_LC_LOF_can.txt")
                set val(name), file("${name}_HC_LOF_noncan_1p.txt")
                set val(name), file("${name}_HC_LOF_can_1p.txt")




                publishDir "${outputdir}/${name}", mode: 'copy'
                
        script:
        """
        Counting_Filtering_vepfrqs_args.py --vepfile $mvepfrq --list ${gene_list} --output ${name}
        """
}
