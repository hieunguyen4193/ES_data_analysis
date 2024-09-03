new_dest_folder="/volume1/CRC1382_NAS01/CRC1382/OFFICIAL/raw_data/230605_Stange/cellranger_output_new";
mkdir -p ${new_dest_folder};

for sample in adult_GF d10_SPF d20_SPF d4_LPS d7_GF adult_SPF d15_SPF d4_GF d4_SPF d7_SPF SC11 SC12;do \
raw_output=/volume1/CRC1382_NAS01/CRC1382/OFFICIAL/raw_data/230605_Stange/cellranger_output/${sample}/outs/per_sample_outs/${sample}/count;

mkdir ${new_dest_folder}/${sample};
mkdir ${new_dest_folder}/${sample}/outs;
mkdir ${new_dest_folder}/${sample}/outs/filtered_feature_bc_matrix;
rsync -avh --progress ${raw_output}/sample_filtered_feature_bc_matrix/* ${new_dest_folder}/${sample}/outs/filtered_feature_bc_matrix;

rsync -avh --progress ${raw_output}/* ${new_dest_folder}/${sample}/outs --exclude=sample_filtered_feature_bc_matrix;

mv ${new_dest_folder}/${sample}/outs/sample_alignments.bam ${new_dest_folder}/${sample}/outs/possorted_genome_bam.bam;
mv ${new_dest_folder}/${sample}/outs/sample_alignments.bam.bai ${new_dest_folder}/${sample}/outs/possorted_genome_bam.bam.bai;done

