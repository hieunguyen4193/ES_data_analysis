export PATH=/home/hieunguyen/samtools/bin:$PATH
path_to_modified_gtf="/media/hieunguyen/HD01/storage/build-mm10/mm10-2020-A_build/gencode.vM23.primary_assembly.annotation.gtf.filtered";

samtools_threads=7;

path_to_cellranger_output="/media/hieunguyen/HD01/cellranger_output_new"

path_to_masked_regions_gtf="/media/hieunguyen/HD01/storage/mm10_rmsk.gtf";

for sample in SC11 SC12;do \
velocyto run10x -m ${path_to_masked_regions_gtf} ${path_to_cellranger_output}/${sample} ${path_to_modified_gtf} -@ ${samtools_threads};done
