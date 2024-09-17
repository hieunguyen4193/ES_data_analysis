#####----------------------------------------------------------------------------#####
##### install cluster profiler
#####----------------------------------------------------------------------------#####
# remove.packages("DOSE")
# remove.packages("GOSemSim")
# remove.packages("yulab.utils")
# remove.packages("clusterProfiler")
devtools::install_github("YuLab-SMU/yulab.utils", upgrade = "never")
remotes::install_github("GuangchuangYu/GOSemSim", upgrade = "never")
remotes::install_github("GuangchuangYu/clusterProfiler", upgrade = "never")