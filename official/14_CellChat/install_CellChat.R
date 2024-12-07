##### install cellchat v2 from my github repo

if ("CellChat" %in% installed.packages()){
  remove.packages("CellChat")
}

devtools::install_github("immunogenomics/presto", upgrade = "never")
devtools::install_github("https://github.com/hieunguyen4193/CellChatv2_862ab34", upgrade = "never")

##### temporary fix --> https://github.com/jinworks/CellChat/issues/202
# cellchat.sample1@data.smooth <- cellchat.sample1@data.project
# cellchat.sample2@data.smooth <- cellchat.sample2@data.project
##### issue: non conformable array: https://github.com/sqjin/CellChat/issues/708

reticulate::install_python(version = '3.8')
reticulate::py_install(packages = 'umap-learn')