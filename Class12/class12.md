---
title: "class12"
author: "Helena Tubb"
date: "5/10/2018"
output: 
  html_document: 
    keep_md: yes
---




```r
library(bio3d)
file.name <- get.pdb("1hsg")
```

```
## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download
```


```r
hiv <- read.pdb(file.name)
```


```r
prot <- trim.pdb(hiv, "protein")
lig <- trim.pdb(hiv, "ligand")

write.pdb(prot, file="1hsg_protein.pdb")
write.pdb(lig, "1hsg_ligand.pdb")
```




