# UniDesign

## Introduction
UniDesign is a computational framework for protein design, targeting a diversity of protein design and engineering tasks, and it can also be used for protein structure modeling and scoring. UniDesign is extended mainly from the EvoEF2 protein design program with the evolutionary feature taken from EvoDesign, but it has many new features. More details can be found in <a href=https://github.com/tommyhuangthu/UniDesign/blob/master/manual.docx>manual.docx</a> within this package.

## Main Applications

#### Protein Design

• Design monomer protein <br>
• Design protein-protein interaction <br>
• Design protein-ligand interaction <br>
• Design protein-nucleic acid interaction <br>
• Design enzymes <br>
<br>
All protein design tasks are conducted with UniDesign command <i>ProteinDesign</i> with specific options. See <a href=https://github.com/tommyhuangthu/UniDesign/blob/master/manual.docx>manual.docx</a> for details and examples.

#### Protein Structure Modeling

• Protein side-chain packing (command <i>ProteinDesign</i> with option --wildtype) <br>
• Repair incomplete protein sidechains (command <i>RepairStructure</i>) <br>
• Protein minimization to remove sidechain clashes (command <i>Minimize</i> <br>
• Build mutant structural models (command <i>BuildMutant</i>) <br>
• Add polar hydrogen atoms (command <i>AddPolarHydrogen</i>) <br>
• Optimize hydrogen atom's position (command <i>OptimizeHydrogen</i>) <br>

#### Protein Scoring

• Compute protein stability (command <i>ComputeStability</i>) <br>
• Compute protein interchain binding interaction (command <i>ComputeBinding</i>) <br>


## Installation, Usage and Tutorial
Please refer to the <a href=https://github.com/tommyhuangthu/UniDesign/blob/master/manual.docx>manual.docx</a> for details.

## Copyright
Copyright (c) Xiaoqiang Huang. UniDesign is free to academic users. For suggestions, please contact xiaoqiah@umich.edu or xiaoqiah@outlook.com.

## References
Huang, X., Zhou, J., Yang, D., Zhang, J., Xia, X., Chen, Y. E., and Xu, J. Decoding CRISPR–Cas PAM recognition with UniDesign. Briefings in Bioinformatics. 2023, 24(3):bbad133. doi: 10.1093/bib/bbad133.
