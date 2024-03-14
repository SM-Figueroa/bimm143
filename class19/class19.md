# class19: Cancer mutation mini-project
Sawyer (PID: A16335277)

``` r
# load in fasta file in current dir
library(bio3d)

seqs <- read.fasta("A16335277_mutant_seq.fa")
```

``` r
# find differences in two sequences in fasta file

inds = which(conserv(seqs) < 1)
print(inds)
```

    [1] 736 841 864 894

``` r
seqs$ali[,inds]
```

                 [,1] [,2] [,3] [,4]
    wt_healthy   "G"  "P"  "Y"  "G" 
    mutant_tumor "V"  "E"  "R"  "Y" 

``` r
# select 724-1005 amino acids in mutant sequence
seqs$ali[,724:1005]
```

                 [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
    wt_healthy   "L"  "V"  "L"  "G"  "K"  "T"  "L"  "G"  "E"  "G"   "E"   "F"  
    mutant_tumor "L"  "V"  "L"  "G"  "K"  "T"  "L"  "G"  "E"  "G"   "E"   "F"  
                 [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23]
    wt_healthy   "G"   "K"   "V"   "V"   "K"   "A"   "T"   "A"   "F"   "H"   "L"  
    mutant_tumor "V"   "K"   "V"   "V"   "K"   "A"   "T"   "A"   "F"   "H"   "L"  
                 [,24] [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34]
    wt_healthy   "K"   "G"   "R"   "A"   "G"   "Y"   "T"   "T"   "V"   "A"   "V"  
    mutant_tumor "K"   "G"   "R"   "A"   "G"   "Y"   "T"   "T"   "V"   "A"   "V"  
                 [,35] [,36] [,37] [,38] [,39] [,40] [,41] [,42] [,43] [,44] [,45]
    wt_healthy   "K"   "M"   "L"   "K"   "E"   "N"   "A"   "S"   "P"   "S"   "E"  
    mutant_tumor "K"   "M"   "L"   "K"   "E"   "N"   "A"   "S"   "P"   "S"   "E"  
                 [,46] [,47] [,48] [,49] [,50] [,51] [,52] [,53] [,54] [,55] [,56]
    wt_healthy   "L"   "R"   "D"   "L"   "L"   "S"   "E"   "F"   "N"   "V"   "L"  
    mutant_tumor "L"   "R"   "D"   "L"   "L"   "S"   "E"   "F"   "N"   "V"   "L"  
                 [,57] [,58] [,59] [,60] [,61] [,62] [,63] [,64] [,65] [,66] [,67]
    wt_healthy   "K"   "Q"   "V"   "N"   "H"   "P"   "H"   "V"   "I"   "K"   "L"  
    mutant_tumor "K"   "Q"   "V"   "N"   "H"   "P"   "H"   "V"   "I"   "K"   "L"  
                 [,68] [,69] [,70] [,71] [,72] [,73] [,74] [,75] [,76] [,77] [,78]
    wt_healthy   "Y"   "G"   "A"   "C"   "S"   "Q"   "D"   "G"   "P"   "L"   "L"  
    mutant_tumor "Y"   "G"   "A"   "C"   "S"   "Q"   "D"   "G"   "P"   "L"   "L"  
                 [,79] [,80] [,81] [,82] [,83] [,84] [,85] [,86] [,87] [,88] [,89]
    wt_healthy   "L"   "I"   "V"   "E"   "Y"   "A"   "K"   "Y"   "G"   "S"   "L"  
    mutant_tumor "L"   "I"   "V"   "E"   "Y"   "A"   "K"   "Y"   "G"   "S"   "L"  
                 [,90] [,91] [,92] [,93] [,94] [,95] [,96] [,97] [,98] [,99] [,100]
    wt_healthy   "R"   "G"   "F"   "L"   "R"   "E"   "S"   "R"   "K"   "V"   "G"   
    mutant_tumor "R"   "G"   "F"   "L"   "R"   "E"   "S"   "R"   "K"   "V"   "G"   
                 [,101] [,102] [,103] [,104] [,105] [,106] [,107] [,108] [,109]
    wt_healthy   "P"    "G"    "Y"    "L"    "G"    "S"    "G"    "G"    "S"   
    mutant_tumor "P"    "G"    "Y"    "L"    "G"    "S"    "G"    "G"    "S"   
                 [,110] [,111] [,112] [,113] [,114] [,115] [,116] [,117] [,118]
    wt_healthy   "R"    "N"    "S"    "S"    "S"    "L"    "D"    "H"    "P"   
    mutant_tumor "R"    "N"    "S"    "S"    "S"    "L"    "D"    "H"    "E"   
                 [,119] [,120] [,121] [,122] [,123] [,124] [,125] [,126] [,127]
    wt_healthy   "D"    "E"    "R"    "A"    "L"    "T"    "M"    "G"    "D"   
    mutant_tumor "D"    "E"    "R"    "A"    "L"    "T"    "M"    "G"    "D"   
                 [,128] [,129] [,130] [,131] [,132] [,133] [,134] [,135] [,136]
    wt_healthy   "L"    "I"    "S"    "F"    "A"    "W"    "Q"    "I"    "S"   
    mutant_tumor "L"    "I"    "S"    "F"    "A"    "W"    "Q"    "I"    "S"   
                 [,137] [,138] [,139] [,140] [,141] [,142] [,143] [,144] [,145]
    wt_healthy   "Q"    "G"    "M"    "Q"    "Y"    "L"    "A"    "E"    "M"   
    mutant_tumor "Q"    "G"    "M"    "Q"    "R"    "L"    "A"    "E"    "M"   
                 [,146] [,147] [,148] [,149] [,150] [,151] [,152] [,153] [,154]
    wt_healthy   "K"    "L"    "V"    "H"    "R"    "D"    "L"    "A"    "A"   
    mutant_tumor "K"    "L"    "V"    "H"    "R"    "D"    "L"    "A"    "A"   
                 [,155] [,156] [,157] [,158] [,159] [,160] [,161] [,162] [,163]
    wt_healthy   "R"    "N"    "I"    "L"    "V"    "A"    "E"    "G"    "R"   
    mutant_tumor "R"    "N"    "I"    "L"    "V"    "A"    "E"    "G"    "R"   
                 [,164] [,165] [,166] [,167] [,168] [,169] [,170] [,171] [,172]
    wt_healthy   "K"    "M"    "K"    "I"    "S"    "D"    "F"    "G"    "L"   
    mutant_tumor "K"    "M"    "K"    "I"    "S"    "D"    "F"    "Y"    "L"   
                 [,173] [,174] [,175] [,176] [,177] [,178] [,179] [,180] [,181]
    wt_healthy   "S"    "R"    "D"    "V"    "Y"    "E"    "E"    "D"    "S"   
    mutant_tumor "S"    "R"    "D"    "V"    "Y"    "E"    "E"    "D"    "S"   
                 [,182] [,183] [,184] [,185] [,186] [,187] [,188] [,189] [,190]
    wt_healthy   "Y"    "V"    "K"    "R"    "S"    "Q"    "G"    "R"    "I"   
    mutant_tumor "Y"    "V"    "K"    "R"    "S"    "Q"    "G"    "R"    "I"   
                 [,191] [,192] [,193] [,194] [,195] [,196] [,197] [,198] [,199]
    wt_healthy   "P"    "V"    "K"    "W"    "M"    "A"    "I"    "E"    "S"   
    mutant_tumor "P"    "V"    "K"    "W"    "M"    "A"    "I"    "E"    "S"   
                 [,200] [,201] [,202] [,203] [,204] [,205] [,206] [,207] [,208]
    wt_healthy   "L"    "F"    "D"    "H"    "I"    "Y"    "T"    "T"    "Q"   
    mutant_tumor "L"    "F"    "D"    "H"    "I"    "Y"    "T"    "T"    "Q"   
                 [,209] [,210] [,211] [,212] [,213] [,214] [,215] [,216] [,217]
    wt_healthy   "S"    "D"    "V"    "W"    "S"    "F"    "G"    "V"    "L"   
    mutant_tumor "S"    "D"    "V"    "W"    "S"    "F"    "G"    "V"    "L"   
                 [,218] [,219] [,220] [,221] [,222] [,223] [,224] [,225] [,226]
    wt_healthy   "L"    "W"    "E"    "I"    "V"    "T"    "L"    "G"    "G"   
    mutant_tumor "L"    "W"    "E"    "I"    "V"    "T"    "L"    "G"    "G"   
                 [,227] [,228] [,229] [,230] [,231] [,232] [,233] [,234] [,235]
    wt_healthy   "N"    "P"    "Y"    "P"    "G"    "I"    "P"    "P"    "E"   
    mutant_tumor "N"    "P"    "Y"    "P"    "G"    "I"    "P"    "P"    "E"   
                 [,236] [,237] [,238] [,239] [,240] [,241] [,242] [,243] [,244]
    wt_healthy   "R"    "L"    "F"    "N"    "L"    "L"    "K"    "T"    "G"   
    mutant_tumor "R"    "L"    "F"    "N"    "L"    "L"    "K"    "T"    "G"   
                 [,245] [,246] [,247] [,248] [,249] [,250] [,251] [,252] [,253]
    wt_healthy   "H"    "R"    "M"    "E"    "R"    "P"    "D"    "N"    "C"   
    mutant_tumor "H"    "R"    "M"    "E"    "R"    "P"    "D"    "N"    "C"   
                 [,254] [,255] [,256] [,257] [,258] [,259] [,260] [,261] [,262]
    wt_healthy   "S"    "E"    "E"    "M"    "Y"    "R"    "L"    "M"    "L"   
    mutant_tumor "S"    "E"    "E"    "M"    "Y"    "R"    "L"    "M"    "L"   
                 [,263] [,264] [,265] [,266] [,267] [,268] [,269] [,270] [,271]
    wt_healthy   "Q"    "C"    "W"    "K"    "Q"    "E"    "P"    "D"    "K"   
    mutant_tumor "Q"    "C"    "W"    "K"    "Q"    "E"    "P"    "D"    "K"   
                 [,272] [,273] [,274] [,275] [,276] [,277] [,278] [,279] [,280]
    wt_healthy   "R"    "P"    "V"    "F"    "A"    "D"    "I"    "S"    "K"   
    mutant_tumor "R"    "P"    "V"    "F"    "A"    "D"    "I"    "S"    "K"   
                 [,281] [,282]
    wt_healthy   "D"    "L"   
    mutant_tumor "D"    "L"   
