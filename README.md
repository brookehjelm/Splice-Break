# Splice-Break  v1.0.0<br/>
Filters Mapsplice2 alignment output to prioritize mtDNA deletion breakpoints for further inference. <br/>

# Dependencies<br/>        
[Mapsplice](http://www.netlab.uky.edu/p/bioinfo/MapSplice2) <br/>
[samtools](http://www.htslib.org/download) v>1.8<br/>
python2.7.5 <br/>

# Usage<br/>
`./junctionFilter.sh [inputDir] [outputDir] [refDir]`

inputDir: path to alignment output <br/>
outputDir: path for output files <br/>
refDir: path containing the reference files provided

# Output<br/>
1. Coverage.txt<br/>
2. Large_Deletions_NC_012920.1_106-16176.txt<br/>
3. Large_Deletions_NC_012920.1_356-15926.txt<br/>
4. Large_Deletions_NC_012920.1_No-Position-Filter.txt<br/>

# Versioning<br/>
Semantic Versioning 2.0.0

# License<br/>
MIT license

# Contact<br/>
**Brooke Hjelm**<br/>
bhjelm@usc.edu

**Michelle Webb**<br/>
michelgw@usc.edu
