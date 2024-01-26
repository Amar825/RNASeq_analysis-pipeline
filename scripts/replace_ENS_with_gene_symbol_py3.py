import fileinput

# First create the "Gene ID => Gene Symbol" look-up table:
#   Read the file containing Ensembl IDs and gene names line by line.
#   If the line contains TAB-separable columns, associate the gene ID in the first column with the gene name in the second 
genesymbols = {} 
with open("ENS_to_gene-names.txt") as genesfile:
    for line in genesfile.readlines():
        if ('\t' in line):
            id,name = line.strip().split('\t')
            genesymbols[id] = name

# Read a file line by line (from STDIN) and split each line into TAB-separated columns.
# Insert an extra column at the beginning for the gene symbol (by making a copy of the previous first column containing the gene ID).
# If the value in the first column is recognized as an ENSEMBL ID in my look-up table, replace it with the gene symbol.
# If not, then check if this is the header (first line) which must be fixed in an ad-hoc way by adding yet another column
n=0
for line in fileinput.input():
    n=n+1
    columns = line.strip().split('\t')
    columns.insert(0,columns[0])  # duplicate the first column
    if columns[0] in genesymbols:
        columns[0] = genesymbols[columns[0]] # replace ENSEMBL ID in first column with gene symbol
    else:
        if (n==1): # We assume that this is the header row
            columns[0]="gene_id"
            columns.insert(0,"gene_symbol")
    line = '\t'.join(columns)
    print(line)
