# Project summary
Creating a command line tool to extract GlyGen O-GlcNAc and Phosphorylation sites into an ML ready dataset.

# Installation Instruction

To run the program:

- Python must be downloaded on the computer to run the program, please follow instructions on: <https://www.python.org/downloads/>
- In the terminal window run:
  - If the version is Python3 or above:
    - Pip3 install requests
    - Pip3 install pandas
  - If the version is older:
    - Pip install requests
    - Pip install pandas

# Running Instructions
First, ensure that your terminal is in the directory of the program folder (usually within the downloads folder)
To run the program use the following command line:

- If the version is Python3 or above:
  - Python3 main.py (for the initial run)
  - If you would like to restart with no cache add “-n” or “--nocache” after the initial command line (eg. python3 main.py --nocache)
  - The file outputs a tsv by default, if you require a csv, use the command line -o csv after the initial command line (eg. python3 main.py -o csv)
  - The program calls on the production API by default, if you require usage of the beta API use the command line -a beta after the initial command line (eg. python3 main.py -a beta)
  - If you want a combination of the above requirements, they can be added one after another (eg. python3 main.py --nocache -o csv -a beta)
- If the version is an older Python version:
  - Replace python3 with python in the command line (eg. python main.py --nocache -o csv -a beta)

# Output Documentation

| field | description | example | value |
| --- | --- | --- | --- |
| uniprotkb_ac | Canonical Accession of the protein from UniProtKB | P00533-1 | Single |
| protein_name | Name of the protein from UniProtKB | Epidermal growth factor receptor | Single or no value |
| gene_name | HGNC Gene name for the gene coding protein | EGFR | Single or no value |
| site | Position on the protein for the Phosphorylation/ O-linked site that is known to be modified or is unknown to be modified. | 430 | Single (remove entries with no values) |
| amino_acid | Amino acid (aa) that is modified in three letter code | Thr | Single (remove entries with no values) |
| phosphorylation | Is the site O-linked phosphorylated | Y or N | Single |
| glycosylation | Is the site O-linked glycosylated | Y or N | Single |
| peptide_seq_five_before | Peptide sequence 5 aa before the modified residue (Upstream towards N terminal) | IIRGR | Single or no value |
| peptide_seqfive_after | Peptide sequence 5 aa after the modified residue (Downstream towards C terminal) | KQHGQ | Single or no value |
| peptide_seq_ten_before | Peptide sequence 10 aa before the modified residue (towards N terminal) | FENLEIIRGR | Single or no value |
| peptide_seqten_after | Peptide sequence 10 aa after the modified residue (towards C terminal) | KQHGQFSLAV | Single or no value |
| kinase_gene_name | Names of kinase genes that contribute to the phosphorylation | PIM1\|PIM2\|AKT1 | Single or no value or multiple value |
| organism | Scientific name of the organism to which the protein belongs | Homo sapiens | Single |
| taxonomy_id | NCBI Taxonomy ID for the scientific name of the organism | 9606 | Single |
| molecular_function | Gene Ontology Molecular Function name and ID | ATPase binding (GO:0051117)\|ATP binding (GO:0005524) | Single or no value or multiple value |
| biological_process | Gene Ontology Biological Function name and ID | activation of phospholipase C activity (GO:0007202)\|astrocyte activation (GO:0048143) | Single or no value or multiple value |
| cellular_component | Gene Ontology Cellular Component name and ID | apical plasma membrane (GO:0016324)\|basal plasma membrane (GO:0009925) | Single or no value or multiple value |
| domain | Domain Description of the | Protein kinase | Single or no value |
| range | Range of Amino Acids within the Domain | 150-408 | Single or no value |
| status | Whether the site is known or unknown to be glycosylated | Known_glycosite<br><br>unknown_glycosite | Single or no value |

# Examples
