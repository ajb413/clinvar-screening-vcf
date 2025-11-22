# ClinVar Genetic Screening Tool

A genetic screening tool that analyzes DNA variants from VCF files to identify carrier/potentially affected status associated with inherited diseases.

⚠️ **This tool is for informational purposes only. It is NOT a medical diagnosis. Always consult qualified healthcare professionals for medical decisions.**

### Process

1. Downloads a VCF file of known disease-causing genetic mutations from ClinVar
2. Reads the genetic data (VCF file from genetic testing)
3. Compares DNA against the database
4. Reports which disease-causing mutations subject carries

### What It Reports

The tool screens all genes in ClinVar but focuses on clinically significant, pathogenic variants. It filters out:
- Benign or likely benign variants
- Variants of uncertain significance
- Common variants (>1% population frequency)
- Low-quality submissions

Carrier Status, Affected Status, Potential Compound Heterozygote (Uncertain), Compound Heterozygote Affected

### Run

```bash
pip install requests
```

```bash
python clinvar_screening.py genetic_data.vcf
```

**Update the locally stored ClinVar database** (periodically):
```bash
python clinvar_screening.py genetic_data.vcf --update-database
```

### Example Output

```
CARRIER STATUS (2 findings):
----------------------------------------------------------------------
Gene: CFTR
Condition: Cystic fibrosis
Variant: chr7:117559590:G:A
Genotype: 0/1 (Heterozygous)
Clinical Significance: pathogenic
ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/variation/12345/

(VCF file line logged here)

```
