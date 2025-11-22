#!/usr/bin/env python3
"""
VCF Carrier Screening Script with Local ClinVar Database

Usage:
    python clinvar_screening.py sample.vcf
    python clinvar_screening.py sample.vcf --update-database
    python clinvar_screening.py sample.vcf --output-dir results/
"""

import os
import gzip
import shutil
import requests
import argparse
import logging
import time
from datetime import datetime
from typing import Dict, List

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ClinVarParsingError(Exception):
    pass

class VCFValidationError(Exception):
    pass

class ClinVarCarrierScreening:
    def __init__(self, database_dir: str = "./databases"):
        self.database_dir = database_dir
        self.clinvar_vcf_path = os.path.join(database_dir, "clinvar.vcf")
        self.clinvar_url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
        self.pathogenic_variants = {}
        
        os.makedirs(database_dir, exist_ok=True)
        
        # ClinVar classification terms to include
        self.pathogenic_terms = {
            'pathogenic',
            'likely pathogenic',
            'pathogenic/likely pathogenic'
        }
        
        # ClinVar classification terms to exclude
        self.exclude_terms = {
            'benign',
            'likely benign',
            'benign/likely benign',
            'uncertain significance',
            'not provided',
            'no classification for the single variant',
            'no classifications from unflagged records'
        }
        
        # Optional disease filtering (empty by default)
        self.excluded_disease_keywords = set()
        
        # Statistics for logging
        self.stats = {
            'total_processed': 0,
            'pathogenic_match': 0,
            'excluded_term': 0,
            'conflicting_filtered': 0,
            'high_frequency': 0,
            'somatic': 0,
            'disease_keyword': 0,
            'final_kept': 0
        }

    def check_database_age(self) -> bool:
        """Check if database needs updating based on age"""
        if not os.path.exists(self.clinvar_vcf_path):
            return True
        
        file_age_days = (time.time() - os.path.getmtime(self.clinvar_vcf_path)) / (24 * 3600)
        max_age_days = 30  # Maximum database age in days
        
        if file_age_days > max_age_days:
            logger.info(f"Database is {file_age_days:.1f} days old (max: {max_age_days})")
            return True
        return False

    def download_clinvar_database(self, force_update: bool = False):
        """Download and decompress ClinVar VCF database"""
        if os.path.exists(self.clinvar_vcf_path) and not force_update and not self.check_database_age():
            logger.info(f"Using existing database: {self.clinvar_vcf_path}")
            return
        
        logger.info(f"Downloading ClinVar database from {self.clinvar_url}")
        compressed_path = self.clinvar_vcf_path + ".gz"
        
        try:
            response = requests.get(self.clinvar_url, stream=True)
            response.raise_for_status()
            
            total_size = int(response.headers.get('content-length', 0))
            downloaded = 0
            
            with open(compressed_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total_size > 0 and downloaded % (10 * 1024 * 1024) == 0:
                            print(f"Progress: {(downloaded / total_size) * 100:.1f}%", end='\r')
            
            logger.info(f"Downloaded {downloaded / (1024*1024):.1f} MB")
            logger.info("Decompressing...")
            
            with gzip.open(compressed_path, 'rb') as f_in:
                with open(self.clinvar_vcf_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            os.remove(compressed_path)
            logger.info(f"Database ready: {self.clinvar_vcf_path}")
            
        except Exception as e:
            raise ClinVarParsingError(f"Error downloading database: {e}")

    def normalize_chromosome(self, chrom: str) -> str:
        """Add 'chr' prefix if missing"""
        return chrom if chrom.startswith('chr') else f"chr{chrom}"

    def normalize_variant(self, pos: int, ref: str, alt: str) -> tuple:
        """
        Normalize variant by trimming common suffix and prefix.
        Keeps at least one base in both ref and alt.
        Returns: (normalized_pos, normalized_ref, normalized_alt)
        """
        ref = ref.upper()
        alt = alt.upper()
        
        # Trim common suffix
        while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
            ref = ref[:-1]
            alt = alt[:-1]
        
        # Trim common prefix and adjust position
        while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
            ref = ref[1:]
            alt = alt[1:]
            pos += 1
        
        return (pos, ref, alt)

    def create_variant_key(self, chrom: str, pos: int, ref: str, alt: str) -> str:
        """Create standardized variant key for matching"""
        chrom = self.normalize_chromosome(chrom)
        pos, ref, alt = self.normalize_variant(pos, ref, alt)
        return f"{chrom}:{pos}:{ref}:{alt}"

    def parse_info_field(self, info: str) -> Dict[str, str]:
        """Parse VCF INFO field into dictionary"""
        info_dict = {}
        for item in info.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                info_dict[key] = value
        return info_dict

    def is_pathogenic_variant(self, clnsig: str, info_dict: Dict[str, str]) -> bool:
        """Determine if variant should be included based on clinical significance"""
        # Check for pathogenic terms
        if not any(term in clnsig for term in self.pathogenic_terms):
            return False
        
        self.stats['pathogenic_match'] += 1
        
        # Check for excluded terms
        if any(term in clnsig for term in self.exclude_terms):
            self.stats['excluded_term'] += 1
            return False
        
        # Handle conflicting classifications by counting votes
        if 'conflicting' in clnsig:
            clnsigconf = info_dict.get('CLNSIGCONF', '')
            if clnsigconf:
                pathogenic_votes = 0
                benign_votes = 0
                
                for item in clnsigconf.split('|'):
                    if '(' in item and ')' in item:
                        classification = item.split('(')[0].lower().replace('_', ' ')
                        try:
                            count = int(item.split('(')[1].split(')')[0])
                        except ValueError:
                            continue
                        
                        if any(term in classification for term in ['pathogenic', 'likely pathogenic']):
                            pathogenic_votes += count
                        elif any(term in classification for term in ['benign', 'likely benign']):
                            benign_votes += count
                
                if pathogenic_votes >= benign_votes:
                    return True
            
            self.stats['conflicting_filtered'] += 1
            return False
        
        return True

    def extract_population_frequency(self, info_dict: Dict) -> float:
        """Extract maximum population frequency from INFO field"""
        freq_fields = ['AF_ESP', 'AF_EXAC', 'AF_TGP', 'AF']
        max_freq = 0.0
        found = False
        
        for field in freq_fields:
            if field in info_dict:
                try:
                    freqs = [float(f) for f in info_dict[field].split(',') if f != '.']
                    if freqs:
                        max_freq = max(max_freq, max(freqs))
                        found = True
                except (ValueError, TypeError):
                    continue
        
        return max_freq if found else None

    def passes_filters(self, variant_data: Dict) -> bool:
        """Apply filtering criteria to variant"""
        # Review status filter
        review_status = variant_data['review_status'].lower()
        if any(s in review_status for s in ['no assertion criteria', 'no assertion provided', 'no classification']):
            return False
        
        # Population frequency filter (1% threshold)
        pop_freq = variant_data['population_frequency']
        if pop_freq is not None and pop_freq > 0.01:
            self.stats['high_frequency'] += 1
            return False
        
        # Origin filter (germline only)
        if 'somatic' in variant_data['origin'].lower():
            self.stats['somatic'] += 1
            return False
        
        # Disease keyword filter
        condition = variant_data['condition'].lower()
        if any(keyword in condition for keyword in self.excluded_disease_keywords):
            self.stats['disease_keyword'] += 1
            return False
        
        return True

    def parse_clinvar_database(self):
        """Parse ClinVar VCF and extract pathogenic variants"""
        logger.info("Parsing ClinVar database...")
        
        if not os.path.exists(self.clinvar_vcf_path):
            raise FileNotFoundError(f"Database not found: {self.clinvar_vcf_path}")
        
        with open(self.clinvar_vcf_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                self.stats['total_processed'] += 1
                
                if self.stats['total_processed'] % 500000 == 0:
                    logger.info(f"Processed {self.stats['total_processed']:,} variants...")
                
                try:
                    fields = line.strip().split('\t')
                    if len(fields) < 8:
                        continue
                    
                    chrom, pos, clinvar_id, ref, alt, info = fields[0], int(fields[1]), fields[2], fields[3], fields[4], fields[7]
                    
                    # Handle multi-allelic variants
                    alt_alleles = alt.split(',')
                    
                    for single_alt in alt_alleles:
                        variant_key = self.create_variant_key(chrom, pos, ref, single_alt)
                        
                        # Parse INFO field
                        info_dict = self.parse_info_field(info)
                        
                        # Check clinical significance
                        clnsig = info_dict.get('CLNSIG', '').lower().replace('_', ' ')
                        if not clnsig or not self.is_pathogenic_variant(clnsig, info_dict):
                            continue
                        
                        # Extract variant data
                        gene = info_dict.get('GENEINFO', '').split(':')[0] if 'GENEINFO' in info_dict else 'Unknown'
                        condition = info_dict.get('CLNDN', 'Unknown').replace('_', ' ')
                        variant_id = info_dict.get('CLNHGVS', 'Unknown')
                        review_status = info_dict.get('CLNREVSTAT', 'Unknown').replace('_', ' ')
                        origin = info_dict.get('ORIGIN', 'Unknown')
                        pop_freq = self.extract_population_frequency(info_dict)
                        
                        variant_data = {
                            'gene': gene,
                            'condition': condition,
                            'clinical_significance': clnsig,
                            'variant_id': variant_id,
                            'clinvar_variation_id': clinvar_id,
                            'review_status': review_status,
                            'origin': origin,
                            'population_frequency': pop_freq,
                            'chrom': self.normalize_chromosome(chrom),
                            'pos': pos,
                            'ref': ref,
                            'alt': single_alt
                        }
                        
                        # Apply additional filters
                        if not self.passes_filters(variant_data):
                            continue
                        
                        self.pathogenic_variants[variant_key] = variant_data
                        self.stats['final_kept'] += 1
                
                except Exception as e:
                    continue
        
        logger.info(f"Parsing complete:")
        logger.info(f"  Total processed: {self.stats['total_processed']:,}")
        logger.info(f"  Pathogenic kept: {self.stats['final_kept']:,}")

    def parse_vcf_file(self, vcf_path: str) -> List[Dict]:
        """Parse user VCF file and extract variants with genotypes"""
        logger.info(f"Parsing VCF: {vcf_path}")
        
        variants = []
        sample_index = None
        sample_name = None
        
        with open(vcf_path, 'r') as f:
            for line in f:
                if line.startswith('##'):
                    continue
                
                if line.startswith('#CHROM'):
                    header = line.strip().split('\t')
                    if len(header) <= 9:
                        raise VCFValidationError("VCF has no sample columns")
                    
                    # Always use first sample
                    sample_index = 9
                    sample_name = header[9]
                    logger.info(f"Analyzing sample: {sample_name}")
                    continue
                
                if sample_index is None:
                    continue
                
                try:
                    fields = line.strip().split('\t')
                    if len(fields) < sample_index + 1:
                        continue
                    
                    chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
                    format_field = fields[8]
                    sample_field = fields[sample_index]
                    
                    # Extract genotype
                    if 'GT' not in format_field:
                        continue
                    
                    gt_index = format_field.split(':').index('GT')
                    genotype = sample_field.split(':')[gt_index]
                    
                    # Skip reference genotypes
                    if genotype in ['.', './.', '0/0', '0|0']:
                        continue
                    
                    # Handle multi-allelic variants
                    alt_alleles = alt.split(',')
                    
                    for single_alt in alt_alleles:
                        variant_key = self.create_variant_key(chrom, pos, ref, single_alt)
                        variants.append({
                            'variant_key': variant_key,
                            'chrom': chrom,
                            'pos': pos,
                            'ref': ref,
                            'alt': single_alt,
                            'genotype': genotype,
                            'raw_vcf_line': line.strip()
                        })
                
                except (ValueError, IndexError):
                    continue
        
        logger.info(f"Found {len(variants):,} variants with non-reference genotypes")
        return variants

    def interpret_genotype(self, genotype: str) -> str:
        """Interpret genotype for carrier screening"""
        if genotype in ['0/1', '1/0', '0|1', '1|0', '0/2', '2/0', '0|2', '2|0', '1/2', '2/1', '1|2', '2|1']:
            return 'CARRIER'
        if genotype in ['1/1', '1|1', '2/2', '2|2']:
            return 'AFFECTED'
        return 'UNKNOWN'

    def is_phased_genotype(self, genotype: str) -> bool:
        """Check if genotype is phased (uses | instead of /)"""
        return '|' in genotype

    def check_trans_configuration(self, variants: List[Dict]) -> tuple:
        """
        Determine if variants are in trans (compound het) or cis (same chromosome).
        Returns: (is_trans, confidence, explanation)
        """
        # Check if all variants are phased
        all_phased = all(self.is_phased_genotype(v['genotype']) for v in variants)
        
        if not all_phased:
            return (None, 'unknown', 'Unphased genotypes - cannot determine cis/trans configuration')
        
        # For phased genotypes, check which chromosome each mutation is on
        # Genotype format: maternal|paternal
        # 1|0 = mutation on maternal chromosome
        # 0|1 = mutation on paternal chromosome
        
        maternal_mutations = []
        paternal_mutations = []
        
        for v in variants:
            gt = v['genotype']
            alleles = gt.split('|')
            
            # Check which allele is mutant (non-zero)
            if alleles[0] != '0':
                maternal_mutations.append(v['pos'])
            if alleles[1] != '0':
                paternal_mutations.append(v['pos'])
        
        # Trans configuration: mutations on different chromosomes
        # (some on maternal, some on paternal)
        has_maternal = len(maternal_mutations) > 0
        has_paternal = len(paternal_mutations) > 0
        
        if has_maternal and has_paternal:
            return (True, 'high', 'Phased genotypes show variants on different chromosomes (trans)')
        elif has_maternal and not has_paternal:
            return (False, 'high', 'Phased genotypes show variants on same chromosome - maternal (cis)')
        elif has_paternal and not has_maternal:
            return (False, 'high', 'Phased genotypes show variants on same chromosome - paternal (cis)')
        else:
            return (None, 'low', 'Unable to determine configuration from phasing')

    def detect_compound_heterozygotes(self, findings: List[Dict]) -> List[Dict]:
        """
        Detect potential compound heterozygotes (multiple carrier variants in same gene).
        Properly handles phasing to distinguish trans (affected) from cis (carrier).
        """
        gene_variants = {}
        
        for finding in findings:
            gene = finding['gene']
            if gene not in gene_variants:
                gene_variants[gene] = []
            gene_variants[gene].append(finding)
        
        updated_findings = []
        
        for gene, variants in gene_variants.items():
            if len(variants) >= 2:
                logger.info(f"Checking {gene}: found {len(variants)} variants")
                
                carriers = [v for v in variants if v['interpretation'] == 'CARRIER']
                unique_positions = set(v['pos'] for v in carriers)
                
                logger.info(f"  {len(carriers)} are CARRIER status at {len(unique_positions)} unique positions")
                
                # Log each variant's details
                for v in variants:
                    logger.info(f"    Position {v['pos']}: {v['genotype']} - {v['interpretation']}")
                
                # Only consider if we have multiple carriers at different positions
                if len(carriers) >= 2 and len(unique_positions) >= 2:
                    logger.info(f"  Analyzing phasing for {gene}...")
                    
                    # Check phasing to determine if trans or cis
                    is_trans, confidence, explanation = self.check_trans_configuration(carriers)
                    
                    logger.info(f"  Phasing result: is_trans={is_trans}, confidence={confidence}")
                    logger.info(f"  Explanation: {explanation}")
                    
                    if is_trans is True:
                        # Confirmed trans - compound heterozygote (affected)
                        # Add ALL variants with updated interpretation
                        for carrier in carriers:
                            updated_variant = carrier.copy()
                            updated_variant['interpretation'] = 'COMPOUND_HET_AFFECTED'
                            updated_variant['compound_het_info'] = f"Part of compound heterozygote (trans): {len(carriers)} variants in {gene}. {explanation}"
                            updated_variant['phasing_confidence'] = confidence
                            updated_findings.append(updated_variant)
                        logger.info(f"  → Classified as COMPOUND_HET_AFFECTED ({len(carriers)} variants)")
                    elif is_trans is False:
                        # Confirmed cis - double carrier (not affected)
                        # Add ALL variants as regular carriers with note
                        for carrier in carriers:
                            updated_variant = carrier.copy()
                            updated_variant['interpretation'] = 'CARRIER'
                            updated_variant['compound_het_info'] = f"Multiple variants in {gene} on same chromosome (cis): carrier only. {explanation}"
                            updated_variant['phasing_confidence'] = confidence
                            updated_findings.append(updated_variant)
                        logger.info(f"  → Classified as CARRIER (cis) - {len(carriers)} variants")
                    else:
                        # Unknown - unphased data
                        # Add ALL variants as potential compound het
                        for carrier in carriers:
                            updated_variant = carrier.copy()
                            updated_variant['interpretation'] = 'POTENTIAL_COMPOUND_HET'
                            updated_variant['compound_het_info'] = f"UNCERTAIN: Part of {len(carriers)} variants in {gene}. {explanation}. Recommend phasing or parental testing."
                            updated_variant['phasing_confidence'] = confidence
                            updated_findings.append(updated_variant)
                        logger.info(f"  → Classified as POTENTIAL_COMPOUND_HET (unphased) - {len(carriers)} variants")
                else:
                    logger.info(f"  Not enough carriers at different positions for compound het analysis")
                    for variant in variants:
                        variant['compound_het_info'] = f"One of {len(variants)} variants in {gene}"
                        variant['phasing_confidence'] = 'n/a'
                    updated_findings.extend(variants)
            else:
                # Single variant in gene
                variants[0]['compound_het_info'] = f"Single variant in {gene}"
                variants[0]['phasing_confidence'] = 'n/a'
                updated_findings.extend(variants)
        
        return updated_findings

    def analyze_variants(self, vcf_variants: List[Dict]) -> List[Dict]:
        """Match VCF variants against pathogenic ClinVar variants"""
        logger.info("Analyzing variants...")
        
        findings = []
        
        for variant in vcf_variants:
            variant_key = variant['variant_key']
            
            if variant_key in self.pathogenic_variants:
                clinvar_data = self.pathogenic_variants[variant_key]
                interpretation = self.interpret_genotype(variant['genotype'])
                
                if interpretation in ['CARRIER', 'AFFECTED']:
                    finding = {
                        'variant_key': variant_key,
                        'chrom': variant['chrom'],
                        'pos': variant['pos'],
                        'ref': variant['ref'],
                        'alt': variant['alt'],
                        'genotype': variant['genotype'],
                        'raw_vcf_line': variant['raw_vcf_line'],
                        'gene': clinvar_data['gene'],
                        'condition': clinvar_data['condition'],
                        'clinical_significance': clinvar_data['clinical_significance'],
                        'variant_id': clinvar_data['variant_id'],
                        'clinvar_variation_id': clinvar_data['clinvar_variation_id'],
                        'review_status': clinvar_data['review_status'],
                        'population_frequency': clinvar_data['population_frequency'],
                        'interpretation': interpretation,
                        'compound_het_info': ''
                    }
                    findings.append(finding)
        
        logger.info(f"Found {len(findings)} pathogenic variants")
        
        # Detect compound heterozygotes
        findings = self.detect_compound_heterozygotes(findings)
        
        return findings

    def output_results(self, findings: List[Dict], output_dir: str, vcf_path: str):
        """Display and save analysis results"""
        vcf_basename = os.path.splitext(os.path.basename(vcf_path))[0] if vcf_path else "unknown"
        
        if not findings:
            print("\n" + "="*80)
            print("CARRIER SCREENING RESULTS")
            print("="*80)
            print("No pathogenic variants found.")
            return
        
        # Categorize findings
        carriers = [f for f in findings if f['interpretation'] == 'CARRIER']
        affected = [f for f in findings if f['interpretation'] == 'AFFECTED']
        compound_het_affected = [f for f in findings if f['interpretation'] == 'COMPOUND_HET_AFFECTED']
        potential_compound_het = [f for f in findings if f['interpretation'] == 'POTENTIAL_COMPOUND_HET']
        
        # Display results
        print("\n" + "="*80)
        print("CARRIER SCREENING RESULTS")
        print("="*80)
        print(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Total Findings: {len(findings)}")
        print(f"  Carriers: {len(carriers)}")
        print(f"  Affected: {len(affected)}")
        print(f"  Compound Het (Affected): {len(compound_het_affected)}")
        print(f"  Potential Compound Het (Uncertain): {len(potential_compound_het)}")
        
        if carriers:
            print(f"\nCARRIER STATUS ({len(carriers)} findings):")
            print("-" * 70)
            for f in carriers:
                print(f"Gene: {f['gene']}")
                print(f"Condition: {f['condition']}")
                print(f"Variant: {f['variant_key']}")
                print(f"Genotype: {f['genotype']} (Heterozygous)")
                print(f"Clinical Significance: {f['clinical_significance']}")
                if f['population_frequency']:
                    print(f"Population Frequency: {f['population_frequency']:.4f}")
                if f.get('compound_het_info'):
                    print(f"Note: {f['compound_het_info']}")
                if f['clinvar_variation_id'] != 'Unknown':
                    print(f"ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/variation/{f['clinvar_variation_id']}/")
                print()
        
        if potential_compound_het:
            print(f"\nPOTENTIAL COMPOUND HETEROZYGOTES - UNCERTAIN ({len(potential_compound_het)} findings):")
            print("-" * 70)
            print("⚠️  Multiple variants in same gene detected, but phasing unknown.")
            print("   Cannot determine if variants are on same chromosome (cis - carrier)")
            print("   or different chromosomes (trans - affected).")
            print("   Recommendation: Phasing analysis or parental testing needed.\n")
            for f in potential_compound_het:
                print(f"Gene: {f['gene']}")
                print(f"Condition: {f['condition']}")
                print(f"Status: {f['compound_het_info']}")
                print(f"Phasing Confidence: {f.get('phasing_confidence', 'unknown')}")
                if f['clinvar_variation_id'] != 'Unknown':
                    print(f"ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/variation/{f['clinvar_variation_id']}/")
                print()
        
        if compound_het_affected:
            print(f"\nCOMPOUND HETEROZYGOTES - AFFECTED ({len(compound_het_affected)} findings):")
            print("-" * 70)
            print("⚠️  Variants confirmed on different chromosomes (trans configuration).")
            print("   This indicates affected status for recessive conditions.\n")
            for f in compound_het_affected:
                print(f"Gene: {f['gene']}")
                print(f"Condition: {f['condition']}")
                print(f"Info: {f['compound_het_info']}")
                print(f"Phasing Confidence: {f.get('phasing_confidence', 'unknown')}")
                if f['clinvar_variation_id'] != 'Unknown':
                    print(f"ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/variation/{f['clinvar_variation_id']}/")
                print()
        
        if affected:
            print(f"\nHOMOZYGOUS AFFECTED STATUS ({len(affected)} findings):")
            print("-" * 70)
            print("⚠️  Homozygous pathogenic variants detected.")
            print("   Please consult with a genetic counselor.\n")
            for f in affected:
                print(f"Gene: {f['gene']}")
                print(f"Condition: {f['condition']}")
                print(f"Variant: {f['variant_key']}")
                print(f"Genotype: {f['genotype']} (Homozygous)")
                print(f"Clinical Significance: {f['clinical_significance']}")
                if f['clinvar_variation_id'] != 'Unknown':
                    print(f"ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/variation/{f['clinvar_variation_id']}/")
                print()
        
        # Save summary file
        self.save_summary(findings, carriers, affected, compound_het_affected, potential_compound_het, output_dir, vcf_basename)

    def save_summary(self, findings, carriers, affected, compound_het_affected, potential_compound_het, output_dir, vcf_basename):
        """Save summary results to text file"""
        summary_file = os.path.join(output_dir, f'carrier_screening_summary_{vcf_basename}.txt')
        
        with open(summary_file, 'w') as f:
            f.write("CARRIER SCREENING ANALYSIS SUMMARY\n")
            f.write("="*80 + "\n\n")
            f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total Findings: {len(findings)}\n")
            f.write(f"  Carriers: {len(carriers)}\n")
            f.write(f"  Affected (Homozygous): {len(affected)}\n")
            f.write(f"  Compound Het Affected (Trans): {len(compound_het_affected)}\n")
            f.write(f"  Potential Compound Het (Uncertain): {len(potential_compound_het)}\n\n")
            
            if carriers:
                f.write(f"CARRIER STATUS ({len(carriers)} findings):\n")
                f.write("-" * 70 + "\n")
                for finding in carriers:
                    f.write(f"Gene: {finding['gene']}\n")
                    f.write(f"Condition: {finding['condition']}\n")
                    f.write(f"Variant: {finding['variant_key']}\n")
                    f.write(f"Genotype: {finding['genotype']}\n")
                    f.write(f"Clinical Significance: {finding['clinical_significance']}\n")
                    if finding.get('compound_het_info'):
                        f.write(f"Note: {finding['compound_het_info']}\n")
                    if finding['clinvar_variation_id'] != 'Unknown':
                        f.write(f"ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/variation/{finding['clinvar_variation_id']}/\n")
                    f.write(f"\n{finding['raw_vcf_line']}\n\n")
            
            if potential_compound_het:
                f.write(f"\nPOTENTIAL COMPOUND HETEROZYGOTES - UNCERTAIN ({len(potential_compound_het)} findings):\n")
                f.write("-" * 70 + "\n")
                f.write("⚠️  Multiple variants in same gene detected, but phasing unknown.\n")
                f.write("   Cannot determine if cis (carrier) or trans (affected).\n")
                f.write("   Recommendation: Phasing analysis or parental testing.\n\n")
                for finding in potential_compound_het:
                    f.write(f"Gene: {finding['gene']}\n")
                    f.write(f"Condition: {finding['condition']}\n")
                    f.write(f"Status: {finding['compound_het_info']}\n")
                    f.write(f"Phasing Confidence: {finding.get('phasing_confidence', 'unknown')}\n")
                    if finding['clinvar_variation_id'] != 'Unknown':
                        f.write(f"ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/variation/{finding['clinvar_variation_id']}/\n")
                    f.write(f"\n{finding['raw_vcf_line']}\n\n")
            
            if compound_het_affected:
                f.write(f"\nCOMPOUND HETEROZYGOTES - AFFECTED ({len(compound_het_affected)} findings):\n")
                f.write("-" * 70 + "\n")
                f.write("⚠️  Variants confirmed on different chromosomes (trans).\n\n")
                for finding in compound_het_affected:
                    f.write(f"Gene: {finding['gene']}\n")
                    f.write(f"Condition: {finding['condition']}\n")
                    f.write(f"Info: {finding['compound_het_info']}\n")
                    f.write(f"Phasing Confidence: {finding.get('phasing_confidence', 'unknown')}\n")
                    if finding['clinvar_variation_id'] != 'Unknown':
                        f.write(f"ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/variation/{finding['clinvar_variation_id']}/\n")
                    f.write(f"\n{finding['raw_vcf_line']}\n\n")
            
            if affected:
                f.write(f"\nHOMOZYGOUS AFFECTED STATUS ({len(affected)} findings):\n")
                f.write("-" * 70 + "\n")
                for finding in affected:
                    f.write(f"Gene: {finding['gene']}\n")
                    f.write(f"Condition: {finding['condition']}\n")
                    f.write(f"Variant: {finding['variant_key']}\n")
                    f.write(f"Genotype: {finding['genotype']}\n")
                    f.write(f"Clinical Significance: {finding['clinical_significance']}\n")
                    if finding['clinvar_variation_id'] != 'Unknown':
                        f.write(f"ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/variation/{finding['clinvar_variation_id']}/\n")
                    f.write(f"\n{finding['raw_vcf_line']}\n\n")
            
            if not findings:
                f.write("No pathogenic variants found.\n")
            
            # Gene summary
            f.write("\nGENES BY CATEGORY:\n")
            f.write("-"*30 + "\n")
            f.write(f"CARRIER: {', '.join(sorted(set(f['gene'] for f in carriers))) or 'None'}\n")
            f.write(f"AFFECTED (Homozygous): {', '.join(sorted(set(f['gene'] for f in affected))) or 'None'}\n")
            f.write(f"COMPOUND_HET_AFFECTED (Trans): {', '.join(sorted(set(f['gene'] for f in compound_het_affected))) or 'None'}\n")
            f.write(f"POTENTIAL_COMPOUND_HET (Uncertain): {', '.join(sorted(set(f['gene'] for f in potential_compound_het))) or 'None'}\n")
        
        print(f"\nSummary saved to: {summary_file}")

    def run_analysis(self, vcf_path: str, force_update: bool = False, output_dir: str = "."):
        """Run complete carrier screening analysis"""
        logger.info("Starting carrier screening analysis...")
        
        # Setup database
        self.download_clinvar_database(force_update)
        self.parse_clinvar_database()
        
        # Parse and analyze VCF
        vcf_variants = self.parse_vcf_file(vcf_path)
        findings = self.analyze_variants(vcf_variants)
        
        # Output results
        self.output_results(findings, output_dir, vcf_path)
        
        return findings


def main():
    parser = argparse.ArgumentParser(
        description='VCF Carrier Screening Tool with Local ClinVar Database',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python clinvar_screening.py sample.vcf
  python clinvar_screening.py sample.vcf --update-database
  python clinvar_screening.py sample.vcf --output-dir results/
        """
    )
    
    parser.add_argument('vcf_file', help='Path to VCF file for analysis')
    parser.add_argument('--database-dir', default='./databases', 
                       help='Database directory (default: ./databases)')
    parser.add_argument('--output-dir', default='.', 
                       help='Output directory (default: current directory)')
    parser.add_argument('--update-database', action='store_true',
                       help='Force update of ClinVar database')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    try:
        screener = ClinVarCarrierScreening(args.database_dir)
        screener.run_analysis(
            vcf_path=args.vcf_file,
            force_update=args.update_database,
            output_dir=args.output_dir
        )
        logger.info("Analysis completed successfully!")
        return 0
        
    except (ClinVarParsingError, VCFValidationError, FileNotFoundError) as e:
        logger.error(f"Analysis failed: {e}")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return 1


if __name__ == "__main__":
    exit(main())