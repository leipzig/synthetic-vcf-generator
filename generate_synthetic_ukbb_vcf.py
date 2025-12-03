#!/usr/bin/env python3
"""Generate a single-sample exome VCF with UKBB model fields"""

import sys
import uuid
import multiprocessing
import subprocess
import random
import csv
from pathlib import Path
from collections import defaultdict
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    def tqdm(iterable, **kwargs):
        return iterable
from synthetic_vcf_generator.vcf_model import VCFModel
from synthetic_vcf_generator.virtual_vcf import VirtualVCF
from synthetic_vcf_generator.bed_parser import BEDIntervals

def generate_ukbb_vcf(
    model_vcf_path: str,
    output_path: str,
    bed_file_path: str,
    num_rows: int = None,
    seed: int = None,
    sample_uuid: str = None,
    sex: str = None,
):
    """
    Generate a single-sample exome VCF using UKBB model fields.
    
    Args:
        model_vcf_path: Path to UKBB model VCF file
        output_path: Path to output VCF file
        bed_file_path: Path to BED file defining exome regions
        num_rows: Number of variants to generate per chromosome (default: auto-calculated)
        seed: Random seed for reproducibility
        sample_uuid: UUID for sample name (if None, generates a new one)
        sex: Sex of sample ("male" or "female", if None, randomly assigned)
    
    Returns:
        tuple: (variant_count, sample_uuid, sex, variant_counts_by_chrom, is_exome)
    """
    # Check if output should be suppressed (for batch processing)
    suppress_output = hasattr(generate_ukbb_vcf, '_suppress_output') and generate_ukbb_vcf._suppress_output
    
    # Parse the BED file to get chromosomes and intervals
    if not suppress_output:
        print(f"Loading BED file from {bed_file_path}")
    bed_intervals = BEDIntervals(bed_file_path)
    chromosomes = sorted(bed_intervals.get_chromosomes())
    if not suppress_output:
        print(f"Found {len(chromosomes)} chromosomes: {', '.join(chromosomes[:5])}...")
    
    # Calculate total exome size
    total_bases = bed_intervals.get_total_bases()
    if not suppress_output:
        print(f"Total exome coverage: {total_bases:,} bases")
    
    # Auto-calculate num_rows if not provided
    # Aim for roughly 1 variant per 1000 bases in exome regions
    if num_rows is None:
        num_rows = max(100, total_bases // 1000)
        if not suppress_output:
            print(f"Auto-calculated {num_rows} variants per chromosome")
    
    # Parse the UKBB model
    model = VCFModel(model_vcf_path)
    
    # Generate UUID for sample name if not provided
    if sample_uuid is None:
        sample_uuid = str(uuid.uuid4())
    
    # Assign sex if not provided (use UUID hash for independent randomness)
    if sex is None:
        # Use UUID hash to ensure each sample gets independent sex assignment
        # This ensures variety even when seeds are similar
        uuid_hash = hash(sample_uuid)
        sex_random = random.Random(uuid_hash)
        sex = sex_random.choice(["male", "female"])
    
    # Determine if this is an exome (has BED file) or genome
    is_exome = bed_file_path is not None
    
    # Track variant counts per chromosome
    variant_counts_by_chrom = defaultdict(int)
    if not suppress_output:
        print(f"Generating exome VCF with {num_rows} variants per chromosome")
    vcf = VirtualVCF(
        num_rows=num_rows,
        num_samples=1,
        chromosomes=chromosomes,
        sample_prefix="",
        id_type="uuid",
        random_seed=seed,
        phased=True,
        large_format=True,
        bed_intervals=bed_intervals,
    )
    
    # Generate the VCF with custom header
    if not suppress_output:
        print(f"Writing to {output_path}")
    with open(output_path, 'w') as f:
        # Write fileformat
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##source=VirtualVCF-UKBB-Model 0.1.0\n")
        
        # Write model header lines (FORMAT and INFO)
        for line in model.get_header_lines():
            f.write(line + "\n")
        
        # Write contig lines
        for chrom in chromosomes:
            f.write(f"##contig=<ID={chrom}>\n")
        
        # Write column header
        format_string = model.get_format_string()
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_uuid}\n")
        
        # Generate and write variant rows
        variant_count = 0
        skipped_count = 0
        for line in vcf:
            if line.startswith("#"):
                continue  # Skip header lines from VirtualVCF
            
            # Parse the line to extract basic fields
            parts = line.strip().split('\t')
            if len(parts) < 10:
                skipped_count += 1
                continue
            
            # Check if position is in BED intervals (VirtualVCF should filter, but double-check)
            chrom = parts[0]
            pos = int(parts[1])
            if not bed_intervals.contains_position(chrom, pos):
                skipped_count += 1
                continue
            
            # Sex-specific filtering for chrX and chrY
            if chrom == "chrY" and sex == "female":
                # Females should not have chrY variants
                skipped_count += 1
                continue
            if chrom == "chrX" and sex == "male":
                # Males have only one X chromosome, so skip homozygous ref variants
                # We'll adjust genotypes below
                pass
            
            var_id = parts[2]
            ref = parts[3]
            alt = parts[4]
            qual = parts[5]
            filt = parts[6]
            info = parts[7]
            fmt = parts[8]
            sample_data = parts[9] if len(parts) > 9 else ""
            
            # Generate UKBB-style INFO fields
            # Extract some values from existing info or generate new ones
            info_fields = {}
            for field in info.split(';'):
                if '=' in field:
                    k, v = field.split('=', 1)
                    info_fields[k] = v
            
            # Build UKBB INFO string with all required fields
            ukbb_info = []
            
            # AC, AF, AN - allele count/frequency
            if 'AF' in info_fields:
                af_val = float(info_fields['AF'])
                ac_val = int(af_val * 2)  # For single sample, AN=2
                ukbb_info.append(f"AC={ac_val}")
                ukbb_info.append(f"AF={af_val:.2f}")
            else:
                ukbb_info.append("AC=1")
                ukbb_info.append("AF=0.50")
            ukbb_info.append("AN=2")
            
            # DP - depth
            if 'DP' in info_fields:
                ukbb_info.append(f"DP={info_fields['DP']}")
            else:
                ukbb_info.append("DP=80")
            
            # Quality metrics
            ukbb_info.append("MQ=58.2")
            ukbb_info.append("MQRankSum=0.72")
            ukbb_info.append("QD=25.3")
            ukbb_info.append("ReadPosRankSum=0.15")
            ukbb_info.append("FS=3.2")
            ukbb_info.append("SOR=1.1")
            
            # UKBB-specific fields
            ukbb_info.append(f"ALLELE_ID=A{variant_count+1}")
            ukbb_info.append("TARGETED")
            
            # HML - heterozygote model (het, hom, hom_alt)
            # Adjust genotype for sex-specific chromosomes
            gt = sample_data.split(':')[0] if sample_data else "0/0"
            
            # For males on chrX: convert to hemizygous (only one allele)
            # Males have only one X chromosome, so they can't be heterozygous
            if chrom == "chrX" and sex == "male":
                if gt in ["0|0", "0/0"]:
                    gt = "0"  # Hemizygous reference
                elif gt in ["1|1", "1/1"]:
                    gt = "1"  # Hemizygous alternate
                elif gt in ["0|1", "1|0", "0/1", "1/0"]:
                    # Males can't be heterozygous on X - convert to hemizygous alternate
                    # (more realistic than reference since X-inactivation doesn't apply to males)
                    gt = "1"
            
            # For males on chrY: should be hemizygous (only one allele)
            if chrom == "chrY" and sex == "male":
                if gt in ["0|0", "0/0"]:
                    gt = "0"
                elif gt in ["1|1", "1/1"]:
                    gt = "1"
                elif gt in ["0|1", "1|0", "0/1", "1/0"]:
                    # chrY should rarely be heterozygous, convert to homozygous alt
                    gt = "1"
            
            if gt in ["0/0", "0|0", "0"]:
                hml = "hom"
            elif gt in ["1/1", "1|1", "1"]:
                hml = "hom_alt"
            else:
                hml = "het"
            ukbb_info.append(f"HML={hml}")
            
            ukbb_info.append("JIDS=S1")
            ukbb_info.append("MOSAIC=0.02")
            ukbb_info.append("SoftClipRatio=0.01")
            ukbb_info.append("FractionInformativeReads=0.95")
            
            # Optional flags
            if variant_count % 3 == 0:
                ukbb_info.append("Recombinant")
            
            ukbb_info.append(f"PS={100 + variant_count * 10}")
            ukbb_info.append(f"EVENT=E{variant_count+1}")
            
            # EVENTTYPE
            if len(alt) > len(ref):
                event_type = "INS"
            elif len(alt) < len(ref):
                event_type = "DEL"
            else:
                event_type = "SNV"
            ukbb_info.append(f"EVENTTYPE={event_type}")
            
            # Generate UKBB-style FORMAT data
            # Use the format string from model
            format_fields = format_string.split(':')
            
            # Generate values for each format field
            format_values = []
            for field in format_fields:
                if field == "GT":
                    # Ensure GT is in proper format (use | for phased, / for unphased)
                    # For hemizygous (single allele), use | notation
                    if gt in ["0", "1"]:
                        format_values.append(gt)  # Hemizygous: single allele
                    else:
                        format_values.append(gt)  # Diploid: already formatted
                elif field == "AD":
                    # Allelic depths
                    if gt in ["0/0", "0|0", "0"]:
                        format_values.append("45,5")
                    elif gt in ["1/1", "1|1", "1"]:
                        format_values.append("5,45")
                    else:
                        format_values.append("35,35")
                elif field == "DP":
                    format_values.append("80")
                elif field == "AF":
                    if gt in ["0/0", "0|0", "0"]:
                        format_values.append("0.93,0.07")
                    elif gt in ["1/1", "1|1", "1"]:
                        format_values.append("0.07,0.93")
                    else:
                        format_values.append("0.44,0.56")
                elif field == "GQ":
                    format_values.append("99")
                elif field == "PL":
                    if gt in ["0/0", "0|0", "0"]:
                        format_values.append("0,120,890")
                    elif gt in ["1/1", "1|1", "1"]:
                        format_values.append("1120,145,0")
                    else:
                        format_values.append("0,120,890")
                elif field == "GP":
                    format_values.append("10,0,990")
                elif field == "F1R2":
                    format_values.append("32,28")
                elif field == "F2R1":
                    format_values.append("25,20")
                elif field == "SB":
                    format_values.append("0,1,39,20")
                elif field == "MB":
                    format_values.append("0.12")
                elif field == "SQ":
                    format_values.append("28.5")
                elif field == "PRI":
                    format_values.append("50")
                elif field == "PS":
                    format_values.append(str(100 + variant_count * 10))
                elif field == "QL":
                    format_values.append("0.88")
                elif field.startswith("J"):
                    # Joint fields - simplified
                    if field == "JAD":
                        format_values.append("70,50")
                    elif field == "JAF":
                        format_values.append("0.15,0.32")
                    elif field == "JDP":
                        format_values.append("120")
                    elif field == "JGQ":
                        format_values.append("98")
                    elif field == "JGT":
                        format_values.append(gt)
                    elif field == "JPL":
                        format_values.append("0,145,1120")
                    elif field == "JQL":
                        format_values.append("0.88")
                    elif field == "JVQL":
                        format_values.append("0.88")
                else:
                    format_values.append(".")
            
            # Write the variant line
            info_str = ";".join(ukbb_info)
            format_data = ":".join(format_values)
            f.write(f"{chrom}\t{pos}\t{var_id}\t{ref}\t{alt}\t{qual}\t{filt}\t{info_str}\t{format_string}\t{format_data}\n")
            variant_count += 1
            variant_counts_by_chrom[chrom] += 1
    
    if not suppress_output:
        print(f"✓ Generated {variant_count} variants in {output_path}")
        if skipped_count > 0:
            print(f"  (Skipped {skipped_count} variants outside exome regions)")
    
    return variant_count, sample_uuid, sex, dict(variant_counts_by_chrom), is_exome


def _generate_single_vcf(args):
    """Worker function for multiprocessing batch generation."""
    model_path, bed_path, output_dir, num_rows, seed, index = args
    sample_uuid = str(uuid.uuid4())
    output_path = Path(output_dir) / f"{sample_uuid}.vcf"
    
    try:
        # Suppress print statements in worker processes
        generate_ukbb_vcf._suppress_output = True
        
        # Use index as seed for variety and sex assignment
        # If base seed is provided, combine with index for unique seeds
        sample_seed = (seed + index) if seed is not None else index
        
        variant_count, _, sex, variant_counts_by_chrom, is_exome = generate_ukbb_vcf(
            model_vcf_path=model_path,
            output_path=str(output_path),
            bed_file_path=bed_path,
            num_rows=num_rows,
            seed=sample_seed,
            sample_uuid=sample_uuid,
            sex=None,  # Will be randomly assigned based on UUID hash
        )
        
        # Compress with bgzip
        try:
            subprocess.run(
                ['bgzip', '-f', str(output_path)],
                capture_output=True,
                check=True
            )
            compressed_path = output_path.with_suffix('.vcf.gz')
            
            # Create tabix index
            subprocess.run(
                ['tabix', '-p', 'vcf', str(compressed_path)],
                capture_output=True,
                check=True
            )
        except subprocess.CalledProcessError as e:
            return (False, sample_uuid, None, None, None, None, sample_seed, f"Compression/indexing error: {e.stderr}")
        except Exception as e:
            return (False, sample_uuid, None, None, None, None, sample_seed, f"Compression/indexing exception: {str(e)}")
        
        return (True, sample_uuid, variant_count, sex, variant_counts_by_chrom, is_exome, sample_seed, None)
    except Exception as e:
        return (False, sample_uuid, None, None, None, None, None, str(e))


def generate_batch_ukbb_vcf(
    model_vcf_path: str,
    output_dir: str,
    bed_file_path: str,
    num_vcfs: int = 100000,
    num_rows: int = None,
    seed: int = None,
    num_threads: int = None,
):
    """
    Generate multiple single-sample exome VCF files, each named by sample UUID.
    
    Args:
        model_vcf_path: Path to UKBB model VCF file
        output_dir: Directory to write VCF files
        bed_file_path: Path to BED file defining exome regions
        num_vcfs: Number of VCF files to generate
        num_rows: Number of variants to generate per chromosome (default: auto-calculated)
        seed: Base random seed for reproducibility
        num_threads: Number of parallel workers (default: CPU count)
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    if num_threads is None:
        num_threads = multiprocessing.cpu_count()
    
    # Calculate num_rows if not provided (need to load BED file)
    if num_rows is None:
        bed_intervals = BEDIntervals(bed_file_path)
        total_bases = bed_intervals.get_total_bases()
        num_rows = max(100, total_bases // 1000)
    
    print(f"Generating {num_vcfs:,} single-sample exome VCF files")
    print(f"Output directory: {output_path}")
    print(f"Using {num_threads} parallel workers")
    print(f"Base seed: {seed}")
    print(f"Variants per chromosome: {num_rows}")
    
    # Write batch metadata file
    batch_metadata = {
        "batch_id": str(uuid.uuid4()),
        "num_samples": num_vcfs,
        "model_vcf_path": str(model_vcf_path),
        "bed_file_path": str(bed_file_path),
        "num_rows_per_chromosome": num_rows,
        "base_seed": seed,
        "num_threads": num_threads,
        "output_directory": str(output_path),
        "type": "exome" if bed_file_path else "genome",
    }
    
    import json
    metadata_path = output_path / "batch_metadata.json"
    with open(metadata_path, 'w') as f:
        json.dump(batch_metadata, f, indent=2)
    print(f"Batch metadata written to {metadata_path}")
    
    # Prepare arguments for workers
    args_list = [
        (model_vcf_path, bed_file_path, output_path, num_rows, seed, i)
        for i in range(num_vcfs)
    ]
    
    # Generate VCFs in parallel
    success_count = 0
    failed_count = 0
    sample_sheet_data = []
    
    with multiprocessing.Pool(num_threads) as pool:
        if HAS_TQDM:
            results = list(tqdm(
                pool.imap(_generate_single_vcf, args_list),
                total=num_vcfs,
                desc="Generating VCFs",
                unit="file"
            ))
        else:
            # Fallback without tqdm
            results = []
            for i, result in enumerate(pool.imap(_generate_single_vcf, args_list)):
                results.append(result)
                if (i + 1) % 1000 == 0:
                    print(f"Generated {i + 1:,}/{num_vcfs:,} VCF files...")
    
    for result in results:
        success, uuid_val, variant_count, sex, variant_counts_by_chrom, is_exome, sample_seed, error = result
        if success:
            success_count += 1
            # Add to sample sheet
            sample_sheet_data.append({
                'sample_uuid': uuid_val,
                'vcf_file': f"{uuid_val}.vcf.gz",
                'num_variants': variant_count,
                'sex': sex,
                'type': 'exome' if is_exome else 'genome',
                'chrX_variants': variant_counts_by_chrom.get('chrX', 0),
                'chrY_variants': variant_counts_by_chrom.get('chrY', 0),
                'seed': sample_seed,
            })
        else:
            failed_count += 1
            print(f"Failed to generate {uuid_val}: {error}")
    
    # Write sample sheet
    sample_sheet_path = output_path / "sample_sheet.csv"
    if sample_sheet_data:
        with open(sample_sheet_path, 'w', newline='') as f:
            fieldnames = ['sample_uuid', 'vcf_file', 'num_variants', 'sex', 'type', 'chrX_variants', 'chrY_variants', 'seed']
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(sample_sheet_data)
        print(f"\n✓ Sample sheet written to {sample_sheet_path}")
    
    print(f"\n✓ Generated {success_count:,} VCF files successfully")
    if failed_count > 0:
        print(f"✗ Failed to generate {failed_count:,} VCF files")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate UKBB exome VCF files")
    parser.add_argument("--batch", action="store_true", help="Generate batch of VCF files")
    parser.add_argument("--num-vcfs", type=int, default=100000, help="Number of VCF files to generate (batch mode)")
    parser.add_argument("--output-dir", type=str, default="ukbb-exomes", help="Output directory for batch generation")
    parser.add_argument("--num-rows", type=int, default=None, help="Number of variants per chromosome")
    parser.add_argument("--num-threads", type=int, default=None, help="Number of parallel workers")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    
    args = parser.parse_args()
    
    model_path = "ukbb-model.vcf"
    bed_file_path = "S31285117_Covered.bed"
    
    if args.batch:
        # Batch generation mode
        generate_batch_ukbb_vcf(
            model_vcf_path=model_path,
            output_dir=args.output_dir,
            bed_file_path=bed_file_path,
            num_vcfs=args.num_vcfs,
            num_rows=args.num_rows,
            seed=args.seed,
            num_threads=args.num_threads,
        )
    else:
        # Single file generation mode
        output_path = "ukbb-synthetic-exome.vcf"
        variant_count, sample_uuid, sex, variant_counts_by_chrom, is_exome = generate_ukbb_vcf(
            model_vcf_path=model_path,
            output_path=output_path,
            bed_file_path=bed_file_path,
            num_rows=args.num_rows,
            seed=args.seed,
        )
        print(f"✓ Generated {variant_count} variants in {output_path}")
        print(f"  Sample UUID: {sample_uuid}")
        print(f"  Sex: {sex}")
        print(f"  Type: {'exome' if is_exome else 'genome'}")
        print(f"  chrX variants: {variant_counts_by_chrom.get('chrX', 0)}")
        print(f"  chrY variants: {variant_counts_by_chrom.get('chrY', 0)}")

