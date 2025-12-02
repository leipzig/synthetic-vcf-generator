from __future__ import annotations

import json
import random
import uuid
from collections import deque
from pathlib import Path
from typing import List, Literal

import fastrand

from synthetic_vcf_generator import vcf_reference, version
from synthetic_vcf_generator.bed_parser import BEDIntervals


class VirtualVCF:
    def __init__(
        self,
        num_rows: int,
        num_samples: int,
        chromosomes: List[str],
        sample_prefix: str | None = "sample_",
        id_type: Literal["count", "padded_count", "uuid"] | None = "padded_count",
        random_seed: int | None = None,
        phased: bool | None = True,
        large_format: bool | None = True,
        reference_dir: str | Path | None = None,
        bed_intervals: BEDIntervals | None = None,
    ):
        """
        Initialize VirtualVCF object.

        Args:
            num_rows (int): Number of rows (variants) to generate per chromosome.
            num_samples (int): Number of samples.
            chromosomes (List[str]): List of chromosome IDs to include in the VCF.
            sample_prefix (str, optional): Prefix for sample names. Defaults to "sample_".
            id_type (str): Type of unique ID to use for samples.
            random_seed (int, optional): Random seed for reproducibility. Defaults to None.
            phased (bool, optional): Phased or unphased genotypes. Defaults to True.
            large_format (bool, optional): Use large format VCF. Defaults to True.
            reference_dir (str or Path, optional): Path to reference file directory.
            bed_intervals (BEDIntervals, optional): BED intervals to restrict variant generation.
            bed_intervals (BEDIntervals, optional): BED intervals to restrict variant generation.

        Raises:
            ValueError: If num_samples or num_rows is less than 1.
        """
        if num_samples < 1:
            raise ValueError(
                f"Number of samples must be greater than 0: {num_samples=}"
            )
        if num_rows < 1:
            raise ValueError(f"Number of rows must be greater than 0: {num_rows=}")
        if id_type not in {"count", "padded_count", "uuid"}:
            raise ValueError(f"Unexpected sample ID type: {id_type=}")

        self.num_rows = num_rows
        self.num_samples = num_samples
        self.max_rotation = num_samples // 10 if num_samples >= 10 else num_samples
        self.chromosomes = chromosomes
        self.sample_prefix = sample_prefix
        self.id_type = id_type
        self.phased = phased
        self.random = random.Random(random_seed)  # seed for reproducibility
        if random_seed is not None:
            fastrand.pcg32_seed(random_seed)
        self.large_format = large_format
        self.fmt = "GT:AD:DP:GQ:PL" if large_format else "GT"
        self.info = f"DP=10;AF=0.5;NS={num_samples}"
        self.reference_dir = Path(reference_dir) if reference_dir else None
        self.reference_files = {}
        self.reference_metadata = {}
        self.bed_intervals = bed_intervals
        self.alleles = ["A", "C", "G", "T"]

        # Setup
        self._setup_reference_data()
        self._generate_samples()

    def _generate_vcf_header(self):
        """
        Generates the VCF header.
        """
        # Lines included in every header
        # TODO: include length, etc when using reference_dir
        contigs = "\n".join([f"##contig=<ID={c}>" for c in self.chromosomes])
        header_lines = [
            "##fileformat=VCFv4.2",
            f"##source=VirtualVCF {version}",
            '##FILTER=<ID=PASS,Description="All filters passed">',
            '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">',
            contigs,
            f"##reference=ftp://ftp.example.com/{self.reference_metadata.get('source_reference_file', 'sample.fa')}",
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">',
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        ]

        # Add large format lines
        if self.large_format:
            header_lines += [
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">',
                '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
                '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype Likelihoods">',
            ]

        # Add column line
        columns = [
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
        ]
        for i in range(1, self.num_samples + 1):
            if self.id_type == "count":
                sample_name = f"{self.sample_prefix}{i}"
            elif self.id_type == "padded_count":
                sample_name = f"{self.sample_prefix}{i:07d}"
            elif self.id_type == "uuid":
                sample_name = f"{self.sample_prefix}{uuid.uuid4()}"
            columns.append(sample_name)
        header_lines += ["\t".join(columns)]

        # Join the lines into a single header
        header = "\n".join(header_lines) + "\n"
        return header

    def _get_ref_alt_at_position(self, position, reference_data):
        """
        Retrieves the reference value at a given position if it exists in reference data
        or returns the allele at the given index.
        """
        ref_index = fastrand.pcg32randint(0, 3)
        allele_index = fastrand.pcg32randint(1, 3)
        if reference_data:
            ref = reference_data.get_ref_at_pos(position - 1)
        else:
            ref = self.alleles[ref_index]
        if ref in self.alleles:
            alt = self.alleles[self.alleles.index(ref) - allele_index]
        else:
            alt = self.alleles[ref_index - allele_index]
        return ref, alt

    def _generate_vcf_row(self, chromosome, position, reference_data):
        """
        Generates a VCF row.
        """
        # Generate random values for each field in the VCF row
        # TODO: make ID strategy configurable via CLI
        vid = "."
        ref, alt = self._get_ref_alt_at_position(position, reference_data)
        qual = f"{fastrand.pcg32randint(10, 100)}"
        # TODO: add support for other filters and make configurable via CLI
        # filter = self.random.choice(["PASS"])
        filter = "PASS"
        # TODO: copmute info rather than using static values

        # Random select a sample rotation
        self.available_samples.rotate(fastrand.pcg32randint(1, self.max_rotation))
        samples = "\t".join(self.available_samples.copy())

        # Make the row
        row = f"{chromosome}\t{position}\t{vid}\t{ref}\t{alt}\t{qual}\t{filter}\t{self.info}\t{self.fmt}\t{samples}\n"

        return row

    def _setup_reference_data(self):
        """
        Load refrence metadata and maps chromosomes to reference files.
        """
        if self.reference_dir:
            with open(
                self.reference_dir / vcf_reference.METADATA_FILE_NAME
            ) as metadata_file:
                self.reference_metadata = json.load(metadata_file)

            for chromosome in self.chromosomes:
                if chromosome not in self.reference_metadata["reference_files"].keys():
                    raise ValueError(
                        f'"{chromosome}" does not exist in the reference data'
                    )
                chromosome_file = (
                    self.reference_dir
                    / self.reference_metadata["reference_files"][chromosome]
                )
                self.reference_files[chromosome] = chromosome_file

    def _generate_samples(self):
        """
        Generates sample format fields to choose from when generating data.
        """
        if self.phased:
            self.sample_values = ["0|0", "1|0", "0|1", "1|1"]
            self.sample_value_weights = [
                self.num_samples * 10,
                int(self.num_samples / 500),
                int(self.num_samples / 500),
                int(self.num_samples / 300),
            ]
        else:
            self.sample_values = ["0/0", "0/1", "1/1"]
            self.sample_value_weights = [
                self.num_samples * 10,
                int(self.num_samples / 250),
                int(self.num_samples / 300),
            ]

        if self.large_format:
            extra_data = [
                "0,30:30:89:913,89,0",
                "0,10:10:49:413,33,0",
                "0,20:20:55:489,89,0",
                "0,40:00:66:726,85,0",
            ]
            self.sample_values = [
                f"{sv}:{self.random.choice(extra_data)}" for sv in self.sample_values
            ]

        self.available_samples = deque(
            self.random.choices(
                self.sample_values,
                weights=self.sample_value_weights,
                k=self.num_samples,
            )
        )
        # Generate sample rotations for later random selection
        # self.sample_rotations = []
        # rotation = self.available_samples.copy()
        # for i in range(self.max_rotation):
        #     rotation.rotate(1)
        #     self.sample_rotations.append("\t".join(rotation))

        # Check so that at lest one sample in avail_samples is not 0|0 or 0/0
        if all(
            [
                sample.startswith(self.sample_values[0])
                for sample in self.available_samples
            ]
        ):
            if self.phased:
                self.available_samples[0] = self.available_samples[0].replace(
                    "0", "1", 1
                )
            else:
                self.available_samples[0] = self.available_samples[0].replace(
                    "0/0", "0/1", 1
                )

    def __enter__(self):
        """
        Enters the context.
        """
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """
        Exits the context.
        """
        pass

    def __iter__(self):
        """
        Iterates over VirtualVCF object.
        """
        # Start with the header
        yield self._generate_vcf_header()
        # Generate VCF rows for each chromosome
        for chromosome in self.chromosomes:
            # TODO: this is a magic number
            chromosome_length = self.num_rows * 100
            # Load the chromosome's reference file
            reference_data = None
            if self.reference_dir:
                chromosome_file = self.reference_files[chromosome]
                reference_data = vcf_reference.load_reference_data(chromosome_file)
                reference_data.open()
                chromosome_length = reference_data.ref_length()
            # Generate and sort positions
            positions = self.random.sample(range(1, chromosome_length), self.num_rows)
            positions.sort()
            # Filter positions to BED intervals if provided
            if self.bed_intervals:
                positions = self.bed_intervals.get_valid_positions(chromosome, positions)
                # If no valid positions after filtering, skip this chromosome
                if not positions:
                    if reference_data:
                        reference_data.close()
                    continue
            # Generate a VCF row for each position
            for p in positions:
                yield self._generate_vcf_row(chromosome, p, reference_data)
            if reference_data:
                reference_data.close()
