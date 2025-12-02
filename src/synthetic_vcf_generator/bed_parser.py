"""BED file parser for filtering VCF generation to specific genomic regions."""

from pathlib import Path
from typing import Dict, List, Tuple


class BEDIntervals:
    """Parser and storage for BED file intervals."""

    def __init__(self, bed_file_path: str):
        """
        Parse a BED file and store intervals by chromosome.

        Args:
            bed_file_path: Path to BED file (tab-separated: chr, start, end, ...)
        """
        self.intervals: Dict[str, List[Tuple[int, int]]] = {}
        self._parse_bed_file(bed_file_path)

    def _parse_bed_file(self, bed_file_path: str) -> None:
        """
        Parse BED file and store intervals.

        BED format: chr start end [name] [score] [strand] ...
        BED uses 0-based half-open coordinates [start, end)
        """
        bed_path = Path(bed_file_path)
        if not bed_path.exists():
            raise FileNotFoundError(f"BED file not found: {bed_file_path}")

        with open(bed_path, "r") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                # Skip empty lines and comments
                if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                    continue

                fields = line.split("\t")
                if len(fields) < 3:
                    raise ValueError(
                        f"Invalid BED format at line {line_num}: expected at least 3 fields, got {len(fields)}"
                    )

                chrom = fields[0]
                try:
                    start = int(fields[1])
                    end = int(fields[2])
                except ValueError as e:
                    raise ValueError(f"Invalid coordinates at line {line_num}: {e}")

                if start < 0 or end < 0:
                    raise ValueError(f"Negative coordinates at line {line_num}: start={start}, end={end}")
                if start >= end:
                    raise ValueError(f"Invalid interval at line {line_num}: start={start} >= end={end}")

                # Store interval for this chromosome
                if chrom not in self.intervals:
                    self.intervals[chrom] = []
                self.intervals[chrom].append((start, end))

        # Sort intervals for each chromosome by start position
        for chrom in self.intervals:
            self.intervals[chrom].sort()

    def get_chromosomes(self) -> List[str]:
        """Return list of chromosomes in BED file."""
        return list(self.intervals.keys())

    def get_intervals(self, chromosome: str) -> List[Tuple[int, int]]:
        """
        Get intervals for a specific chromosome.

        Args:
            chromosome: Chromosome name

        Returns:
            List of (start, end) tuples in BED format (0-based half-open)
        """
        return self.intervals.get(chromosome, [])

    def contains_position(self, chromosome: str, position: int) -> bool:
        """
        Check if a position (1-based VCF coordinate) is within any interval.

        Args:
            chromosome: Chromosome name
            position: 1-based position (VCF format)

        Returns:
            True if position is within any BED interval for this chromosome
        """
        if chromosome not in self.intervals:
            return False

        # Convert VCF 1-based position to BED 0-based coordinate
        bed_pos = position - 1

        # Binary search for efficiency
        intervals = self.intervals[chromosome]
        for start, end in intervals:
            if start <= bed_pos < end:
                return True
            if start > bed_pos:
                # Since sorted, no need to check further
                break

        return False

    def get_valid_positions(self, chromosome: str, positions: List[int]) -> List[int]:
        """
        Filter a list of positions to only those within BED intervals.

        Args:
            chromosome: Chromosome name
            positions: List of 1-based VCF positions

        Returns:
            Filtered list of positions that fall within BED intervals
        """
        if chromosome not in self.intervals:
            return []

        return [pos for pos in positions if self.contains_position(chromosome, pos)]

    def get_total_bases(self, chromosome: str = None) -> int:
        """
        Calculate total bases covered by intervals.

        Args:
            chromosome: If specified, only count bases for this chromosome.
                       If None, count all chromosomes.

        Returns:
            Total number of bases covered
        """
        if chromosome:
            intervals = self.intervals.get(chromosome, [])
            return sum(end - start for start, end in intervals)
        else:
            return sum(
                sum(end - start for start, end in intervals)
                for intervals in self.intervals.values()
            )
