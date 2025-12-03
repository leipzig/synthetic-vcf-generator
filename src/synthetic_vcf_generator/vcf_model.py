"""VCF model parser for extracting FORMAT and INFO field definitions."""

import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional


class VCFModel:
    """Parser for VCF model files to extract field definitions and generate realistic values."""

    def __init__(self, model_vcf_path: str):
        """
        Parse a model VCF file and extract field definitions.

        Args:
            model_vcf_path: Path to model VCF file
        """
        self.model_path = Path(model_vcf_path)
        self.format_fields: Dict[str, Dict] = {}
        self.info_fields: Dict[str, Dict] = {}
        self.header_lines: List[str] = []
        self._parse_model()

    def _parse_model(self):
        """Parse model VCF to extract field definitions."""
        if not self.model_path.exists():
            raise FileNotFoundError(f"Model VCF not found: {self.model_path}")

        with open(self.model_path, 'r') as f:
            for line in f:
                if not line.startswith('##'):
                    break
                
                self.header_lines.append(line.rstrip())
                
                # Parse FORMAT fields
                if line.startswith('##FORMAT='):
                    field_info = self._parse_vcf_meta_line(line)
                    if field_info:
                        self.format_fields[field_info['ID']] = field_info
                
                # Parse INFO fields
                elif line.startswith('##INFO='):
                    field_info = self._parse_vcf_meta_line(line)
                    if field_info:
                        self.info_fields[field_info['ID']] = field_info

    def _parse_vcf_meta_line(self, line: str) -> Optional[Dict]:
        """
        Parse a VCF metadata line like ##FORMAT=<...> or ##INFO=<...>

        Args:
            line: Metadata line from VCF header

        Returns:
            Dictionary with field properties or None if parsing fails
        """
        # Extract content between < >
        match = re.search(r'<(.+)>', line)
        if not match:
            return None

        content = match.group(1)
        field_info = {}

        # Parse key=value pairs, handling quoted values
        pattern = r'(\w+)=("(?:[^"\\]|\\.)*"|[^,]+)'
        for match in re.finditer(pattern, content):
            key = match.group(1)
            value = match.group(2)
            # Remove quotes if present
            if value.startswith('"') and value.endswith('"'):
                value = value[1:-1]
            field_info[key] = value

        return field_info

    def get_format_string(self) -> str:
        """
        Get FORMAT field string for VCF header.

        Returns:
            Colon-separated FORMAT field names
        """
        # Maintain standard order for common fields
        priority_fields = ['GT', 'AD', 'DP', 'AF', 'GQ', 'PL', 'GP']
        
        # Get all field IDs
        all_fields = list(self.format_fields.keys())
        
        # Sort: priority fields first, then the rest alphabetically
        format_list = []
        for field in priority_fields:
            if field in all_fields:
                format_list.append(field)
                all_fields.remove(field)
        
        format_list.extend(sorted(all_fields))
        return ':'.join(format_list)

    def get_header_lines(self) -> List[str]:
        """
        Get all header lines from model VCF.

        Returns:
            List of header lines (without fileformat and source which will be regenerated)
        """
        # Filter out lines we'll regenerate
        filtered = []
        for line in self.header_lines:
            if not any(line.startswith(prefix) for prefix in ['##fileformat', '##source']):
                filtered.append(line)
        return filtered

    def has_field(self, field_type: str, field_id: str) -> bool:
        """
        Check if model has a specific field.

        Args:
            field_type: 'FORMAT' or 'INFO'
            field_id: Field identifier (e.g., 'GT', 'DP')

        Returns:
            True if field exists in model
        """
        if field_type == 'FORMAT':
            return field_id in self.format_fields
        elif field_type == 'INFO':
            return field_id in self.info_fields
        return False

    def get_field_type(self, field_type: str, field_id: str) -> Optional[str]:
        """
        Get the Type of a field (String, Integer, Float, Flag).

        Args:
            field_type: 'FORMAT' or 'INFO'
            field_id: Field identifier

        Returns:
            Type string or None if not found
        """
        fields = self.format_fields if field_type == 'FORMAT' else self.info_fields
        if field_id in fields:
            return fields[field_id].get('Type')
        return None
