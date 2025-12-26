"""Cyrcular external tool wrapper.

Provides typed methods for cyrcular subcommands used by CircleSeeker.
"""

from __future__ import annotations

from pathlib import Path

from circleseeker.external.base import ExternalTool


class Cyrcular(ExternalTool):
    """Wrapper for the `cyrcular` CLI."""

    tool_name = "cyrcular"

    def graph_breakends(
        self,
        bam: Path,
        reference: Path,
        output_candidates: Path,
        output_graph: Path,
        dot_dir: Path,
        *,
        min_read_depth: int = 2,
        min_split_reads: int = 2,
        max_paths_per_component: int = 20,
        max_deletion_length: int = 1000,
        threads: int = 1,
    ) -> None:
        """Run breakend graph detection.

        Produces candidate BCF and graph outputs.
        """
        output_candidates.parent.mkdir(parents=True, exist_ok=True)
        dot_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            self.tool_name,
            "graph",
            "breakends",
            str(bam),
            "--reference",
            str(reference),
            "--min-read-depth",
            str(min_read_depth),
            "--min-split-reads",
            str(min_split_reads),
            "--max-paths-per-component",
            str(max_paths_per_component),
            "--max-deletion-length",
            str(max_deletion_length),
            "-t",
            str(threads),
            "--output",
            str(output_candidates),
            "--graph",
            str(output_graph),
            "--dot",
            f"{dot_dir}/",
        ]
        stdout, stderr = self.run(cmd, capture_output=True)
        if stderr:
            # cyrcular outputs progress to stderr, only log if it looks like an error
            if "WARNING" in stderr or "ERROR" in stderr:
                self.logger.warning(f"cyrcular graph breakends: {stderr[:1000]}")
            else:
                self.logger.debug(f"cyrcular graph breakends progress: {stderr[:500]}")
        self.logger.info(f"Breakend graph created: {output_candidates}")

    def graph_annotate(
        self,
        reference: Path,
        gene_annotation_gff_gz: Path,
        regulatory_annotation_gff_gz: Path,
        graph_input: Path,
        output_graph: Path,
    ) -> None:
        """Annotate graph with optional GFF3 data (can be empty placeholders)."""
        output_graph.parent.mkdir(parents=True, exist_ok=True)

        cmd = [
            self.tool_name,
            "graph",
            "annotate",
            "--reference",
            str(reference),
            "--gene-annotation",
            str(gene_annotation_gff_gz),
            "--regulatory-annotation",
            str(regulatory_annotation_gff_gz),
            "--output",
            str(output_graph),
            str(graph_input),
        ]
        stdout, stderr = self.run(cmd, capture_output=True)
        if stderr:
            if "WARNING" in stderr or "ERROR" in stderr:
                self.logger.warning(f"cyrcular graph annotate: {stderr[:1000]}")
            else:
                self.logger.debug(f"cyrcular graph annotate progress: {stderr[:500]}")
        self.logger.info(f"Annotated graph created: {output_graph}")

    def graph_table(
        self,
        annotated_graph: Path,
        calls_bcf: Path,
        reference: Path,
        circle_table: Path,
        segment_tables_dir: Path,
    ) -> None:
        """Generate overview and per-circle segment tables."""
        circle_table.parent.mkdir(parents=True, exist_ok=True)
        segment_tables_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            self.tool_name,
            "graph",
            "table",
            str(annotated_graph),
            str(calls_bcf),
            "--reference",
            str(reference),
            "--circle-table",
            str(circle_table),
            "--segment-tables",
            f"{segment_tables_dir}/",
        ]
        stdout, stderr = self.run(cmd, capture_output=True)
        if stderr:
            # cyrcular table often outputs WARNING messages that are informational
            if "WARNING:" in stderr:
                # Extract and log each warning separately at debug level
                for line in stderr.split("\n"):
                    if "WARNING:" in line:
                        self.logger.debug(f"cyrcular: {line.strip()}")
            elif "ERROR" in stderr:
                self.logger.warning(f"cyrcular graph table: {stderr[:1000]}")
            else:
                self.logger.debug(f"cyrcular graph table progress: {stderr[:500]}")
        self.logger.info(f"Circle table created: {circle_table}")
