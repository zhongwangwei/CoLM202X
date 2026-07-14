#!/usr/bin/env python3
"""Static regression checks for tracer collective and gather costs."""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]


def routine(source: str, name: str) -> str:
    return source.split(f"SUBROUTINE {name}", 1)[1].split(
        f"END SUBROUTINE {name}", 1
    )[0]


class TracerMPIStaticChecks(unittest.TestCase):
    def test_balance_normal_path_has_three_collectives(self) -> None:
        source = (
            ROOT / "main/TRACER/MOD_Tracer_Conservation.F90"
        ).read_text(encoding="utf-8")
        report = routine(source, "tracer_balance_report")
        unconditional = report.split(
            "IF (resid_hard_nbad_total > 0) THEN", 1
        )[0]
        self.assertEqual(
            len(re.findall(r"CALL\s+mpi_(?:reduce|bcast)", unconditional, re.I)),
            3,
        )
        self.assertIn("maxima_local, maxima_total, 3", unconditional)
        self.assertIn("counts_local, counts_total, 3", unconditional)
        self.assertIn("mpi_bcast(counts_total, 3", unconditional)

    def test_maxloc_collectives_run_only_for_failures(self) -> None:
        source = (
            ROOT / "main/TRACER/MOD_Tracer_Conservation.F90"
        ).read_text(encoding="utf-8")
        report = routine(source, "tracer_balance_report")
        self.assertRegex(
            report,
            r"IF \(resid_hard_nbad_total > 0\) THEN[\s\S]*?MPI_MAXLOC",
        )
        self.assertRegex(report, r"IF \(nbad_total > 0\) THEN[\s\S]*?MPI_MAXLOC")
        self.assertEqual(
            len(
                re.findall(
                    r"CALL\s+mpi_allreduce[\s\S]*?MPI_MAXLOC",
                    report,
                    re.I,
                )
            ),
            2,
        )

    def test_gathers_keep_only_the_field_boundary_barrier(self) -> None:
        source = (ROOT / "main/HYDRO/MOD_Vector_ReadWrite.F90").read_text(
            encoding="utf-8"
        )
        for name in (
            "vector_gather_to_master",
            "vector_gather_matrix_to_master",
        ):
            gather = routine(source, name)
            self.assertEqual(
                len(re.findall(r"CALL\s+mpi_barrier", gather, re.I)), 1, name
            )
            self.assertGreater(gather.index("CALL mpi_barrier"), gather.index("CALL mpi_recv"))

    def test_custom_history_gathers_keep_trailing_field_barrier(self) -> None:
        source = (ROOT / "main/TRACER/MOD_Tracer_Hist.F90").read_text(
            encoding="utf-8"
        )
        for name in (
            "write_history_tracer_vector_2d",
            "write_history_tracer_ratio_vector_3d",
        ):
            gather = routine(source, name)
            self.assertEqual(
                len(re.findall(r"CALL\s+mpi_barrier", gather, re.I)), 1, name
            )
            self.assertGreater(
                gather.index("CALL mpi_barrier"), gather.index("CALL mpi_recv")
            )


if __name__ == "__main__":
    unittest.main()
