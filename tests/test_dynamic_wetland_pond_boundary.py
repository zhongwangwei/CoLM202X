"""Ponded wetlands and dry lakes must not lose sub-tolerance surface fluxes."""

from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


SOURCE = Path(__file__).resolve().parents[1] / "main/MOD_SoilSnowHydrology.F90"


def test_pond_boundary_preserves_small_evaporation_and_rain(tmp_path):
    source = SOURCE.read_text()
    start = source.index("      ! The Richards rainfall boundary can discard a flux below its")
    end = source.index("\n#if (defined CaMa_Flood) || (defined GridRiverLakeFlow)", start)
    boundary = source[start:end]
    assert "richards_water_tolerance, wblc)" in source
    program = f"""
program check_wetland_pond
   implicit none
   integer, parameter :: r8 = kind(1.d0)
   real(r8), parameter :: richards_water_tolerance = 1.e-3_r8
   real(r8) :: qgtop, wdsrf, deltim, original, storage_before
   integer :: patchtype, case_id
   logical :: DEF_USE_Dynamic_Wetland, is_dry_lake
   deltim = 1800._r8
   DEF_USE_Dynamic_Wetland = .true.
   do case_id = 1, 10
      patchtype = 2
      is_dry_lake = .false.
      wdsrf = 52.8_r8
      select case (case_id)
      case (1) ! Reproduced negative VSF boundary: 9.16e-5 mm vanished.
         qgtop = -9.1581171e-5_r8/deltim
      case (2) ! Reproduced positive boundary: 6.26e-6 mm vanished.
         qgtop = 6.2631257e-6_r8/deltim
      case (3) ! Resolvable rainfall still enters Richards.
         qgtop = 0.01_r8/deltim
      case (4) ! Resolved evaporation retains the Richards/root-ET order.
         wdsrf = 1._r8
         qgtop = -2._r8/deltim
      case (5) ! Non-wetland remains unchanged.
         patchtype = 0
         qgtop = -9.1581171e-5_r8/deltim
      case (6) ! Dry wetland remains unchanged.
         wdsrf = 0._r8
         qgtop = -9.1581171e-5_r8/deltim
      case (7) ! A tiny demand exceeding the pond still uses Richards.
         wdsrf = 5.e-5_r8
         qgtop = -9.1581171e-5_r8/deltim
      case (8) ! Reproduced dry-lake negative VSF boundary: 1.10e-4 mm vanished.
         patchtype = 4
         is_dry_lake = .true.
         wdsrf = 96.337_r8
         qgtop = -1.098072136e-4_r8/deltim
      case (9) ! Dry-lake small positive boundary is also conserved.
         patchtype = 4
         is_dry_lake = .true.
         qgtop = 6.2631257e-6_r8/deltim
      case (10) ! A filled lake does not take the dry-lake soil path.
         patchtype = 4
         qgtop = -1.098072136e-4_r8/deltim
      end select
      original = qgtop
      storage_before = wdsrf
{boundary}
      if (abs((wdsrf-storage_before) + (qgtop-original)*deltim) > 1.e-12_r8) stop 1
      select case (case_id)
      case (1, 2, 8, 9)
         if (abs(qgtop) > 1.e-20_r8) stop 2
      case (3, 4, 5, 6, 7, 10)
         if (abs(qgtop-original) > 1.e-20_r8) stop 3
      end select
   end do
end program
"""
    compiler = require_runnable_fortran_compiler(tmp_path)
    source_path = tmp_path / "check_wetland_pond.f90"
    exe = tmp_path / "check_wetland_pond"
    source_path.write_text(program)
    built = subprocess.run(
        [compiler, "-fcheck=all", str(source_path), "-o", str(exe)],
        capture_output=True,
        text=True,
    )
    assert built.returncode == 0, built.stderr
    ran = subprocess.run([str(exe)], capture_output=True, text=True)
    assert ran.returncode == 0, ran.stderr + ran.stdout
