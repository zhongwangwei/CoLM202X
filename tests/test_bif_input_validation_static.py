from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


ROOT = Path(__file__).resolve().parents[1]
NETWORK = (ROOT / "main/HYDRO/MOD_Grid_RiverLakeNetwork.F90").read_text()


def test_bif_parameter_shapes_and_finite_values_are_validated() -> None:
    assert "bifurcation parameter dimensions are inconsistent" in NETWORK
    assert "ieee_is_finite(bif_dist_all)" in NETWORK
    assert "ieee_is_finite(bif_elev_all)" in NETWORK
    assert "ieee_is_finite(bif_wdth_all)" in NETWORK
    assert "ieee_is_finite(bif_mann_all)" in NETWORK


def test_bif_active_geometry_requires_positive_distance_width_and_manning() -> None:
    assert "any(bif_dist_all <= 0._r8)" in NETWORK
    finite_distance = NETWORK.index("any(.not. ieee_is_finite(bif_dist_all))")
    positive_distance = NETWORK.index("any(bif_dist_all <= 0._r8)")
    assert finite_distance < positive_distance
    assert "ieee_is_finite(bif_dist_all)) .or." not in NETWORK
    assert "active bifurcation layer requires a positive Manning coefficient" in NETWORK
    assert "bifurcation pathway has no active positive-width layer" in NETWORK


def test_sparse_bifurcation_layers_runtime(tmp_path: Path) -> None:
    """Compile the production post-read validator, without NetCDF/MPI setup."""
    compiler = require_runnable_fortran_compiler(tmp_path)
    reader = NETWORK.split("   SUBROUTINE read_bifurcation_global_arrays (", 1)[1]
    reader = reader.split("   END SUBROUTINE read_bifurcation_global_arrays", 1)[0]
    declarations = reader.split("      CALL ncio_inquire_length", 1)[0]
    declarations = declarations[declarations.index("   integer ::"):]
    validation = reader[reader.index("      IF (size(bif_down_all)"):]
    source = tmp_path / "validate_bif.f90"
    source.write_text("""
program validate_bif
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   integer, parameter :: r8 = kind(1.d0), npthlev_bif = 5, totalnpthout = 2, totalnumucat = 2
   integer :: bif_upst_all(2) = [1, 2], bif_down_all(2) = [2, 1]
   real(r8) :: bif_dist_all(2) = 100._r8
   real(r8) :: bif_wdth_all(5,2), bif_elev_all(5,2), bif_mann_all(5)
   real(r8) :: original_width(5,2), original_elev(5,2)
   read(*,*) bif_wdth_all(:,1)
   read(*,*) bif_elev_all(:,1)
   read(*,*) bif_mann_all
   ! The second path also checks that the previous active layer is reset.
   bif_wdth_all(:,2) = [1._r8, 0._r8, 0._r8, 0._r8, 0._r8]
   bif_elev_all(:,2) = [-100._r8, 1.e20_r8, 1.e20_r8, 1.e20_r8, 1.e20_r8]
   original_width = bif_wdth_all
   original_elev = bif_elev_all
   call validate()
   if (any(bif_wdth_all /= original_width)) error stop 'widths were altered'
   if (any(bif_elev_all /= original_elev)) error stop 'elevations were altered'
   print '(A)', 'PASS'
contains
   subroutine validate()
""" + declarations + validation + """
   end subroutine validate
   subroutine CoLM_stop(message)
      character(len=*), intent(in) :: message
      print '(A)', message
      stop 1
   end subroutine CoLM_stop
end program validate_bif
""", encoding="utf-8")
    executable = tmp_path / "validate_bif"
    compiled = subprocess.run(
        [compiler, "-fcheck=all", "-ffpe-trap=invalid,zero,overflow",
         "-ffree-line-length-0", str(source), "-o", str(executable)],
        cwd=tmp_path, capture_output=True, text=True, timeout=60,
    )
    assert compiled.returncode == 0, compiled.stdout + compiled.stderr

    def check(widths, elevations, manning, error=""):
        result = subprocess.run(
            [str(executable)],
            input="\n".join(" ".join(map(str, values)) for values in (widths, elevations, manning)) + "\n",
            capture_output=True, text=True, timeout=10,
        )
        assert result.returncode == (1 if error else 0), result.stdout + result.stderr
        assert (error or "PASS") in result.stdout

    profiles = [[0, 0, 0, 46.97, 93.95], [0, 46.64, 140.89, 0, 94.73]]
    profiles += [[(mask >> level) & 1 for level in range(5)] for mask in range(1, 32)]
    for widths in profiles:
        elevations = [float(level) if width else 1.e20 for level, width in enumerate(widths)]
        check(widths, elevations, [0.03] * 5)
    check([0, 1, 0, 1, 0], [1.e20, 2, 1.e20, 2, 1.e20], [0.03] * 5)
    check([0, 1, 0, 1, 0], [1.e20, 2, 1.e20, 1, 1.e20], [0.03] * 5,
          "active bifurcation layer elevation must be non-decreasing")
    check([0] * 5, [1.e20] * 5, [0.03] * 5, "bifurcation pathway has no active positive-width layer")
    check([0, -1, 0, 1, 0], list(range(5)), [0.03] * 5, "bifurcation width must be non-negative")
    check([0, 0, 0, 1, 1], list(range(5)), [0.03, 0, 0, 0.03, 0.03])
    check([0, 0, 0, 1, 1], list(range(5)), [0.03, 0.03, 0.03, 0, 0.03],
          "active bifurcation layer requires a positive Manning coefficient")
