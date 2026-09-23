from pathlib import Path
import subprocess

import pytest

from fortran_test_support import SMOKE_TIMEOUT, require_runnable_fortran_compiler


ROOT = Path(__file__).resolve().parents[1]
LEAF_SOURCES = (
    ROOT / "main/MOD_LeafTemperature.F90",
    ROOT / "extends/interception/MOD_LeafTemperature_Extended.F90",
)
PC_SOURCES = (
    (ROOT / "main/MOD_LeafTemperaturePC.F90", "hvap", "ldew(i)"),
    (
        ROOT / "extends/interception/MOD_LeafTemperaturePC_Extended.F90",
        "htvpl(i)",
        "ldew_VIC_evap(i)",
    ),
)


def test_scalar_negative_transpiration_is_routed_before_capacity_clip():
    for source, latent_heat in zip(LEAF_SOURCES, ("hvap", "htvpl")):
        text = source.read_text(encoding="utf-8")
        final_flux = text.split(
            "evplwet = evplwet + evplwet_dtl*dtl(it-1)", 1
        )[1].split("taux", 1)[0]
        assert "IF (etr < 0._r8) THEN" in final_flux
        assert final_flux.index("evplwet = evplwet + etr") < final_flux.index("etr = 0._r8")
        assert "etrsun = 0._r8" in final_flux
        assert "etrsha = 0._r8" in final_flux
        assert "IF (DEF_USE_PLANTHYDRAULICS) rootflux = 0._r8" in final_flux
        assert final_flux.index("etr = 0._r8") < final_flux.index("elwdif  = max(0., evplwet-elwmax)")
        assert f"fsenl = fsenl + {latent_heat}*elwdif" in final_flux


def test_extended_rootflux_balance_precedes_final_condensation_reclassification():
    text = LEAF_SOURCES[1].read_text(encoding="utf-8")
    final_update = text.split("etr   = etr + etr_dtl*dtl(it-1)", 1)[1].split(
        "elwmax  = ldew_VIC_evap/deltim", 1
    )[0]
    balance = "CALL balance_phs_rootflux"
    reclassify = "IF (etr < 0._r8) THEN"
    assert final_update.index(balance) < final_update.index(reclassify)
    assert balance not in final_update.split(reclassify, 1)[1]


@pytest.mark.parametrize("source", LEAF_SOURCES)
def test_scalar_reclassification_production_block(tmp_path, source):
    compiler = require_runnable_fortran_compiler(tmp_path)
    text = source.read_text(encoding="utf-8")
    marker = "! Dry-leaf transpiration cannot supply condensation to the roots."
    start = text.index(marker)
    end = text.index("      ENDIF", start) + len("      ENDIF")
    block = text[start:end]
    driver = tmp_path / "scalar_negative_etr_probe.f90"
    executable = tmp_path / "scalar_negative_etr_probe"
    driver.write_text(
        f"""program scalar_negative_etr_probe
  implicit none
  integer, parameter :: r8 = kind(1.d0)
  logical :: DEF_USE_PLANTHYDRAULICS
  real(r8) :: etr, evplwet, etrsun, etrsha, rootflux(2)

  DEF_USE_PLANTHYDRAULICS = .false.
  etr = -2._r8
  evplwet = 12._r8
  etrsun = 1._r8
  etrsha = 1._r8
  rootflux = [1._r8, 2._r8]
{block}
  if (etr /= 0._r8 .or. evplwet /= 10._r8) error stop 1
  if (etrsun /= 0._r8 .or. etrsha /= 0._r8) error stop 2
  if (any(rootflux /= [1._r8, 2._r8])) error stop 3

  DEF_USE_PLANTHYDRAULICS = .true.
  etr = -2._r8
  evplwet = 12._r8
  etrsun = 1._r8
  etrsha = 1._r8
  rootflux = [1._r8, 2._r8]
{block}
  if (etr /= 0._r8 .or. evplwet /= 10._r8) error stop 4
  if (etrsun /= 0._r8 .or. etrsha /= 0._r8) error stop 5
  if (any(rootflux /= 0._r8)) error stop 6

  etr = 3._r8
  evplwet = 12._r8
  etrsun = 1._r8
  etrsha = 1._r8
  rootflux = [5._r8, 6._r8]
{block}
  if (etr /= 3._r8 .or. evplwet /= 12._r8) error stop 7
  if (etrsun /= 1._r8 .or. etrsha /= 1._r8) error stop 8
  if (any(rootflux /= [5._r8, 6._r8])) error stop 9
  print '(A)', 'SCALAR_NEGATIVE_ETR_OK'
end program scalar_negative_etr_probe
""",
        encoding="utf-8",
    )
    built = subprocess.run(
        [compiler, "-ffree-line-length-none", str(driver), "-o", str(executable)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=SMOKE_TIMEOUT,
    )
    assert built.returncode == 0, built.stdout + built.stderr
    ran = subprocess.run(
        [str(executable)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=SMOKE_TIMEOUT,
    )
    assert ran.returncode == 0, ran.stdout + ran.stderr
    assert ran.stdout.strip() == "SCALAR_NEGATIVE_ETR_OK"


def test_pc_negative_transpiration_is_routed_before_capacity_clip():
    for source, latent_heat, canopy_storage in PC_SOURCES:
        text = source.read_text(encoding="utf-8")
        final_flux = text.split(
            "evplwet(i) = evplwet(i) + evplwet_dtl(i)*dtl(it-1,i)", 1
        )[1].split("! precipitation sensible heat", 1)[0]
        condition = "IF (etr(i) < 0._r8) THEN"
        assert condition in final_flux
        assert final_flux.index("evplwet(i) = evplwet(i) + etr(i)") < final_flux.index(
            "etr(i) = 0._r8"
        )
        assert "etrsun(i) = 0._r8" in final_flux
        assert "etrsha(i) = 0._r8" in final_flux
        assert "IF (DEF_USE_PLANTHYDRAULICS) rootflux(:,i) = 0._r8" in final_flux
        assert final_flux.index("etr(i) = 0._r8") < final_flux.index(
            f"elwmax = {canopy_storage}/deltim"
        )
        assert f"fsenl(i) = fsenl(i) + {latent_heat}*elwdif" in final_flux


@pytest.mark.parametrize("source,unused_latent_heat,unused_storage", PC_SOURCES)
def test_pc_reclassification_block_runs_per_pft(
    tmp_path, source, unused_latent_heat, unused_storage
):
    compiler = require_runnable_fortran_compiler(tmp_path)
    text = source.read_text(encoding="utf-8")
    marker = "! Dry-leaf transpiration cannot supply condensation to the roots."
    start = text.index(marker)
    end = text.index("            ENDIF", start) + len("            ENDIF")
    block = text[start:end]
    driver = tmp_path / "pc_negative_etr_probe.f90"
    executable = tmp_path / "pc_negative_etr_probe"
    driver.write_text(
        f"""program pc_negative_etr_probe
  implicit none
  integer, parameter :: r8 = kind(1.d0)
  integer :: i
  logical :: DEF_USE_PLANTHYDRAULICS
  real(r8) :: etr(2), evplwet(2), etrsun(2), etrsha(2), rootflux(2,2)

  DEF_USE_PLANTHYDRAULICS = .false.
  etr = [-2._r8, 3._r8]
  evplwet = [12._r8, 12._r8]
  etrsun = [1._r8, 1._r8]
  etrsha = [1._r8, 1._r8]
  rootflux = reshape([1._r8, 2._r8, 3._r8, 4._r8], shape(rootflux))
  do i = 1, 2
{block}
  end do
  if (etr(1) /= 0._r8 .or. evplwet(1) /= 10._r8) error stop 1
  if (etrsun(1) /= 0._r8 .or. etrsha(1) /= 0._r8) error stop 2
  if (etr(2) /= 3._r8 .or. evplwet(2) /= 12._r8) error stop 3
  if (etrsun(2) /= 1._r8 .or. etrsha(2) /= 1._r8) error stop 4

  DEF_USE_PLANTHYDRAULICS = .true.
  i = 1
  etr(1) = -2._r8
  evplwet(1) = 12._r8
  etrsun(1) = 1._r8
  etrsha(1) = 1._r8
{block}
  if (etr(1) /= 0._r8 .or. evplwet(1) /= 10._r8) error stop 5
  if (etrsun(1) /= 0._r8 .or. etrsha(1) /= 0._r8) error stop 6
  if (any(rootflux(:,1) /= 0._r8)) error stop 7

  i = 2
  etr(2) = 3._r8
  evplwet(2) = 12._r8
  etrsun(2) = 1._r8
  etrsha(2) = 1._r8
  rootflux(:,2) = [5._r8, 6._r8]
{block}
  if (etr(2) /= 3._r8 .or. evplwet(2) /= 12._r8) error stop 8
  if (etrsun(2) /= 1._r8 .or. etrsha(2) /= 1._r8) error stop 9
  if (any(rootflux(:,2) /= [5._r8, 6._r8])) error stop 10
  print '(A)', 'PC_NEGATIVE_ETR_OK'
end program pc_negative_etr_probe
""",
        encoding="utf-8",
    )
    built = subprocess.run(
        [compiler, "-ffree-line-length-none", str(driver), "-o", str(executable)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=SMOKE_TIMEOUT,
    )
    assert built.returncode == 0, built.stdout + built.stderr
    ran = subprocess.run(
        [str(executable)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=SMOKE_TIMEOUT,
    )
    assert ran.returncode == 0, ran.stdout + ran.stderr
    assert ran.stdout.strip() == "PC_NEGATIVE_ETR_OK"


def test_ordinary_pc_phase_storage_sync_runs_production_fragment(tmp_path):
    compiler = require_runnable_fortran_compiler(tmp_path)
    source = (ROOT / "main/MOD_LeafTemperaturePC.F90").read_text(encoding="utf-8")
    marker = "! Keep phase-resolved storage synchronized for downstream users"
    start = source.index(marker)
    end = source.index("            IF ( DEF_VEG_SNOW ) THEN", start)
    block = source[start:end]
    assert source.index("ldew(i) = max(0., ldew(i)-evplwet(i)*deltim)") < start

    driver = tmp_path / "pc_phase_storage_probe.f90"
    executable = tmp_path / "pc_phase_storage_probe"
    driver.write_text(
        f"""program pc_phase_storage_probe
  implicit none
  integer, parameter :: r8 = kind(1.d0)
  integer :: i
  logical :: DEF_VEG_SNOW
  real(r8), parameter :: tfrz = 273.15_r8
  real(r8) :: tl(5), ldew(5), ldew_rain(5), ldew_snow(5)

  tl = [tfrz+1._r8, tfrz-1._r8, tfrz+1._r8, tfrz-1._r8, tfrz+1._r8]
  ldew = [2._r8, 3._r8, 0._r8, 0._r8, 9._r8]
  ldew_rain = [-1._r8, -1._r8, -1._r8, -1._r8, 7._r8]
  ldew_snow = [-1._r8, -1._r8, -1._r8, -1._r8, 8._r8]

  DEF_VEG_SNOW = .false.
  do i = 1, 4
{block}
  end do
  if (ldew_rain(1) /= 2._r8 .or. ldew_snow(1) /= 0._r8) error stop 1
  if (ldew_rain(2) /= 0._r8 .or. ldew_snow(2) /= 3._r8) error stop 2
  if (any(ldew_rain(3:4) /= 0._r8) .or. any(ldew_snow(3:4) /= 0._r8)) error stop 3
  if (any(ldew_rain(1:4) + ldew_snow(1:4) /= ldew(1:4))) error stop 4

  DEF_VEG_SNOW = .true.
  i = 5
{block}
  if (ldew_rain(5) /= 7._r8 .or. ldew_snow(5) /= 8._r8) error stop 5
  print '(A)', 'PC_PHASE_STORAGE_OK'
end program pc_phase_storage_probe
""",
        encoding="utf-8",
    )
    built = subprocess.run(
        [compiler, "-ffree-line-length-none", str(driver), "-o", str(executable)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=SMOKE_TIMEOUT,
    )
    assert built.returncode == 0, built.stdout + built.stderr
    ran = subprocess.run(
        [str(executable)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=SMOKE_TIMEOUT,
    )
    assert ran.returncode == 0, ran.stdout + ran.stderr
    assert ran.stdout.strip() == "PC_PHASE_STORAGE_OK"


def test_dew_reclassification_preserves_total_canopy_vapor_flux():
    etr = -4.011370168090042e-15
    evplwet = -1.106845234506008e-9
    total_before = etr + evplwet

    evplwet += etr
    etr = 0.0

    assert etr == 0.0
    assert etr + evplwet == total_before


def test_existing_canopy_capacity_clip_preserves_total_energy():
    latent_heat = 2.5e6
    etr, evplwet, elwmax, fsenl = -2.0, 12.0, 10.0, 3.0
    fevpl = etr + evplwet
    energy_before_clip = fsenl + latent_heat * fevpl

    evplwet += etr
    etr = 0.0
    elwdif = max(0.0, evplwet - elwmax)
    evplwet = min(evplwet, elwmax)
    fevpl -= elwdif
    fsenl += latent_heat * elwdif

    assert evplwet <= elwmax
    assert fsenl + latent_heat * fevpl == energy_before_clip
