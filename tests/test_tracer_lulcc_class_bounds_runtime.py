"""Run the production land remap through a lifecycle callback with class 0 present."""

from pathlib import Path
import os
import shutil
import subprocess
import sys

import pytest

from fortran_test_support import netcdf_fortran_flags


ROOT = Path(__file__).resolve().parents[1]


def test_lulcc_class_zero_weights_survive_callback(tmp_path):
    build = Path(os.environ.get("COLM_BLD_DIR", ROOT / ".bld")).resolve()
    library = Path(os.environ.get("COLM_LIB", ROOT / "libcolm.a")).resolve()
    compiler = shutil.which("mpif90")
    if not compiler or not (build / "mod_tracer_defs.mod").exists() or not library.exists():
        pytest.skip("MPI CoLM library and module build required")

    includes, libs = netcdf_fortran_flags()
    flags = ["-fopenmp", "-fdefault-real-8", "-ffree-form", "-fcheck=all", "-cpp",
             "-ffree-line-length-0", "-fallow-argument-mismatch", "-I" + str(ROOT / "include"),
             "-I" + str(tmp_path), "-I" + str(build), *includes, "-J" + str(tmp_path)]
    objects = []
    for name in ("MOD_Tracer_Vars", "MOD_Tracer_Lifecycle"):
        obj = tmp_path / (name + ".o")
        run = subprocess.run([compiler, *flags, "-c", str(ROOT / "main/TRACER" / (name + ".F90")),
                              "-o", str(obj)], capture_output=True, text=True, timeout=120)
        assert run.returncode == 0, run.stdout + run.stderr
        objects.append(obj)

    source = tmp_path / "driver.f90"
    source.write_text("""
program lulcc_class_bounds
 use MOD_Precision
 use MOD_Tracer_Defs, only: ntracers,tracers,STATE_OWNER_GENERIC_WATER,FAMILY_ISOTOPE,REACTION_NONE
 use MOD_Tracer_Vars, only: allocate_Tracer_Vars,save_land_tracer_lulcc_state, &
   remap_land_tracer_lulcc_state,trc_wdsrf
 use MOD_Tracer_Lifecycle, only: tracer_lifecycle_hooks_type,tracer_lifecycle_init, &
   register_tracer_provider,tracer_lifecycle_land_remap_lulcc_state
 implicit none
 real(r8) :: pct(1,0:2)
 type(tracer_lifecycle_hooks_type) :: hooks
 integer :: itrc
 ntracers=1
 allocate(tracers(1))
 tracers(1)%name='H2_18O'
 tracers(1)%category='isotope'
 tracers(1)%state_owner=STATE_OWNER_GENERIC_WATER
 tracers(1)%family_id=FAMILY_ISOTOPE
 tracers(1)%init_delta=0._r8
 call allocate_Tracer_Vars(3,-1,1)
 trc_wdsrf(1,:)=[10._r8,20._r8,30._r8]
 call save_land_tracer_lulcc_state()
 pct(1,:)=[0.2_r8,0.3_r8,0.5_r8]
 call tracer_lifecycle_init()
 hooks%land_remap_lulcc => remap_callback
 call register_tracer_provider('H2_18O','','bounds_test',FAMILY_ISOTOPE, &
   STATE_OWNER_GENERIC_WATER,REACTION_NONE,hooks,itrc)
 if (itrc/=1) error stop 98
 call tracer_lifecycle_land_remap_lulcc_state([0],[1_8],[0,1,2],[1_8,1_8,1_8],pct)
 if (abs(trc_wdsrf(1,1)-23._r8)>1.e-10_r8) then
   print *, 'wrong class-weighted pool:',trc_wdsrf(1,1)
   error stop 99
 endif
contains
 subroutine remap_callback(patchclass_new,eindex_new,patchclass_old,eindex_old, &
   lccpct_patches,new_patch_area,old_patch_area)
   integer,intent(in)::patchclass_new(:),patchclass_old(:)
   integer*8,intent(in)::eindex_new(:),eindex_old(:)
   real(r8),intent(in),optional::lccpct_patches(:,0:),new_patch_area(:),old_patch_area(:)
   call remap_land_tracer_lulcc_state(patchclass_new,eindex_new,patchclass_old,eindex_old, &
     lccpct_patches,new_patch_area,old_patch_area)
 end subroutine
end program
""")
    binary = tmp_path / "lulcc_class_bounds"
    math_libs = ["-framework", "Accelerate"] if sys.platform == "darwin" else ["-llapack", "-lblas"]
    link = subprocess.run([compiler, *flags, str(source), *(str(o) for o in objects),
                           str(library), *libs, *math_libs, "-o", str(binary)],
                          capture_output=True, text=True, timeout=120)
    assert link.returncode == 0, link.stdout + link.stderr
    run = subprocess.run([str(binary)], capture_output=True, text=True, timeout=30)
    assert run.returncode == 0, run.stdout + run.stderr
