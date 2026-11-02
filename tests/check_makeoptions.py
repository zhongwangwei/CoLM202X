import os, subprocess, tempfile
from pathlib import Path
root=Path.cwd()
with tempfile.TemporaryDirectory() as tmp:
    p=Path(tmp)
    for version in (9,14):
        compiler=p/f'gfortran-{version}'
        compiler.write_text(f'#!/bin/sh\necho {version}\n'); compiler.chmod(0o755)
        makefile=p/'check.mk'
        makefile.write_text(f'include {root}/include/Makeoptions.github\nall:\n\t@echo "FC=$(FC);FCMP=$(FCMP);MAJOR=$(GCC_VERSION);FLAGS=$(FOPTS)"\n\t@echo "OMPI=$$OMPI_FC;MPICH=$$MPICH_FC;INTEL=$$I_MPI_F90"\n')
        env=dict(os.environ,FC=str(compiler))
        for key in ('OMPI_FC','MPICH_FC','I_MPI_F90','I_MPI_FC'): env.pop(key,None)
        result=subprocess.check_output(['make','-s','-f',str(makefile),f'FF={compiler}'],env=env,text=True)
        assert f'FC={compiler};FCMP={compiler} -fopenmp;MAJOR={version}' in result,result
        assert ('-fallow-argument-mismatch' in result)==(version>=10),result
        assert f'OMPI={compiler};MPICH={compiler};INTEL={compiler}' in result,result
        result=subprocess.check_output(['make','-s','-f',str(makefile),f'FF={compiler}',f'FC={compiler}','OMPI_FC=custom'],env=env,text=True)
        assert 'OMPI=custom;' in result,result
print('PASS GNU9/14 numeric flags, environment/CLI FC, CaMa and MPI backend propagation')
