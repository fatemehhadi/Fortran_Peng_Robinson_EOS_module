program test_pressure
use mod_hpchem

implicit none

integer, parameter :: PREC=8

real(PREC) :: D,T,massFrac(2)
integer :: nspecs

type(hpchem) :: hpc
!

!... init hectars ...
call hpchem_init(hpc)

nspecs=nSpecies(hpc)  ! number of species
print*,'number of scpecies:',nspecs

! Pressure [Pa], Density [kg/m3], Molar Enthalpy [J/kmol],...
D=50_PREC ![kg/m3]
T=700._PREC ![K]
massFrac(1)=0.4_PREC
massFrac(2)=0.6_PREC

call setState_TDY(hpc,T,D,massFrac)

print*,'p = ',pressure(hpc), 'Pa'
if (pressure(hpc)/=5948368.0652412521_PREC) then
   print*,'Error:pressure is not equal to 5948368.0652412521 Pa'
else
   print*,'Done'
endif 

!----------Answer-----------------
! number of scpecies:           2
! p=    5948368.0652412521      Pa
! Done
!---------------------------------

END
