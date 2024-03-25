module Gaussparam
implicit none 
character(len=4)::gaussChargeSpin
character(len=10)::gaussNproc,gaussMem
character(len=80)::gaussMethod
character(len=20)::gaussTitle 
integer::iounitGaussa,iounitGaussb,iounitfchk,iounitGaussParam
contains

subroutine initializeGaussParam()
implicit none 
integer::j
iounitGaussParam = 100
iounitfchk   = 3000
iounitGaussa = 1000
iounitGaussb = 2000
open(unit=iounitGaussParam,file='gauss.param',status='old')
read(iounitGaussParam,'(A)')gaussNproc
read(iounitGaussParam,'(A)')gaussMem
read(iounitGaussParam,*)
read(iounitGaussParam,'(A)')gaussMethod
read(iounitGaussParam,*)
read(iounitGaussParam,'(A)')gaussTitle 
read(iounitGaussParam,*)
read(iounitGaussParam,'(A)')gaussChargeSpin 
endsubroutine 
end module  
