module  constants
implicit none 
real(8),parameter:: electronMassSI       = 9.1093837015d-31  
real(8),parameter:: electronMassAU       = 1.d0              
real(8),parameter:: BohrToAng            = 0.529177210903d0
real(8),parameter:: AuToJoule            = 4.3597447222071d-18   
real(8),parameter:: LightVel             = 299792458 
real(8),parameter:: AuToFemtoSec         = 2.418884325103622d-2 
real(8),parameter:: KbSI                 = 1.380649d-23     
real(8),parameter:: hbarSI              = 1.054571817d-34 
real(8),parameter:: KbAu                 = 1.d0      
real(8),parameter:: hbarAu              = 1.d0     
real(8),parameter:: protonMassSI         = 1.67262192369d-27 
real(8),parameter:: AmuToAu              = 1822.888486217313d0
real(8),parameter:: HarteeToCminv        = 219474.6315
real(8),parameter:: electronMassSIinv    = 1.d0/electronMassSI         
real(8),parameter:: AngToBohr            = 1.d0/BohrToAng     
real(8),parameter:: JouleToAu            = 1.d0/ AuToJoule    
real(8),parameter:: FemtoSecToAu         = 1.d0/AuToFemtoSec       
real(8),parameter:: KbSIinv              = 1.d0/KbSI    
real(8),parameter:: hbarSIinv            = 1.d0/hbarSI       
real(8),parameter:: protonMassSIinv      = 1.d0/protonMassSI
real(8),parameter:: AuToAmu              = 1.d0/AmuToAu
real(8),parameter:: CminvToHartee        = 1.d0/HarteeToCminv
real(8),parameter:: KelvinToAu           = KbSI * JouleToAu
real(8),parameter:: AuToKelvin           = 1.d0 / KelvinToAu
real(8),parameter:: pi                   = 4.d0*datan(1.d0)  
real(8),parameter:: RadianToDegree       = 180.d0/pi  
real(8),parameter:: DegreeToRadian       = 1.d0/RadianToDegree 
real(8),parameter:: AvogadrosConst       = 6.02214076d+23  
real(8),parameter:: elementaryCharge     = 1.602176634d-19 
real(8),parameter:: ProtonElectronMassRatio = 1836.15267343
real(8),parameter:: MolarGasConstantSI   = 8.314462618
real(8),parameter:: HarteeToEV           =  27.211386245988
end module 
