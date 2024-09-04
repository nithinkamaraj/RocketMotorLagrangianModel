# RocketMotorLagrangianModel


![dg3-ezgif com-speed](https://github.com/user-attachments/assets/f3fda0ff-4979-4369-afb3-876f40e6baf9)
         
         
   The model is built for openFoam 2206    

Add the download the folder the compile it ( use openFOAM 2206)

## Folder Details 

# RocketMotorParcelFoam

* The model is based on InviscidReactingParcelFoam
* dragwork due to parcels is added as a sourceterm in Energy Equation 
* Viscous term is removed ( high speed compressible flow model)
* the model only works for very low parcel volume fraction
  
# Lagrangian 

* changes kinematicParcel.C calculation of dragWork for individual particle 
* new source term Sw matrix is constructed in kinematicCloudI.H


Particle injection models 

New injection models have been added 

* CellMassFlowRateInjectionModel 
              
  calculates number of parcels per timestep for a cellface (part of the boundary face) and injects ramdomly in the cell. note if parcel count is fraction then the particle count accumlates

* MassFlowRateInjectionRPF
                    
   calculates number of parcels per timestep and injects based on uniform random distribution across the patch face












  
