# RocketMotorLagrangianModel


![dg3-ezgif com-speed](https://github.com/user-attachments/assets/f3fda0ff-4979-4369-afb3-876f40e6baf9)
         
         
The model is built for openFoam 2206.    

Download the folder, source bashrc file and compile it.

## Folder Details 

# RocketMotorParcelFoam

* The model is based on InviscidReactingParcelFoam
* Dragwork due to parcels is added as a sourceterm in Energy Equation 
* Viscous term is removed ( high speed compressible flow model)
* The model only works for very low parcel volume fraction
  
# Lagrangian 

Drag Work 
* Changes kinematicParcel.C calculation of dragWork for individual particle 
* New drag source term Sw matrix is constructed in kinematicCloudI.H


Particle injection models 

New injection models have been added 

* CellMassFlowRateInjectionModel 
              
  calculates number of parcels per timestep for a cellface (part of the boundary face) and injects ramdomly in the cell. note if parcel count is < 1 then the particle count accumlates

* MassFlowRateInjectionRPF
                    
   calculates number of parcels per timestep and injects based on uniform random distribution function across the patch face












  
