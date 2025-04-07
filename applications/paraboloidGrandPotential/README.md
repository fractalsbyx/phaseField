# Grand-potential model 
Welcome. This is a manual describing the file structure of the AMMBER implementation of a grand potental-based phase-field model in PRISMS-PF.

# Setting up a simulation
In this AMMBER application, there are four files you need to modify to set up and run a simulation:
- [system.json](system.json)
- [parameters.prm](parameters.prm)
- [BC_AMR.prm](BC_AMR.prm)
- [ICs_and_BCs.cc](ICs_and_BCs.cc)

Most of the model is defined inside system.json, which contains all of the thermodynamic and kinetic parameters of the physical system.

The file parameters.prm is the same as any PRISMS-PF app and controls numerical parameters, I/O parameters, and additional user parameters.

BC_AMR.prm contains the usual PRISMS-PF boundary conditions and adaptive meshing settings for each field. It is kept separately from parameters.prm because it is coupled to the model defined in system.json.

Initial and Dirichlet boundary conditions for each field are defined in ICs_and_BCs.cc. You must recompile after making any changes here.

## system.json

### `dimensions`
Numbers used to non-dimesionalize simulation vaiables.
#### `length_scale`
Real length of 1 simulation length unit, e.g. `"length_scale" : 1.0e-6` means 1 simulation unit is 1.0 micrometers (if your parameters are in meters).
#### `time_scale`
Real duration of 1 simulation time unit, e.g. `"length_scale" : 3600.0` means 1 simulation time unit corresponds to 1 real-time hour (if your parameters are in seconds).

Changes to time scale numerically have the same effect as changes to timestep.
#### `energy_scale`
Real energy density. This has no effect numerically but may be useful for characterization.

### `Vm`
The molar/atomic volume of phases in your system assumed to be constant. This variable has no effect numerically unless `convert_fractional_to_volumetric_energy == true`.
### `l_int`
Width of the diffuse interface.
### `convert_fractional_to_volumetric_energy`
This model uses volumetric free energy, but most free energy representations are given in molar or atomic free energy. If your free energy is provided in molar or atomic units, set this to `true` to scale by `Vm`.
### `solution_component`
Component in the solution that is not explicitly tracked
### `components`
List of chemical components in the system.
### `phases`
Key-value pairs defining the thermodynamic and kinetic parameters of each phase.
>#### `mu_int`
>Interface mobility. Relates driving forces to interface velocity.
>#### `D`
>Diffusivity within the phase.
>#### `sigma`
>Interfacial energy contribution from the phase.
>#### `f_min`
>Minimum free energy in paraboloid representation >[Equation A](#).
>>#### component properties
>>##### `k_well`
>>Curvature of free energy of this phase with respect to the component. [Equation A](#)
>>##### `c_min`
>>Composition of the free energy minimum in this phase. [Equation A](#)
>>##### `x0`
>>Initial composition in this phase.
<br><br>

### `order_parameters`
List of the phase of each order parameter. You can have several order parameters of each phase. They are indexed from 0 in the code and in [ICs_and_BCs.cc](#ics_and_bcscc). For the PRISMS-PF names of the fields, they are automatically assigned the name of the phase + a number. See [BC_AMR.prm](#parametersprm-and-bc_amrprm).



## [parameters.prm](parameters.prm) and [BC_AMR.prm](BC_AMR.prm)
See [Input File](../../doc/doxygen/user_manual/input_file/input_file.md).

In BC_AMR, boundary conditions must be defined for every field in the simulation. 
The names of the fields are as follows:
- `mu_component` for each species present in the system.
- `phase_name_#` for every order parameter starting at 0.

For example, if your order parameters are <br>`[alpha, beta, beta, gamma, alpha]`,<br> the fields will be named <br>`"alpha_0", "alpha_1", "beta_0", "beta_1", "gamma_0"`.

## [ICs_and_BCs.cc](ICs_and_BCs.cc)
See [ICs and BCs](../../doc/doxygen/user_manual/app_files/app_files.md#ics_and_bcscc).
<br>

This application has been set up such that initial compositions constant within each phase and defined by `"x0"` in `system.json`. So only the order parameters need to be set by hand. Set the order parameter initial condition in the vector `eta0`. The index of each order parameter is the same as that defined in `system.json`.

Some abstractions have been provided to make this easier.

The function `interface(x)` returns the correct-width tanh interface profile given a function that matches the level-set near its contour at x = 0. For example,<br> `eta0[0] = interface(r - std::sqrt(x*x + y*y))` creates a circular order-parameter of radius r.

It is also encouraged to define a custom coordiante system (`x,y,z`) from the simulation coordinates (`p[0], p[1], p[2]`) if it is helpful.


# Links
[PRISMS Center Homepage](http://www.prisms-center.org/#/home) <br>
[PRISMS-PF Homepage](https://prisms-center.github.io/phaseField/) <br>
[Code Repository](https://github.com/prisms-center/phaseField) <br>
[User Registration Form](http://goo.gl/forms/GXo7Im8p2Y) <br>
[PRISMS-PF User Forum](https://groups.google.com/forum/#!forum/prisms-pf-users) <br>