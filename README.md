# urbanMicroclimateFoam

An open-source solver for coupled physical processes modeling urban microclimate based on OpenFOAM.

urbanMicroclimateFoam is a multi-region solver consisting of an air subdomain together with subdomains for porous urban building materials. A computational fluid dynamics (CFD) model solves the turbulent, convective air flow and heat and moisture transport in the air subdomain. A coupled heat and moisture (HAM) transport model solves the absorption, transport and storage of heat and moisture in the porous building materials. A radiation model determines the net longwave and shortwave radiative heat fluxes for each surface using a view factor approach.

| Daily variation of surface temperature in a street canyon | Daily variation of air temperature and wind speed at pedestrian height in Münsterhof, Zürich. |
|:---:|:---:|
| <img src="https://gitlab.ethz.ch/openfoam-cbp/solvers/urbanmicroclimatefoam/-/wikis/uploads/9fb2efac6b6827fa7604c3f58960093f/img1_out.gif"  width="75%"> | <img src="https://gitlab.ethz.ch/openfoam-cbp/solvers/urbanmicroclimatefoam/-/wikis/uploads/95d6a8f84991e20b1223c307ca814fc9/img2b_out.gif"  width="45%"> &nbsp; <img src="https://gitlab.ethz.ch/openfoam-cbp/solvers/urbanmicroclimatefoam/-/wikis/uploads/876fcd8fb7d6b18077f9ddeab44db013/img2c_out.gif"  width="45%"> |

The solver is tested for the following OpenFOAM versions:

* OpenFOAM-org (OpenFOAM Foundation) v6

### Usage

You can compile the solver as follows:

	git clone https://gitlab.ethz.ch/openfoam-cbp/solvers/urbanmicroclimatefoam.git
	cd urbanmicroclimatefoam
	./Allwmake
	
### Documentation

Read [urbanMicroclimateFoam wiki](https://gitlab.ethz.ch/openfoam-cbp/solvers/urbanmicroclimatefoam/-/wikis/home) for documentation.

### Tutorial cases

Tutorial case for the solver can be found [here](https://gitlab.ethz.ch/openfoam-cbp/tutorials/urbanmicroclimatefoam).

More information at the Chair of Building Physics: https://carmeliet.ethz.ch
