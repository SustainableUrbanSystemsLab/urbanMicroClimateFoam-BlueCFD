# urbanMicroclimateFoam

An open-source solver for coupled physical processes modeling urban microclimate based on OpenFOAM. Based on original [contributions](https://gitlab.ethz.ch/openfoam-cbp/solvers/urbanmicroclimatefoam) from ETH researchers.

> urbanMicroclimateFoam is a multi-region solver consisting of an air subdomain together with subdomains for porous urban building materials. A computational fluid dynamics (CFD) model solves the turbulent, convective air flow and heat and moisture transport in the air subdomain. A coupled heat and moisture (HAM) transport model solves the absorption, transport and storage of heat and moisture in the porous building materials. A radiation model determines the net longwave and shortwave radiative heat fluxes for each surface using a view factor approach.

---

## Installation with `BlueCFD-2020-1`

1. Start `blueCFD-Core terminal` 
2. Type `git clone https://github.com/Eddy3D-Dev/urbanMicroClimateFoam-BlueCFD`
3. Type `cd urbanMicroClimateFoam-BlueCFD`
4. Type `git checkout of-org_v8.0-bluecfd`
5. Type `./Allwclean && ./Allwmake`
6. It should produce no errors and end compilation in

```
C:/PROGRA~1/BLUECF~1/ofuser-of8/platforms/mingw_w64GccDPInt32Opt/bin/urbanMicroclimateFoam.pdb: cannot load PDB helper DLL
Error occurred with cv2pdb, have stripped binary as a workaround.
```

