/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "HamstadBrick.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace buildingMaterialModels
{
    defineTypeNameAndDebug(HamstadBrick, 0);

    addToRunTimeSelectionTable
    (
        buildingMaterialModel,
        HamstadBrick,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buildingMaterialModels::HamstadBrick::HamstadBrick
(
    const word& name,
    const dictionary& buildingMaterialProperties,
    const word& cellZoneModel
    //volScalarField& h,
    //volScalarField& theta,
    //volScalarField& kr,
    //volScalarField& Ch
)
:
    buildingMaterialModel(name, buildingMaterialProperties, cellZoneModel),// h, theta, kr, Ch),
    HamstadBrickCoeffs_(buildingMaterialProperties.subDict(typeName + "Coeffs")),
    rho_("rho", dimensionSet(1, -3, 0, 0, 0), HamstadBrickCoeffs_.lookup("rho")),
    cap_("cap", dimensionSet(0, 2, -2, -1, 0), HamstadBrickCoeffs_.lookup("cap"))
    /*Ks_(HamstadBrickCoeffs_.lookup("Ks")),
    theta_s_(HamstadBrickCoeffs_.lookup("theta_s")),
    theta_r_(HamstadBrickCoeffs_.lookup("theta_r")),
    alpha_(HamstadBrickCoeffs_.lookup("alpha")),
    beta_(HamstadBrickCoeffs_.lookup("beta")),
    gamma_(HamstadBrickCoeffs_.lookup("gamma")),
    A_(HamstadBrickCoeffs_.lookup("A")),
    Ss_(HamstadBrickCoeffs_.lookup("Ss"))*/
{
    
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::buildingMaterialModels::HamstadBrick::read
(
    const dictionary& buildingMaterialProperties
)
{
    buildingMaterialModel::read(buildingMaterialProperties);

    HamstadBrickCoeffs_ = buildingMaterialProperties.subDict(typeName + "Coeffs");

    HamstadBrickCoeffs_.lookup("rho") >> rho_.value();
    HamstadBrickCoeffs_.lookup("cap") >> cap_.value();

    /*HamstadBrickCoeffs_.lookup("Ks") >> Ks_;
    HamstadBrickCoeffs_.lookup("theta_s") >> theta_s_;
    HamstadBrickCoeffs_.lookup("theta_r") >> theta_r_;
    HamstadBrickCoeffs_.lookup("alpha") >> alpha_;
    HamstadBrickCoeffs_.lookup("beta") >> beta_;
    HamstadBrickCoeffs_.lookup("gamma") >> gamma_;
    HamstadBrickCoeffs_.lookup("A") >> A_;
    HamstadBrickCoeffs_.lookup("Ss") >> Ss_;*/

    return true;
}

//- Correct the buildingMaterial moisture content (cell)
void Foam::buildingMaterialModels::HamstadBrick::update_w_C_cell(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label& celli)
{
    List<scalar> reta; reta.setSize(2); reta[0]=-1.25e-5; reta[1]=-1.80e-5;
    List<scalar> retn; retn.setSize(2); retn[0]=1.65e0; retn[1]=6.00e0;
    List<scalar> retm; retm.setSize(2); retm[0]=0.39394e0; retm[1]=0.83333e0;
    List<scalar> retw; retw.setSize(2); retw[0]=0.300e0; retw[1]=0.700e0;
    scalar w_tmp = 0; scalar tmp = 0; scalar C_tmp = 0; scalar tmp2 = 0;    
    for (int i=0; i<=1; i++)
    {
        tmp = pow( (reta[i]*pc.internalField()[celli]) , retn[i] );
        w_tmp = w_tmp + retw[i] / ( pow( (1 + tmp) , retm[i] ));
        tmp2 = pow( (1 + tmp) , retm[i] );
        C_tmp = C_tmp - retw[i]/tmp2 * retm[i]*retn[i]*tmp/((1 + tmp)*pc.internalField()[celli]);   
    }
    w.internalField()[celli] = w_tmp*157;
    Crel.internalField()[celli] = mag( C_tmp*157 );   
}

//- Correct the buildingMaterial moisture content (boundary)
void Foam::buildingMaterialModels::HamstadBrick::update_w_C_boundary(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label patchi, label patchFacei)
{
    List<scalar> reta; reta.setSize(2); reta[0]=-1.25e-5; reta[1]=-1.80e-5;
    List<scalar> retn; retn.setSize(2); retn[0]=1.65e0; retn[1]=6.00e0;
    List<scalar> retm; retm.setSize(2); retm[0]=0.39394e0; retm[1]=0.83333e0;
    List<scalar> retw; retw.setSize(2); retw[0]=0.300e0; retw[1]=0.700e0;
    scalar w_tmp = 0; scalar tmp = 0; scalar C_tmp = 0; scalar tmp2 = 0;    
    for (int i=0; i<=1; i++)
    {
        tmp = pow( (reta[i]*pc.boundaryField()[patchi][patchFacei]) , retn[i] );
        w_tmp = w_tmp + retw[i] / ( pow( (1 + tmp) , retm[i] ));
        tmp2 = pow( (1 + tmp) , retm[i] );
        C_tmp = C_tmp - retw[i]/tmp2 * retm[i]*retn[i]*tmp/((1 + tmp)*pc.boundaryField()[patchi][patchFacei]);   
    } 
    w.boundaryField()[patchi][patchFacei] = w_tmp*157;
    Crel.boundaryField()[patchi][patchFacei] = mag( C_tmp*157 );  
}

//- Correct the buildingMaterial liquid permeability (cell)
void Foam::buildingMaterialModels::HamstadBrick::update_Krel_cell(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label& celli)
{
    scalar logpc = log10(-pc.internalField()[celli]);
    scalar logKl = 0;
    int i;
    double logpc_M[]={1.0, 
        3.28184696, 3.40184696, 3.52184696, 3.64184696, 3.76184696, 3.88184696,
        4.00184696, 4.12184696, 4.24184696, 4.36184696, 4.48184696, 4.60184696, 4.72184696, 4.84184696, 4.96184696,
        5.08184696, 5.20184696, 5.32184696, 5.44184696, 5.56184696, 5.68184696, 5.80184696, 5.92184696,
        6.04184696, 6.16184696, 6.28184696, 6.40184696, 6.52184696, 6.64184696, 6.76184696, 6.88184696,
        7.00184696, 7.12184696, 7.24184696, 7.36184696, 7.48184696, 7.60184696, 7.72184696, 7.84184696, 7.96184696,
        8.08184696, 8.20184696, 8.32184696, 8.44184696, 8.56184696, 8.68184696, 8.80184696, 8.92184696,
        9.04184696, 9.16184696, 
        12.0};
    double logKl_M[]={-8.71974955,
        -8.71974955, -8.71974955, -8.71974955, -8.71974955, -8.71974955, -8.71974955,
        -8.7214439, -8.7214439, -8.7214439, -8.78646748, -8.92297423, -9.37365288, -10.14654938, -11.13719525, -11.77298037,
        -11.88446026, -11.95855648, -12.00081456, -12.09588942, -12.26136355, -12.39405661, -12.46363481, -12.55777242,
        -12.70786679, -13.00242575, -13.21407102, -13.51572314, -13.83399581, -14.10417193, -14.41224594, -14.73604922,
        -15.15797775, -15.34007327, -15.53780407, -15.80537097, -16.02314218, -31.48, -31.16, -30.84, -30.52,
        -30.2, -29.88, -29.56, -29.24, -28.92, -28.6, -28.28, -27.96,
        -27.64, -27.32,
        -27.0};

    if (logpc < scalar(1.0))
    {
        i = 0;
        logKl = logKl_M[i] + (((logKl_M[i+1] - logKl_M[i])/(logpc_M[i+1] - logpc_M[i]))*(logpc - logpc_M[i])) ;
    }
    else if (logpc >= scalar(12.0))
    {
        i = 50;
        logKl = logKl_M[i] + (((logKl_M[i+1] - logKl_M[i])/(logpc_M[i+1] - logpc_M[i]))*(logpc - logpc_M[i])) ;
    }
    else
    {
        for (i=0; i<=50; ++i)
        {
            if ( (logpc_M[i] <= logpc) && (logpc < logpc_M[i+1]) )
            {
                logKl = logKl_M[i] + (((logKl_M[i+1] - logKl_M[i])/(logpc_M[i+1] - logpc_M[i]))*(logpc - logpc_M[i])) ;
                break;
            }
        }
    }
    Krel.internalField()[celli] = pow(10,logKl);
}

//- Correct the buildingMaterial liquid permeability (boundary)
void Foam::buildingMaterialModels::HamstadBrick::update_Krel_boundary(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label patchi, label patchFacei)
{
    scalar logpc = log10(-pc.boundaryField()[patchi][patchFacei]);
    scalar logKl = 0;
    int i;
    double logpc_M[]={1.0, 
        3.28184696, 3.40184696, 3.52184696, 3.64184696, 3.76184696, 3.88184696,
        4.00184696, 4.12184696, 4.24184696, 4.36184696, 4.48184696, 4.60184696, 4.72184696, 4.84184696, 4.96184696,
        5.08184696, 5.20184696, 5.32184696, 5.44184696, 5.56184696, 5.68184696, 5.80184696, 5.92184696,
        6.04184696, 6.16184696, 6.28184696, 6.40184696, 6.52184696, 6.64184696, 6.76184696, 6.88184696,
        7.00184696, 7.12184696, 7.24184696, 7.36184696, 7.48184696, 7.60184696, 7.72184696, 7.84184696, 7.96184696,
        8.08184696, 8.20184696, 8.32184696, 8.44184696, 8.56184696, 8.68184696, 8.80184696, 8.92184696,
        9.04184696, 9.16184696, 
        12.0};
    double logKl_M[]={-8.71974955,
        -8.71974955, -8.71974955, -8.71974955, -8.71974955, -8.71974955, -8.71974955,
        -8.7214439, -8.7214439, -8.7214439, -8.78646748, -8.92297423, -9.37365288, -10.14654938, -11.13719525, -11.77298037,
        -11.88446026, -11.95855648, -12.00081456, -12.09588942, -12.26136355, -12.39405661, -12.46363481, -12.55777242,
        -12.70786679, -13.00242575, -13.21407102, -13.51572314, -13.83399581, -14.10417193, -14.41224594, -14.73604922,
        -15.15797775, -15.34007327, -15.53780407, -15.80537097, -16.02314218, -31.48, -31.16, -30.84, -30.52,
        -30.2, -29.88, -29.56, -29.24, -28.92, -28.6, -28.28, -27.96,
        -27.64, -27.32,
        -27.0};

    if (logpc < scalar(1.0))
    {
        i = 0;
        logKl = logKl_M[i] + (((logKl_M[i+1] - logKl_M[i])/(logpc_M[i+1] - logpc_M[i]))*(logpc - logpc_M[i])) ;
    }
    else if (logpc >= scalar(12.0))
    {
        i = 50;
        logKl = logKl_M[i] + (((logKl_M[i+1] - logKl_M[i])/(logpc_M[i+1] - logpc_M[i]))*(logpc - logpc_M[i])) ;
    }
    else
    {
        for (i=0; i<=50; ++i)
        {
            if ( (logpc_M[i] <= logpc) && (logpc < logpc_M[i+1]) )
            {
                logKl = logKl_M[i] + (((logKl_M[i+1] - logKl_M[i])/(logpc_M[i+1] - logpc_M[i]))*(logpc - logpc_M[i])) ;
                break;
            }
        }
    }
    Krel.boundaryField()[patchi][patchFacei] = pow(10,logKl);
}

//- Correct the buildingMaterial vapor permeability (cell)
void Foam::buildingMaterialModels::HamstadBrick::update_Kv_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, label& celli)
{
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451*1000/(18.01534); 

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]
    scalar relhum = Foam::exp(pc.internalField()[celli]/(rho_l*R_v*T.internalField()[celli])); // relative humidity [-]
    
    scalar tmp = 1 - (w.internalField()[celli]/1.57e2); 
    scalar delta = 2.61e-5 * tmp/(R_v*T.internalField()[celli]*30*(0.503*tmp*tmp + 0.497)); // Water vapour diffusion coefficient "for brick" [s]
    
    K_v.internalField()[celli] = (delta*p_vsat*relhum)/(rho_l*R_v*T.internalField()[celli]);
}

//- Correct the buildingMaterial vapor permeability (boundary)
void Foam::buildingMaterialModels::HamstadBrick::update_Kv_boundary(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, label patchi, label patchFacei)
{
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451*1000/(18.01534); 

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.boundaryField()[patchi][patchFacei] - 5.976*Foam::log(T.boundaryField()[patchi][patchFacei])); // saturation vapour pressure [Pa]
    scalar relhum = Foam::exp(pc.boundaryField()[patchi][patchFacei]/(rho_l*R_v*T.boundaryField()[patchi][patchFacei])); // relative humidity [-]
    
    scalar tmp = 1 - (w.boundaryField()[patchi][patchFacei]/1.57e2); 
    scalar delta = 2.61e-5 * tmp/(R_v*T.boundaryField()[patchi][patchFacei]*30*(0.503*tmp*tmp + 0.497)); // Water vapour diffusion coefficient "for brick" [s]
    
    K_v.boundaryField()[patchi][patchFacei] = (delta*p_vsat*relhum)/(rho_l*R_v*T.boundaryField()[patchi][patchFacei]);
}

//- Correct the buildingMaterial K_pt (cell)
void Foam::buildingMaterialModels::HamstadBrick::update_Kpt_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_pt, label& celli)
{
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451*1000/(18.01534); 
    scalar L_v = 2.5e6;

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]
    //scalar dpsatdt = (7.06627e3/(T.internalField()[celli]*T.internalField()[celli]) - 5.976/T.internalField()[celli]) * p_vsat; // saturation vapour pressure [Pa]
        
    scalar relhum = Foam::exp(pc.internalField()[celli]/(rho_l*R_v*T.internalField()[celli])); // relative humidity [-]
    
    scalar tmp = 1 - (w.internalField()[celli]/1.57e2); 
    scalar delta = 2.61e-5 * tmp/(R_v*T.internalField()[celli]*30*(0.503*tmp*tmp + 0.497)); // Water vapour diffusion coefficient "for brick" [s]

    K_pt.internalField()[celli] = ( (delta*p_vsat*relhum)/(rho_l*R_v*pow(T.internalField()[celli],2)) ) * (rho_l*L_v - pc.internalField()[celli]);
}

//- Correct the buildingMaterial K_pt (boundary)
void Foam::buildingMaterialModels::HamstadBrick::update_Kpt_boundary(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_pt, label patchi, label patchFacei)
{
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451*1000/(18.01534); 
    scalar L_v = 2.5e6;

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.boundaryField()[patchi][patchFacei] - 5.976*Foam::log(T.boundaryField()[patchi][patchFacei])); // saturation vapour pressure [Pa]
    //scalar dpsatdt = (7.06627e3/(T.boundaryField()[patchi][patchFacei]*T.boundaryField()[patchi][patchFacei]) - 5.976/T.boundaryField()[patchi][patchFacei]) * p_vsat; // saturation vapour pressure [Pa]
        
    scalar relhum = Foam::exp(pc.boundaryField()[patchi][patchFacei]/(rho_l*R_v*T.boundaryField()[patchi][patchFacei])); // relative humidity [-]
    
    scalar tmp = 1 - (w.boundaryField()[patchi][patchFacei]/1.57e2); 
    scalar delta = 2.61e-5 * tmp/(R_v*T.boundaryField()[patchi][patchFacei]*30*(0.503*tmp*tmp + 0.497)); // Water vapour diffusion coefficient "for brick" [s]

    K_pt.boundaryField()[patchi][patchFacei] = ( (delta*p_vsat*relhum)/(rho_l*R_v*pow(T.boundaryField()[patchi][patchFacei],2)) ) * (rho_l*L_v - pc.boundaryField()[patchi][patchFacei]);
}

//- Correct the buildingMaterial lambda (cell)
void Foam::buildingMaterialModels::HamstadBrick::update_lambda_cell(const volScalarField& w, volScalarField& lambda, label& celli)
{

    lambda.internalField()[celli] = 0.5+4.5e-3*w.internalField()[celli];
}

//- Correct the buildingMaterial lambda (boundary)
void Foam::buildingMaterialModels::HamstadBrick::update_lambda_boundary(const volScalarField& w, volScalarField& lambda, label patchi, label patchFacei)
{

    lambda.boundaryField()[patchi][patchFacei] = 0.5+4.5e-3*w.boundaryField()[patchi][patchFacei];
}

//*********************************************************** //
