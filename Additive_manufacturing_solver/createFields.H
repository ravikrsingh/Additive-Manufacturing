    // Reading fields
    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
        volScalarField Told
    (
        IOobject
        (
            "Told",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T
    );
        Info<< "\nReading g" << endl;
    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    
    Info<< "Generating/reading field gradT\n" << endl;

    volVectorField gradT

    (

        IOobject

        (

            "gradT",

            runTime.timeName(),

            mesh,

            IOobject::READ_IF_PRESENT,

            IOobject::AUTO_WRITE

        ),

        mesh,

        dimensionSet(0,-1,0,1,0,0,0)

    );

    gradT=fvc::grad(T);

    Info<< "Reading field alpha\n" << endl;
    volScalarField alpha
    (
        IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

    // Reading transport properties
    Info<< "Reading thermophysical properties\n" << endl;
    #include "readTransportProperties.H"

    // Calculating fit-parameters for phase change function
    Info<< "Calculating phase change properties\n" << endl;
    
//     dimensionedScalar Tmelt
//     (
//         "Tmelt", 
//         (Tl+Ts)/2.0
//     );

    // Kinematic density for buoyancy force
    volScalarField rhok
    (
        IOobject
        (
            "rhok",
            runTime.timeName(),
            mesh
        ),
        1.0 - Foam::max(beta*(T - Tl), 0.0)
    );

    // D'arcy-type source term field
    volScalarField DC
    (
        IOobject
        (
            "DC",
            runTime.timeName(),
            mesh
        ),
        DCl*Foam::pow(1.0-alpha,2)/(Foam::pow(alpha,3)+DCs)
    );

    // Thermal conductivity field
    volScalarField lambda
    (
        IOobject
        (
            "lambda",
            runTime.timeName(),
            mesh
        ),
        alpha*(-6.66+1.46e-2*T*sig3)*lambdaL+(1.0-alpha)*(1.57 + 1.6e-2*T*sig3- 1e-6*T*T*sig3*sig3)*lambdaS  
    );

    // Heat capacity field
    volScalarField cp
    (
        IOobject
        (
            "cp",
            runTime.timeName(),
            mesh
        ),
        alpha*cpL+(1.0-alpha)*cpS*(492.4 + 0.025*T*sig3-4.18e-6*T*T*sig3*sig3) ////cps=0.52T+271
    );
    

    // Kinematic viscosity field
    volScalarField nu
    (
        IOobject
        (
            "nu",
            runTime.timeName(),
            mesh
        ),
     alpha*nuL+(1.0-alpha)*nuS
    );

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rhok*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    }

        // Calculating fit-parameters for phase change function
   

//dimensionedScalar ky(transportProperties.lookup("ky"));
    Info<< "Reading Lame's coefficients\n" << endl;
    
        volScalarField muE
        (
            IOobject
            (
                "muE",
                runTime.timeName(),
                mesh
            ),
         (1-alpha)*muE1*(0.0000234*T*T*sig3*sig3-0.085*T*sig3+78.3)*1e9
       
        );
    
        
        volScalarField lambdaE
        (
            IOobject
            (
                "lambdaE",
                runTime.timeName(),
                mesh
            ),
            (1-alpha)*lambdaE1*(0.0000359*T*T*sig3*sig3-0.128*T*sig3+117.3)*1e9
        );
        

    Info<< "Reading field D\n" << endl;
volVectorField dD
(
    IOobject
    (
        "dD",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,

        IOobject::AUTO_WRITE
    ), 
    mesh,
    dimensionSet(0,1,0,0,0,0,0)

);

volTensorField graddD
(
    IOobject
    (
        "graddD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ), 
    fvc::grad(dD)
);

volVectorField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ), 
    dD
);

volSymmTensorField deps_eigen
    (
        IOobject
        (
            "deps_eigen",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),

         mesh,
         dimensionedSymmTensor("deps_eigen", dimless, symmTensor::zero)
    );



Info<< "Reading incremental plastic strain field deps_p\n" << endl;
    volSymmTensorField deps_p
    (
        IOobject
        (
            "deps_p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("deps_p", dimless, symmTensor::zero)
    );

    
    
    Info<< "Calculating stress field sigmaD\n" << endl;
volSymmTensorField sigmaD
(
    IOobject
    (
        "sigmaD",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    //2*muEAvg*(symm(fvc::grad(dD))-symm(deps_eigen))  +lambdaEAvg*I*(tr(fvc::grad(dD))-tr(deps_eigen))
);
volScalarField sqrtJ2
(
    IOobject
    (
        "sqrtJ2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sqrt ((1.0/2.0)*magSqr(dev(sigmaD)))
    //2*muEAvg*(symm(fvc::grad(dD))-symm(deps_eigen))  +lambdaEAvg*I*(tr(fvc::grad(dD))-tr(deps_eigen))
);



mesh.setFluxRequired(dD.name());