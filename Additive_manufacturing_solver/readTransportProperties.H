   IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

	dimensionedScalar pi = constant::mathematical::pi;

    // solid -> phase 1
    // liquid -> phase 2

    // Reading density rho
    dimensionedScalar rho(transportProperties.lookup("rho"));


    // Reading thermal conductivity lambda
    dimensionedScalar lambdaS(transportProperties.lookup("lambdaS"));
    dimensionedScalar lambdaL(transportProperties.lookup("lambdaL"));

    // Reading heat capacity cp of solid and liquid phase
    dimensionedScalar cpS(transportProperties.lookup("cpS"));
    dimensionedScalar cpL(transportProperties.lookup("cpL"));

    // Reading kinematic viscosity of solid and liquid
    dimensionedScalar nuS(transportProperties.lookup("nuS"));
    dimensionedScalar nuL(transportProperties.lookup("nuL"));

    // Reading latent heat of fusion hs
    dimensionedScalar hs(transportProperties.lookup("hs"));

    // Reading solid bound of melting temperature Ts
    dimensionedScalar Ts(transportProperties.lookup("Ts"));

    // Reading liquid bound of melting temperature Tl
    dimensionedScalar Tl(transportProperties.lookup("Tl"));

    // Reading volume expansion factor beta
    dimensionedScalar beta(transportProperties.lookup("beta"));

   // Reading large D'arcy-type source term constant DCl
    dimensionedScalar DCl(transportProperties.lookup("DCl"));

   // Reading small D'arcy-type source term constant DCs
    dimensionedScalar DCs(transportProperties.lookup("DCs"));
    
    // Reading melting temperature
    dimensionedScalar Tmelt(transportProperties.lookup("Tmelt"));
    dimensionedScalar constant(transportProperties.lookup("constant"));
    
    // Reading thermal expansion coefficient
    dimensionedSymmTensor threeKalpha(transportProperties.lookup("threeKalpha"));
    
    // Reading yield strength to avoid dividing by 0
    dimensionedScalar ky3(transportProperties.lookup("ky3"));
    
    // Reading Lame's parameters of solid and liquid
    dimensionedScalar muE1(transportProperties.lookup("muE1"));
    dimensionedScalar lambdaE1(transportProperties.lookup("lambdaE1"));
    
    
    // Reading dimensionless parametre
    dimensionedScalar sig3(transportProperties.lookup("sig3"));
    
    // Reading yield strength for dimensionalising
    dimensionedScalar ky(transportProperties.lookup("ky"));
    
    // Reading reference temperature for thermal expansion
    dimensionedScalar Tref(transportProperties.lookup("Tref"));
    
    
//     dimensionedScalar r(transportProperties.lookup("r"));
//     
//     dimensionedScalar P(transportProperties.lookup("P"));
//     dimensionedScalar dimx(transportProperties.lookup("dimx"));
//     dimensionedScalar n(transportProperties.lookup("n"));
//     dimensionedScalar v(transportProperties.lookup("v"));
//     dimensionedScalar H(transportProperties.lookup("H"));

           
