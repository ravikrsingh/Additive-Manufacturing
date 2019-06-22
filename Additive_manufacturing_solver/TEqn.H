// Solving the energy equation

{
    fvScalarMatrix TEqn
    (
        fvm::ddt(cp, T)
      + fvm::div(phi*fvc::interpolate(cp), T)
      + hs*4.0*exp(-pow(4.0*(T-Tmelt)/(Tl-Ts),2))/Foam::sqrt(pi)/(Tl-Ts)*fvm::ddt(T)
      + hs*4.0*exp(-pow(4.0*(T-Tmelt)/(Tl-Ts),2))/Foam::sqrt(pi)/(Tl-Ts)*(U & fvc::grad(T))
//       +hs*(fvc::div(phi, alpha)+(fvc::ddt(alpha)))
       - fvm::laplacian(lambda/rho, T)

    );

    TEqn.relax();
    TEqn.solve();
    alpha = 0.5*Foam::erf(4.0*(T-Tmelt)/(Tl-Ts))+scalar(0.5);
    /*volScalarField alphaNew = alpha + alphaRelax*(T-((Tl-Ts)*alpha+Ts))*cp/hs; // added
    alpha = max(scalar(0), min(alphaNew, scalar(1)));*/ // added


    cp = alpha*cpL+(1.0-alpha)*cpS*(492.4 + 0.025*T*sig3-4.18e-6*T*T*sig3*sig3);
    lambda = alpha*(-6.66+1.46e-2*T*sig3)*lambdaL+(1.0-alpha)*(1.57 + 1.6e-2*T*sig3- 1e-6*T*T*sig3*sig3)*lambdaS;  
    nu = alpha*nuL+(1.0-alpha)*nuS;
    rhok = 1.0 - Foam::max(beta*(T - Tl),0.0);
    DC = DCl*Foam::pow(1.0-alpha,2)/(Foam::pow(alpha,3)+DCs);
}
Info<< "Min/max T:" << min(T()).value() << ' '
    << max(T()).value() << endl;