
   Info<< "\nCalculating stress field\n" << endl;
        graddD = fvc::grad(dD);
        volSymmTensorField sigma_old = (1-alpha)*sigmaD;
        Info<< "\nMax/Min sigma_old\t" <<mag(max(sigma_old()).value())<<'\t'
                <<mag(min(sigma_old()).value())<<endl;
        sigmaD = sigma_old +(1-alpha)*(2*muE*(symm(graddD)-deps_eigen-deps_p) +lambdaE*I*(tr(symm(graddD))-tr(deps_eigen)-tr(deps_p)));
               

       sqrtJ2 = sqrt ((1.0/2.0)*magSqr(dev(sigmaD)));

      //volScalarField ky4 = Foam::max(ky*(-0.5843*T*sig3+1160)*1e6,ky3);
        volScalarField ky2 = Foam::max(ky*(-0.8486*T*sig3+1252)*1e6,ky*(-0.18*T*sig3+316)*1e6);  //Yield-strength
       volScalarField ky4 = Foam::max(ky2,ky3);
       volScalarField fac = sqrtJ2/(ky4);
       Info<< "\nJ2 sqrt\n" <<max(sqrtJ2()).value()<<endl;
       Info<< "\nfac\n" <<max(fac)<<endl;

       forAll(fac,cellI)
       {
           if(mag(fac[cellI]) > 1.0)
           {
               sigmaD[cellI] = 1.0/3.0 *I*tr(sigmaD[cellI]) + dev(sigmaD[cellI])/fac[cellI];
               symmTensor dsigma = sigmaD[cellI]-sigma_old[cellI];
                symmTensor deps_e = 1.0/3.0 *I*tr(dsigma)/(3.0*lambdaE[cellI]+2.0*muE[cellI]) + dev(dsigma)/(2.0*muE[cellI]);
               symmTensor deps   = symm(graddD[cellI]);
               deps_p[cellI] += deps - deps_e - deps_eigen[cellI];
              
           }
             else
             {
                 deps_p[cellI] =symmTensor::zero;
             }
       }

       Info<< "\nAfter plasticity\n" <<endl;
       sqrtJ2 = sqrt ((1.0/2.0)*magSqr(dev(sigmaD)));

       Info<< "\nJ2 sqrt\n" <<max(sqrtJ2()).value()<<endl;

        Info<< "\nMax/Min deps_p magnitude\t" <<mag(max(deps_p()).value())<<'\t'
                <<mag(min(deps_p()).value())<<endl;


