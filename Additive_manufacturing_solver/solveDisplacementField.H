
//muE=muE1*Foam::max((69.5-0.0535*T*sig3)*1e9,4e9);//*(0.0000234*T*T*sig3*sig3- 0.085*T*sig3+77.3)*1e9;
muE= muE1*Foam::max(1e9*(1-alpha)*(122.7-0.0565*T*sig3)/(2*1.3),4e6); //poisson's ratio=0.3
//lambdaE=lambdaE1*Foam::max((104-0.08*T*sig3)*1e9,6e9);//*(0.0000359*T*T*sig3*sig3-0.128*T*sig3+114)*1e9;
lambdaE=lambdaE1*Foam::max(1e9*(1-alpha)*(122.7-0.0565*T*sig3)*0.3/(1.3*0.4),6e6);

deps_eigen = (1-alpha)*fvc::ddt(T)*runTime.time().deltaT()*Foam::min((4.42e-9*T*sig3+ 7.5e-6),12.5e-6)*threeKalpha;
Info<< "\nCalculating dd field\n" << endl;

do
        {


            {
                fvVectorMatrix dDEqn
                (

                   fvm::laplacian((2*muE +lambdaE), dD, "laplacian(dDdD,dD)")
                    ==
                     -fvc::div
                    (
                      (1-alpha)*(muE*graddD.T() + lambdaE*(I*tr(graddD)) -(muE + lambdaE)*graddD), "div(sigmaExp)"

                    )
                     + fvc::div
                     (
                        (1-alpha)*(2*muE*deps_p + lambdaE*(I*tr(deps_p))),"div(sigmaP)"
                      )
                    +fvc::div
                   (
                        2*muE*deps_eigen + lambdaE*I*tr(deps_eigen),"div(sigmaEigen)"//tr(deps_eigen)=3*beta*deltaT
                    )
                );



                initialResidual = dDEqn.solve().max().initialResidual();
            }
            graddD = fvc::grad(dD);

        } while (initialResidual > convergenceTolerance && ++iCorr < nCorr);

        
        D = D*(1-alpha) + dD;
        Info<< "\nCalculating stress1 field\n" << endl;
       