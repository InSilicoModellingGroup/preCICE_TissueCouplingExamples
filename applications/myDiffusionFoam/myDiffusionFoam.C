/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Transient scalar diffusion solver with optional two-way sink feedback
    and optional first-order extracellular decay.

    Legacy (o2Cube) equation:
        d(O2)/dt - laplacian(D, O2) = -O2Sink

    Multi-species mode (tissue_modelling), selected when
    transportProperties contains an entry "species":
        d(Yi)/dt - laplacian(Di, Yi) + ki*Yi = -YiSink

    where ki is optional extracellular decay (k_H2O2, k_NO2, ...; default 0).

    Non-negativity limiter after each solve:
        Yi^(n+1) <- max(Yi^(n+1), 0)
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include <cmath>
#include <cstdint>
#include <cstring>

namespace
{
inline bool sinkValueNeedsClampNoFpe(const Foam::scalar val)
{
    std::uint64_t bits = 0;
    std::memcpy(&bits, &val, sizeof(bits));

    constexpr std::uint64_t signMask = 0x8000000000000000ULL;
    constexpr std::uint64_t expMask = 0x7FF0000000000000ULL;

    const bool nonFinite = ((bits & expMask) == expMask);
    const bool negative = ((bits & signMask) != 0ULL);

    return nonFinite || negative;
}

void clampSinkField(Foam::volScalarField& sink, const Foam::word& name)
{
    sink.correctBoundaryConditions();

    {
        auto& internalField = sink.primitiveFieldRef();
        bool foundBad = false;
        forAll(internalField, cellI)
        {
            if (sinkValueNeedsClampNoFpe(internalField[cellI]))
            {
                internalField[cellI] = 0.0;
                foundBad = true;
            }
        }
        if (foundBad)
        {
            WarningInFunction
                << name << " contained non-finite or negative values "
                   "(possible Aitken NaN). Replaced with 0." << nl;
        }
    }

    sink.max(Foam::dimensionedScalar("zeroSink", sink.dimensions(), 0.0));
}

void solveSpecies
(
    Foam::volScalarField& Y,
    Foam::volScalarField& YSink,
    const Foam::dimensionedScalar& D,
    const Foam::dimensionedScalar& kDecay,
    const bool twoWayCoupling,
    const Foam::word& name
)
{
    using namespace Foam;

    if (!twoWayCoupling)
    {
        YSink.primitiveFieldRef() = 0.0;
        YSink.correctBoundaryConditions();
    }

    clampSinkField(YSink, name + "Sink");

    Info<< name << "Sink stats: min=" << gMin(YSink.primitiveField())
        << " max=" << gMax(YSink.primitiveField())
        << " integral=" << fvc::domainIntegrate(YSink).value() << nl;
    if (kDecay.value() > SMALL)
    {
        Info<< name << " extracellular decay k=" << kDecay.value() << " 1/s" << nl;
    }

    // First-order extracellular decay: + Sp(k, Y)  <=>  -k*Y on the RHS.
    fvScalarMatrix YEqn
    (
        fvm::ddt(Y)
      - fvm::laplacian(D, Y)
      + fvm::Sp(kDecay, Y)
     == -YSink
    );

    YEqn.solve();

    Y.correctBoundaryConditions();
    Y.max(dimensionedScalar("zeroY", Y.dimensions(), 0.0));

    Info<< name << " stats: min=" << gMin(Y.primitiveField())
        << " max=" << gMax(Y.primitiveField())
        << " integral=" << fvc::domainIntegrate(Y).value() << nl << endl;
}
}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient scalar diffusion solver (O2 legacy or multi-species)"
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    const bool twoWayCoupling =
        transportProperties.getOrDefault<bool>("twoWayCoupling", true);

    const scalar configuredWriteInterval =
        runTime.controlDict().getOrDefault<scalar>("writeInterval", 0.01);

    const scalar writeTolerance =
        max(SMALL, 1e-9 * max(scalar(1), configuredWriteInterval));

    scalar nextCoupledWriteTime = GREAT;
    if (twoWayCoupling && configuredWriteInterval > SMALL)
    {
        nextCoupledWriteTime =
            configuredWriteInterval
          * (std::floor((runTime.value() + writeTolerance)/configuredWriteInterval) + 1.0);
    }

    // ------------------------------------------------------------------
    // Multi-species mode: transportProperties { species (H2O2 NO2); ... }
    // ------------------------------------------------------------------
    const bool multiSpecies = transportProperties.found("species");

    PtrList<volScalarField> speciesFields;
    PtrList<volScalarField> speciesSinks;
    PtrList<dimensionedScalar> speciesD;
    PtrList<dimensionedScalar> speciesK;
    wordList speciesNames;

    // Legacy single-field storage (kept alive for the whole run)
    autoPtr<volScalarField> O2Ptr;
    autoPtr<volScalarField> O2SinkPtr;
    autoPtr<dimensionedScalar> DLegacyPtr;
    autoPtr<dimensionedScalar> kLegacyPtr;

    if (multiSpecies)
    {
        speciesNames = wordList(transportProperties.lookup("species"));
        if (speciesNames.empty())
        {
            FatalErrorInFunction
                << "transportProperties 'species' list is empty" << nl
                << exit(FatalError);
        }

        speciesFields.setSize(speciesNames.size());
        speciesSinks.setSize(speciesNames.size());
        speciesD.setSize(speciesNames.size());
        speciesK.setSize(speciesNames.size());

        forAll(speciesNames, i)
        {
            const word& name = speciesNames[i];
            const word DName = "D_" + name;
            const word kName = "k_" + name;
            // Tissue cases use explicit species sink names (H2O2Sink, NO2Sink).
            // The legacy single-field O2 mode below retains O2Sink.
            const word sinkName = name + "Sink";

            if (transportProperties.found(DName))
            {
                speciesD.set
                (
                    i,
                    new dimensionedScalar(DName, dimViscosity, transportProperties)
                );
            }
            else if (transportProperties.found("D"))
            {
                speciesD.set
                (
                    i,
                    new dimensionedScalar("D", dimViscosity, transportProperties)
                );
                Info<< "Using shared D for species " << name << nl;
            }
            else
            {
                FatalErrorInFunction
                    << "Missing diffusivity '" << DName
                    << "' (or shared 'D') for species " << name << nl
                    << exit(FatalError);
            }

            // Optional first-order extracellular decay [1/s]; default 0.
            if (transportProperties.found(kName))
            {
                speciesK.set
                (
                    i,
                    new dimensionedScalar
                    (
                        kName,
                        dimless/dimTime,
                        transportProperties
                    )
                );
            }
            else
            {
                speciesK.set
                (
                    i,
                    new dimensionedScalar(kName, dimless/dimTime, 0.0)
                );
            }

            Info<< "Reading field " << name << nl;
            speciesFields.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        name,
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );

            Info<< "Reading field " << sinkName << nl;
            speciesSinks.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        sinkName,
                        runTime.timeName(),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedScalar
                    (
                        sinkName,
                        speciesFields[i].dimensions()/dimTime,
                        0.0
                    )
                )
            );
        }

        Info<< "\nStarting multi-species diffusion time loop ("
            << speciesNames.size() << " species)\n" << endl;
    }
    else
    {
        DLegacyPtr.reset
        (
            new dimensionedScalar("D", dimViscosity, transportProperties)
        );
        kLegacyPtr.reset
        (
            new dimensionedScalar("k", dimless/dimTime, 0.0)
        );
        if (transportProperties.found("k"))
        {
            kLegacyPtr.reset
            (
                new dimensionedScalar("k", dimless/dimTime, transportProperties)
            );
        }

        Info<< "Reading field O2\n" << endl;
        O2Ptr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "O2",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );

        Info<< "Reading field O2Sink\n" << endl;
        O2SinkPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "O2Sink",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "O2Sink",
                    O2Ptr().dimensions()/dimTime,
                    0.0
                )
            )
        );

        Info<< "\nStarting O2 diffusion time loop\n" << endl;
    }

    while (runTime.loop())
    {
        if (runTime.controlDict().getOrDefault<bool>("adjustTimeStep", false))
        {
            runTime.setDeltaT(runTime.deltaTValue(), true);
        }

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (multiSpecies)
        {
            forAll(speciesNames, i)
            {
                solveSpecies
                (
                    speciesFields[i],
                    speciesSinks[i],
                    speciesD[i],
                    speciesK[i],
                    twoWayCoupling,
                    speciesNames[i]
                );
            }
        }
        else
        {
            solveSpecies
            (
                O2Ptr(),
                O2SinkPtr(),
                DLegacyPtr(),
                kLegacyPtr(),
                twoWayCoupling,
                "O2"
            );
        }

        if (twoWayCoupling && configuredWriteInterval > SMALL)
        {
            if ((runTime.value() + writeTolerance) >= nextCoupledWriteTime)
            {
                Info<< "Writing coupled output at time = "
                    << runTime.timeName() << nl << endl;
                runTime.writeNow();

                while ((runTime.value() + writeTolerance) >= nextCoupledWriteTime)
                {
                    nextCoupledWriteTime += configuredWriteInterval;
                }
            }
        }
        else
        {
            runTime.write();
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
