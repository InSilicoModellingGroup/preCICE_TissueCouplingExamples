#include "Adapter.H"
#include "Utilities.H"

#include "IOstreams.H"
#include "Time.H"
#include "fvMesh.H"

#include <algorithm>
#include <exception>
#include <precice/Participant.hpp>

using namespace Foam;

preciceAdapter::Adapter::Adapter(const Time& runTime, const fvMesh& mesh)
: runTime_(runTime),
  mesh_(mesh)
{
    adapterInfo("Loaded the OpenFOAM-preCICE adapter - v1.3.1.", "info");

    // Seed timesteps with the solver's current deltaT to avoid uninitialized zero dt in the first advance()
    timestepSolver_ = runTime_.deltaTValue();
    timestepStored_ = timestepSolver_;

    return;
}

bool preciceAdapter::Adapter::configFileRead()
{
    Info << "Adapter::configFileRead() - START" << endl;
    // We need a try-catch here, as if reading preciceDict fails,
    // the respective exception will be reduced to a warning.
    // See also comment in preciceAdapter::Adapter::configure().
    try
    {
        SETUP_TIMER();
        adapterInfo("Reading preciceDict...", "info");

        // TODO: static is just a quick workaround to be able
        // to find the dictionary also out of scope (e.g. in KappaEffective).
        // We need a better solution.
        static IOdictionary preciceDict(
            IOobject(
                "preciceDict",
                runTime_.system(),
                mesh_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE));

        // Read and display the preCICE configuration file name
        preciceConfigFilename_ = preciceDict.get<fileName>("preciceConfig");
        DEBUG(adapterInfo("  precice-config-file : " + preciceConfigFilename_));

        // Read and display the participant name
        participantName_ = preciceDict.get<word>("participant");
        DEBUG(adapterInfo("  participant name    : " + participantName_));

        // Read and display the list of modules
        DEBUG(adapterInfo("  modules requested   : "));
        auto modules_ = preciceDict.get<wordList>("modules");
        for (const auto& module : modules_)
        {
            DEBUG(adapterInfo("  - " + module + "\n"));

            if (module == "FA" || module == "FP")
            {
                FAenabled_ = true;
            }
            else
            {
                adapterInfo("This report subset supports only the FA module.",
                            "error-deferred");
                return false;
            }
        }

        // Every interface is a subdictionary of "interfaces",
        // each with an arbitrary name. Read all of them and create
        // a list (here: pointer) of dictionaries.
        const auto* interfaceDictPtr = preciceDict.findDict("interfaces");
        DEBUG(adapterInfo("  interfaces : "));

        // Check if we found any interfaces
        // and get the details of each interface
        if (!interfaceDictPtr)
        {
            adapterInfo("  Empty list of interfaces", "warning");
            return false;
        }
        else
        {
            for (const entry& interfaceDictEntry : *interfaceDictPtr)
            {
                if (interfaceDictEntry.isDict())
                {
                    const dictionary& interfaceDict = interfaceDictEntry.dict();
                    struct InterfaceConfig interfaceConfig;

                    interfaceConfig.meshName = interfaceDict.get<word>("mesh");
                    DEBUG(adapterInfo("  - mesh         : " + interfaceConfig.meshName));

                    // By default, assume "faceCenters" as locationsType
                    interfaceConfig.locationsType = interfaceDict.lookupOrDefault<word>("locations", "faceCenters");
                    DEBUG(adapterInfo("    locations    : " + interfaceConfig.locationsType));

                    // By default, assume that no mesh connectivity is required (i.e. no nearest-projection mapping)
                    interfaceConfig.meshConnectivity = interfaceDict.lookupOrDefault<bool>("connectivity", false);
                    // Mesh connectivity only makes sense in case of faceNodes, check and raise a warning otherwise
                    if (interfaceConfig.meshConnectivity && (interfaceConfig.locationsType == "faceCenters" || interfaceConfig.locationsType == "volumeCenters" || interfaceConfig.locationsType == "volumeCentres"))
                    {
                        DEBUG(adapterInfo("Mesh connectivity is not supported for faceCenters or volumeCenters. \n"
                                          "Please configure the desired interface with the locationsType faceNodes. \n"
                                          "Have a look in the adapter documentation for detailed information.",
                                          "warning"));
                        return false;
                    }
                    DEBUG(adapterInfo("    connectivity : " + std::to_string(interfaceConfig.meshConnectivity)));

                    DEBUG(adapterInfo("    patches      : "));
                    auto patches = interfaceDict.get<wordList>("patches");
                    for (auto patch : patches)
                    {
                        interfaceConfig.patchNames.push_back(patch);
                        DEBUG(adapterInfo("      - " + patch));
                    }

                    DEBUG(adapterInfo("    cellSets      : "));
                    auto cellSets = interfaceDict.lookupOrDefault<wordList>("cellSets", wordList());

                    for (auto cellSet : cellSets)
                    {
                        interfaceConfig.cellSetNames.push_back(cellSet);
                        DEBUG(adapterInfo("      - " + cellSet));
                    }

                    if (!interfaceConfig.cellSetNames.empty() && !(interfaceConfig.locationsType == "volumeCenters" || interfaceConfig.locationsType == "volumeCentres"))
                    {
                        adapterInfo("Cell sets are not supported for locationType != volumeCenters. \n"
                                    "Please configure the desired interface with the locationsType volumeCenters. \n"
                                    "Have a look in the adapter documentation for detailed information.",
                                    "warning");
                        return false;
                    }

                    DEBUG(adapterInfo("    writeData    : "));
                    auto writeData = interfaceDict.get<wordList>("writeData");
                    for (auto writeDatum : writeData)
                    {
                        interfaceConfig.writeData.push_back(writeDatum);
                        DEBUG(adapterInfo("      - " + writeDatum));
                    }

                    DEBUG(adapterInfo("    readData     : "));
                    auto readData = interfaceDict.get<wordList>("readData");
                    for (auto readDatum : readData)
                    {
                        interfaceConfig.readData.push_back(readDatum);
                        DEBUG(adapterInfo("      - " + readDatum));
                    }
                    interfacesConfig_.push_back(interfaceConfig);
                }
            }
        }

        if (FAenabled_)
        {
            Info << "[PRINT] Adapter::configure() - FA module IS enabled" << endl;
            Info << "[PRINT] Adapter::configure() - About to create FA::FluidParticle" << endl;
            FA_ = new FA::FluidParticle(mesh_);
            Info << "[PRINT] Adapter::configure() - Successfully created FA::FluidParticle" << endl;
            Info << "[PRINT] Adapter::configure() - About to call FA_->configure()" << endl;
            if (!FA_->configure(preciceDict))
            {
                Info << "[PRINT] Adapter::configure() - FA_->configure() returned false" << endl;
                return false;
            }
            Info << "[PRINT] Adapter::configure() - Successfully called FA_->configure()" << endl;
        }

        // NOTE: Create your module and read any options specific to it here

        if (!FAenabled_)
        {
            adapterInfo("No module is enabled.", "error-deferred");
            return false;
        }

        // TODO: Loading modules should be implemented in more general way,
        // in order to avoid code duplication. See issue #16 on GitHub.

        ACCUMULATE_TIMER(timeInConfigRead_);
    }
    catch (const Foam::error& e)
    {
        adapterInfo(e.message(), "error-deferred");
        Info << "Adapter::configFileRead() - Caught unknown exception" << endl; // <-- PRINT
        return false;
    }
    Info << "Adapter::configFileRead() - END" << endl; // <-- PRINT STATEMENT ADDED
    return true;
}

void preciceAdapter::Adapter::configure()
{
    Info << "[PRINT] Adapter::configure() - START" << endl; // <-- ADDED
    // Read the adapter's configuration file
    if (!configFileRead())
    {
        Info << "[PRINT] Adapter::configure() - configFileRead() returned false" << endl; // <-- ADDED
        errorsInConfigure = true;
        return;
    }

     // Check if configFileRead set errorsInConfigure (it might catch errors and return true but set flag)
    if (errorsInConfigure) {
         Info << "[PRINT] Adapter::configure() - errorsInConfigure was set true during configFileRead()" << endl; // <-- ADDED
         return; // Exit before trying to create participant
    }

    Info << "[PRINT] Adapter::configure() - configFileRead() completed OK" << endl; // <-- ADDED

    try
    {
        // Check the timestep type (fixed vs adjustable)
        DEBUG(adapterInfo("Checking the timestep type (fixed vs adjustable)..."));
        adjustableTimestep_ = runTime_.controlDict().lookupOrDefault("adjustTimeStep", false);
        if (adjustableTimestep_) { DEBUG(adapterInfo("  Timestep type: adjustable.")); }
        else { DEBUG(adapterInfo("  Timestep type: fixed.")); }

        // --- Construct preCICE Participant ---
        SETUP_TIMER();
        DEBUG(adapterInfo("Creating the preCICE solver interface..."));
        DEBUG(adapterInfo("  Number of processes: " + std::to_string(Pstream::nProcs())));
        DEBUG(adapterInfo("  MPI rank: " + std::to_string(Pstream::myProcNo())));

        Info << "[PRINT] Adapter::configure() - About to create precice::Participant" << endl; // <-- ADDED
        Info << "[PRINT] Adapter::configure() - Participant Name = " << participantName_ << endl; // <-- ADDED
        Info << "[PRINT] Adapter::configure() - Config File = " << preciceConfigFilename_ << endl; // <-- ADDED
        Info << "[PRINT] Adapter::configure() - Rank = " << Pstream::myProcNo() << ", Size = " << Pstream::nProcs() << endl; // <-- ADDED

        precice_ = new precice::Participant(participantName_.c_str(), preciceConfigFilename_.c_str(), Pstream::myProcNo(), Pstream::nProcs()); // <-- CRASH POINT (using .c_str() for safety)

        // If the line above crashes, the next line will NOT be printed
        Info << "[PRINT] Adapter::configure() - Successfully created precice::Participant" << endl; // <-- ADDED

        DEBUG(adapterInfo("  preCICE solver interface was created."));
        ACCUMULATE_TIMER(timeInPreciceConstruct_);
         // --- End Construct preCICE Participant ---


         // --- Create Interfaces ---
         REUSE_TIMER();
         DEBUG(adapterInfo("Creating interfaces..."));
         Info << "[PRINT] Adapter::configure() - Checking if interfaces need creation (count: " << interfacesConfig_.size() << ")" << endl; // <-- ADDED
         for (uint i = 0; i < interfacesConfig_.size(); i++)
         {
             Info << "[PRINT] Adapter::configure() - Starting interface creation " << i << endl; // <-- ADDED
             std::string namePointDisplacement = "default";
             std::string nameCellDisplacement = "default";
             bool restartFromDeformed = false;

             Interface* interface = new Interface(*precice_, mesh_, interfacesConfig_.at(i).meshName, interfacesConfig_.at(i).locationsType, interfacesConfig_.at(i).patchNames, interfacesConfig_.at(i).cellSetNames, interfacesConfig_.at(i).meshConnectivity, restartFromDeformed, namePointDisplacement, nameCellDisplacement);
             interfaces_.push_back(interface);
             DEBUG(adapterInfo("Interface created on mesh " + interfacesConfig_.at(i).meshName));
             Info << "[PRINT] Adapter::configure() - Interface object " << i << " created for mesh " << interfacesConfig_.at(i).meshName << endl; // <-- ADDED

             // Special handling for wall surface interfaces - register wall geometry writers
            if (interfacesConfig_.at(i).meshName == "WallSurfaceMesh" &&
                interfacesConfig_.at(i).locationsType == "faceCenters" && FAenabled_)
            {
                Info << "[PRINT] Adapter::configure() - Detected wall surface interface, setting up geometry writers..." << endl;
                // Register wall surface writers (face centers + normals + areas)
                if (FA_->addWallSurfaceWriter(interfacesConfig_.at(i).meshName,
                                               interfacesConfig_.at(i).patchNames,
                                               interface))
                {
                    Info << "[PRINT] Adapter::configure() - Wall surface writers registered successfully" << endl;
                }
                else
                {
                    adapterInfo("Failed to register wall surface writers", "warning");
                }
            }
            else if (interfacesConfig_.at(i).meshName == "BoundaryTriangleMesh" &&
                     (interfacesConfig_.at(i).locationsType == "triangleCenters" ||
                      interfacesConfig_.at(i).locationsType == "triangleCentres") &&
                     FAenabled_)
            {
                if (!FA_->addBoundaryTriangleWriters(
                        interfacesConfig_.at(i).patchNames, interface))
                {
                    adapterInfo("Failed to register boundary-triangle writers",
                                "error-deferred");
                }
            }

             DEBUG(adapterInfo("Adding coupling data writers..."));
            Info << "[PRINT] Adapter::configure() - Adding writers for interface " << i << endl; // <-- ADDED

            // Skip writeData loop for WallSurfaceMesh since writers were already registered by addWallSurfaceWriter()
            bool skipWriteData = (interfacesConfig_.at(i).meshName == "WallSurfaceMesh" &&
                                  interfacesConfig_.at(i).locationsType == "faceCenters" && FAenabled_)
                              || (interfacesConfig_.at(i).meshName == "BoundaryTriangleMesh" &&
                                  (interfacesConfig_.at(i).locationsType == "triangleCenters" ||
                                   interfacesConfig_.at(i).locationsType == "triangleCentres") &&
                                  FAenabled_);

            if (!skipWriteData)
            {
                for (uint j = 0; j < interfacesConfig_.at(i).writeData.size(); j++)
                {
                    std::string dataName = interfacesConfig_.at(i).writeData.at(j);
                    Info << "[PRINT] Adapter::configure() - Trying to add writer for: " << dataName << endl; // <-- ADDED
                    unsigned int inModules = 0;
                    if (FAenabled_ && FA_->addWriters(dataName, interface)) { inModules++; Info << "[PRINT] Added FA writer for " << dataName << endl; } // Your module

                    if (inModules == 0) { adapterInfo(/* ... unknown writer error ... */ "error-deferred"); Info << "[PRINT] ERROR: No writer module found for " << dataName << endl;}
                    else if (inModules > 1) { adapterInfo(/* ... multiple writer error ... */ "error-deferred"); Info << "[PRINT] ERROR: Multiple writer modules found for " << dataName << endl;}
                    else { Info << "[PRINT] Adapter::configure() - Successfully added writer for: " << dataName << endl; } // <-- ADDED
                }
            }
            else
            {
                Info << "[PRINT] Adapter::configure() - Skipping writeData loop for WallSurfaceMesh (already registered)" << endl;
            }

             DEBUG(adapterInfo("Adding coupling data readers..."));
             Info << "[PRINT] Adapter::configure() - Adding readers for interface " << i << endl; // <-- ADDED
             for (uint j = 0; j < interfacesConfig_.at(i).readData.size(); j++)
             {
                 std::string dataName = interfacesConfig_.at(i).readData.at(j);
                 Info << "[PRINT] Adapter::configure() - Trying to add reader for: " << dataName << endl; // <-- ADDED
                 unsigned int inModules = 0;
                 if (FAenabled_ && FA_->addReaders(dataName, interface)) { inModules++; Info << "[PRINT] Added FA reader for " << dataName << endl; } // Your module

                 if (inModules == 0) { adapterInfo(/* ... unknown reader error ... */ "error-deferred"); Info << "[PRINT] ERROR: No reader module found for " << dataName << endl; }
                 else if (inModules > 1) { adapterInfo(/* ... multiple reader error ... */ "error-deferred"); Info << "[PRINT] ERROR: Multiple reader modules found for " << dataName << endl; }
                 else { Info << "[PRINT] Adapter::configure() - Successfully added reader for: " << dataName << endl; } // <-- ADDED
             }

             Info << "[PRINT] Adapter::configure() - Creating buffer for interface " << i << endl; // <-- ADDED
             interface->createBuffer();
             Info << "[PRINT] Adapter::configure() - Buffer created for interface " << i << endl; // <-- ADDED
             Info << "[PRINT] Adapter::configure() - Finished configuration for interface " << i << endl; // <-- ADDED
         }
         ACCUMULATE_TIMER(timeInMeshSetup_);
         Info << "[PRINT] Adapter::configure() - Finished interface creation loop" << endl; // <-- ADDED
        // --- End Create Interfaces ---


        // --- Initialize preCICE ---
        if (precice_) { // Check if participant was created successfully
             Info << "[PRINT] Adapter::configure() - Calling initialize()" << endl; // <-- ADDED
             initialize(); // Calls precice_->initialize()
             Info << "[PRINT] Adapter::configure() - initialize() finished" << endl; // <-- ADDED

             if (requiresWritingCheckpoint())
             {
                 Info << "[PRINT] Adapter::configure() - Setting up implicit-coupling checkpointing." << endl;
                 setupCheckpointing();
                 checkpointing_ = true;
                 writeCheckpoint();
             }
             if (!adjustableTimestep_) { adjustSolverTimeStepAndReadData(); }
        } else {
             adapterInfo("precice_ pointer is null after constructor attempt, cannot initialize.", "error-deferred");
             errorsInConfigure = true;
             Info << "[PRINT] Adapter::configure() - preCICE pointer null before initialize()" << endl; // <-- ADDED
             return;
        }
        // --- End Initialize preCICE ---

        // --- Set endTime ---
        adapterInfo(/* ... message ... */ "info");
        const_cast<Time&>(runTime_).setEndTime(GREAT);
        // --- End Set endTime ---

    }
    catch (const Foam::error& e) {
        adapterInfo(e.message(), "error-deferred"); errorsInConfigure = true;
        Info << "[PRINT] Adapter::configure() - Caught Foam::error: " << e.what() << endl; // <-- PRINT + details
    }
    catch (const std::exception& e_std) {
         adapterInfo(e_std.what(), "error-deferred"); errorsInConfigure = true;
         Info << "[PRINT] Adapter::configure() - Caught std::exception: " << e_std.what() << endl; // <-- PRINT
    }
    catch (...) {
         adapterInfo("Caught unknown exception during preCICE participant creation or configuration", "error-deferred"); errorsInConfigure = true;
         Info << "[PRINT] Adapter::configure() - Caught unknown exception" << endl; // <-- PRINT
    }

    Info << "[PRINT] Adapter::configure() - END" << endl; // <-- PRINT STATEMENT ADDED
    return;
}

void preciceAdapter::Adapter::execute()
{
    Info << "[PRINT] Adapter::execute() - START" << endl;
    if (errorsInConfigure)
    {
        Info << "[PRINT] Adapter::execute() - Found errorsInConfigure flag set, calling error." << endl;
        adapterInfo("There was a problem while configuring the adapter. See the log for details.", "error");
        return;
    }

    // Check if preCICE was initialized before proceeding
    if (!precice_ || !preciceInitialized_) {
         Info << "[PRINT] Adapter::execute() - preCICE not initialized, skipping execute." << endl;
         if(!errorsInConfigure) adapterInfo("Execute called but preCICE not initialized!", "error");
         return;
    }

    Info << "[PRINT] Adapter::execute() - Calling writeCouplingData()" << endl;
    writeCouplingData();
    Info << "[PRINT] Adapter::execute() - Finished writeCouplingData()" << endl;

    Info << "[PRINT] Adapter::execute() - Calling advance()" << endl;
    advance();
    Info << "[PRINT] Adapter::execute() - Finished advance()" << endl;

    // Implicit rollback must complete before reading fresh coupling data.
    // Reading before checkpoint restore would be overwritten by the restored sink fields.
    Info << "[PRINT] Adapter::execute() - Checking checkpoints..." << endl;
    if (requiresReadingCheckpoint())
    {
        Info << "[PRINT] Adapter::execute() - Reading checkpoint..." << endl;
        readCheckpoint();
        Info << "[PRINT] Adapter::execute() - Finished reading checkpoint." << endl;
    }

    if (requiresWritingCheckpoint())
    {
         if (!checkpointing_)
         {
             Info << "[PRINT] Adapter::execute() - Late checkpoint setup." << endl;
             setupCheckpointing();
             checkpointing_ = true;
         }
         Info << "[PRINT] Adapter::execute() - Writing checkpoint..." << endl;
         writeCheckpoint();
         Info << "[PRINT] Adapter::execute() - Finished writing checkpoint." << endl;
    }

    // Always refresh/read for the next solver step (fixed and adjustable dt).
    Info << "[PRINT] Adapter::execute() - Refreshing timestep and reading coupling data for next step." << endl;
    adjustSolverTimeStepAndReadData();
    Info << "[PRINT] Adapter::execute() - Finished refresh/read." << endl;

    SETUP_TIMER();
    if (checkpointing_ && isCouplingTimeWindowComplete())
    {
         Info << "[PRINT] Adapter::execute() - Time window complete, checking if results need writing..." << endl;
         if (runTime_.timePath().type() == fileName::DIRECTORY)
         {
             Info << "Writing coupled output at time = "
                  << runTime_.timeName() << nl << endl;
             const_cast<Time&>(runTime_).writeNow();
             Info << "[PRINT] Adapter::execute() - Called writeNow()" << endl;
         }
    }
    ACCUMULATE_TIMER(timeInWriteResults_);

    if (!isCouplingOngoing())
    {
         Info << "[PRINT] Adapter::execute() - Coupling determined not ongoing." << endl;
         adapterInfo("The coupling completed.", "info");
         finalize();
         const_cast<Time&>(runTime_).setEndTime(runTime_.value());
         adapterInfo("...", "info");
         // Removed: runTime_.completed() call since Foam::Time does not have a completed() member.
         Info << "[PRINT] Adapter::execute() - Set solver endTime, coupling finished." << endl;
    }
    Info << "[PRINT] Adapter::execute() - END" << endl;
    return;
}



void preciceAdapter::Adapter::adjustTimeStep()
{
    Info << "[PRINT] Adapter::adjustTimeStep() - START" << endl; // <-- ADDED
    // Only adjust if preCICE is initialized
    if (!precice_ || !preciceInitialized_) {
         Info << "[PRINT] Adapter::adjustTimeStep() - preCICE not initialized, skipping." << endl; // <-- ADDED
         return;
    }
    adjustSolverTimeStepAndReadData();
    Info << "[PRINT] Adapter::adjustTimeStep() - END" << endl; // <-- ADDED
    return;
}

void preciceAdapter::Adapter::readCouplingData(double relativeReadTime)
{
    Info << "[PRINT] Adapter::readCouplingData() - START" << endl; // <-- ADDED
     // Only read if preCICE is initialized and interfaces exist
    if (!precice_ || !preciceInitialized_ || interfaces_.empty()) {
         Info << "[PRINT] Adapter::readCouplingData() - preCICE not initialized or no interfaces, skipping." << endl; // <-- ADDED
         return;
    }
    SETUP_TIMER();
    DEBUG(adapterInfo("Reading coupling data..."));

    for (uint i = 0; i < interfaces_.size(); i++)
    {
        Info << "[PRINT] Adapter::readCouplingData() - Calling readCouplingData for interface " << i << endl; // <-- ADDED
        interfaces_.at(i)->readCouplingData(relativeReadTime);
        Info << "[PRINT] Adapter::readCouplingData() - Finished readCouplingData for interface " << i << endl; // <-- ADDED
    }

    ACCUMULATE_TIMER(timeInRead_);
    Info << "[PRINT] Adapter::readCouplingData() - END" << endl; // <-- ADDED
    return;
}

void preciceAdapter::Adapter::writeCouplingData()
{
    Info << "[PRINT] Adapter::writeCouplingData() - START" << endl;
    // Allow an initial write before initialize() if preCICE explicitly requests it.
    // This is required for initialize="true" exchanges (e.g. startup wall data).
    bool allowInitialWrite = false;
    if (precice_ && !preciceInitialized_)
    {
        try
        {
            allowInitialWrite = precice_->requiresInitialData();
        }
        catch (...)
        {
            allowInitialWrite = false;
        }
    }

    // Write if either:
    //  1) preCICE is fully initialized (normal runtime), or
    //  2) preCICE requests initial data before initialize().
    if (!precice_ || interfaces_.empty() || (!preciceInitialized_ && !allowInitialWrite)) {
         Info << "[PRINT] Adapter::writeCouplingData() - preCICE not initialized or no interfaces, skipping." << endl;
         return;
    }
    SETUP_TIMER();
    DEBUG(adapterInfo("Writing coupling data..."));

    // The Interface::writeCouplingData() method doesn't return a value!
    // Don't try to capture return values
    for (uint i = 0; i < interfaces_.size(); i++)
    {
         Info << "[PRINT] Adapter::writeCouplingData() - Calling writeCouplingData for interface " << i << endl;

         // For debugging, check if we're writing temperature in this interface
         for (const auto& dataName : interfacesConfig_.at(i).writeData) {
             if (dataName == "T") {
                 Info << "[PRINT] Adapter::writeCouplingData() - Interface " << i
                      << " will write Temperature data" << endl;
                 break;
             }
         }

         // Call the method - it returns void, don't try to capture a return value
         interfaces_.at(i)->writeCouplingData();

         Info << "[PRINT] Adapter::writeCouplingData() - Interface " << i
              << " completed writing data" << endl;
    }

    Info << "[PRINT] Adapter::writeCouplingData() - Finished writing data for all interfaces" << endl;

    ACCUMULATE_TIMER(timeInWrite_);
    Info << "[PRINT] Adapter::writeCouplingData() - END" << endl;
    return;
}

void preciceAdapter::Adapter::initialize()
{
    Info << "[PRINT] Adapter::initialize() - START" << endl; // <-- ADDED
    // Should only be called if precice_ is valid (checked in configure)
    if (!precice_) {
        Info << "[PRINT] Adapter::initialize() - ERROR: precice_ is null!" << endl; // <-- ADDED
        errorsInConfigure = true; // Mark configuration as failed
        adapterInfo("Cannot initialize, precice_ is null", "error");
        return;
    }

    DEBUG(adapterInfo("Initializing the preCICE solver interface..."));
    SETUP_TIMER();

    if (precice_->requiresInitialData()) {
        Info << "[PRINT] Adapter::initialize() - preCICE requires initial data, calling writeCouplingData()" << endl; // <-- ADDED
        // Only data declared by the configured interfaces may be written.
        // rho/mu, when configured, are regular interface writers; attempting
        // to add disabled fields here produced misleading startup warnings.
        writeCouplingData();
        Info << "[PRINT] Adapter::initialize() - Finished writeCouplingData() for initial data" << endl; // <-- ADDED
    } else {
        Info << "[PRINT] Adapter::initialize() - preCICE does not require initial data." << endl; // <-- ADDED
    }


    DEBUG(adapterInfo("Initializing preCICE data..."));
    Info << "[PRINT] Adapter::initialize() - Calling precice_->initialize()" << endl; // <-- ADDED
    try {
        precice_->initialize();
        preciceInitialized_ = true; // Set flag only on success
        Info << "[PRINT] Adapter::initialize() - Successfully called precice_->initialize()" << endl; // <-- ADDED
    }
    catch(const std::exception& e) {
         adapterInfo(std::string("preCICE initialize failed: ") + e.what(), "error");
         Info << "[PRINT] Adapter::initialize() - Caught exception during precice_->initialize(): " << e.what() << endl; // <-- ADDED
         errorsInConfigure = true; // Mark configuration as failed
         // Do not proceed
         return;
    }
    catch(...) {
         adapterInfo("preCICE initialize failed with unknown exception", "error");
         Info << "[PRINT] Adapter::initialize() - Caught unknown exception during precice_->initialize()" << endl; // <-- ADDED
         errorsInConfigure = true; // Mark configuration as failed
         // Do not proceed
         return;
    }

    ACCUMULATE_TIMER(timeInInitialize_);

    adapterInfo("preCICE was configured and initialized", "info");
    Info << "[PRINT] Adapter::initialize() - END" << endl; // <-- ADDED
    return;
}

void preciceAdapter::Adapter::finalize()
{
    Info << "[PRINT] Adapter::finalize() - START" << endl; // <-- ADDED
    // Original check was: if (NULL != precice_ && preciceInitialized_ && !isCouplingOngoing())
    // Let's simplify: Finalize if initialized, teardown regardless if precice_ exists
    if (precice_ && preciceInitialized_)
    {
        DEBUG(adapterInfo("Finalizing the preCICE solver interface..."));
        Info << "[PRINT] Adapter::finalize() - preCICE is initialized, proceeding." << endl; // <-- ADDED

        SETUP_TIMER();
        Info << "[PRINT] Adapter::finalize() - Calling precice_->finalize()" << endl; // <-- ADDED
        try {
            precice_->finalize();
            Info << "[PRINT] Adapter::finalize() - Successfully called precice_->finalize()" << endl; // <-- ADDED
        } catch (const std::exception& e_std) {
             adapterInfo(std::string("Caught std::exception during preCICE finalize: ") + e_std.what(), "error");
             Info << "[PRINT] Adapter::finalize() - Caught exception during precice_->finalize(): " << e_std.what() << endl; // <-- ADDED
        } catch (...) {
             adapterInfo("Caught unknown exception during preCICE finalize", "error");
             Info << "[PRINT] Adapter::finalize() - Caught unknown exception during precice_->finalize()" << endl; // <-- ADDED
        }
        ACCUMULATE_TIMER(timeInFinalize_);

        preciceInitialized_ = false; // Mark as not initialized anymore

        Info << "[PRINT] Adapter::finalize() - Calling teardown()" << endl; // <-- ADDED
        teardown(); // Teardown resources after finalize attempt
    }
    else if (precice_ && !preciceInitialized_)
    {
        // Created but never initialized (e.g., configure failed after participant creation but before initialize)
        Info << "[PRINT] Adapter::finalize() - preCICE created but not initialized, calling teardown()" << endl; // <-- ADDED
        teardown();
    }
    else
    {
        Info << "[PRINT] Adapter::finalize() - preCICE pointer null or already finalized, skipping." << endl; // <-- ADDED
        // adapterInfo("Could not finalize preCICE.", "error"); // This might be too strong if already finalized normally
    }
    Info << "[PRINT] Adapter::finalize() - END" << endl; // <-- ADDED
    return;
}

void preciceAdapter::Adapter::advance()
{
    Info << "[PRINT] Adapter::advance() - START" << endl; // <-- ADDED
     // Only advance if preCICE is initialized
    if (!precice_ || !preciceInitialized_) {
         Info << "[PRINT] Adapter::advance() - preCICE not initialized, skipping." << endl; // <-- ADDED
         return;
    }

    // Guard against non-positive timestep reaching preCICE
    if (timestepSolver_ <= 0.0)
    {
        adapterInfo("advance() received non-positive dt; falling back to current solver deltaT.", "warning");
        timestepSolver_ = std::max(runTime_.deltaTValue(), 1e-16);
        const_cast<Time&>(runTime_).setDeltaT(timestepSolver_, false);
        Info << "[PRINT] Adapter::advance() - Fallback dt set to " << timestepSolver_ << endl;
    }
    DEBUG(adapterInfo("Advancing preCICE..."));
    Info << "[PRINT] Adapter::advance() - Advancing with timestepSolver_ = " << timestepSolver_ << endl; // <-- ADDED

    SETUP_TIMER();
    Info << "[PRINT] Adapter::advance() - Calling precice_->advance()" << endl; // <-- ADDED
    try {
        precice_->advance(timestepSolver_);
        Info << "[PRINT] Adapter::advance() - Successfully called precice_->advance()" << endl; // <-- ADDED
    } catch (const std::exception& e_std) {
         adapterInfo(std::string("Caught std::exception during preCICE advance: ") + e_std.what(), "error");
         Info << "[PRINT] Adapter::advance() - Caught exception during precice_->advance(): " << e_std.what() << endl; // <-- ADDED
         errorsInConfigure = true; // Treat advance errors as configuration errors to stop simulation
    } catch (...) {
         adapterInfo("Caught unknown exception during preCICE advance", "error");
         Info << "[PRINT] Adapter::advance() - Caught unknown exception during precice_->advance()" << endl; // <-- ADDED
         errorsInConfigure = true;
    }
    ACCUMULATE_TIMER(timeInAdvance_);
    Info << "[PRINT] Adapter::advance() - END" << endl; // <-- ADDED
    return;
}

void preciceAdapter::Adapter::adjustSolverTimeStepAndReadData()
{
    Info << "[PRINT] Adapter::adjustSolverTimeStepAndReadData() - START" << endl; // <-- ADDED
    // Only adjust if preCICE is initialized
    if (!precice_ || !preciceInitialized_) {
         Info << "[PRINT] Adapter::adjustSolverTimeStepAndReadData() - preCICE not initialized, skipping." << endl; // <-- ADDED
         return;
    }
    DEBUG(adapterInfo("Adjusting the solver's timestep..."));

    double timestepSolverDetermined = runTime_.deltaTValue(); // Get current dt
    Info << "[PRINT] Adapter::adjustSolverTimeStepAndReadData() - Solver determined dt = " << timestepSolverDetermined << endl; // <-- ADDED

    // --- Handle fixed timestep logic ---
    if (!adjustableTimestep_)
    {
        // Keep the nominal deltaT from controlDict (seeded at construction).
        // Never adopt runTime_.deltaTValue() after subcycling — fp remainder can
        // collapse to ~1e-13 and stall the solver in millions of micro-steps.
        if (useStoredTimestep_ && timestepStored_ > 0.0)
        {
            timestepSolverDetermined = timestepStored_;
            Info << "[PRINT] Adapter::adjustSolverTimeStepAndReadData() - Using nominal fixed solver dt = "
                 << timestepSolverDetermined << endl;
        }
        else if (timestepSolverDetermined > 0.0)
        {
            timestepStored_ = timestepSolverDetermined;
            useStoredTimestep_ = true;
            Info << "[PRINT] Adapter::adjustSolverTimeStepAndReadData() - Seeding fixed solver dt = "
                 << timestepSolverDetermined << endl;
        }
        else
        {
            adapterInfo("Solver reported non-positive dt in fixed mode and no stored fallback is available", "error");
            errorsInConfigure = true;
            return;
        }
    }
    // --- End fixed timestep logic ---


    // --- Get max timestep from preCICE ---
    double maxPreciceDt = 0.0;
    bool couplingOngoing = precice_->isCouplingOngoing(); // Check coupling status
    Info << "[PRINT] Adapter::adjustSolverTimeStepAndReadData() - Is coupling ongoing? " << (couplingOngoing ? "Yes" : "No") << endl; // <-- ADDED

        if (couplingOngoing) {
            Info << "[PRINT] Adapter::adjustSolverTimeStepAndReadData() - Calling precice_->getMaxTimeStepSize()" << endl; // <-- ADDED
            try {
                maxPreciceDt = precice_->getMaxTimeStepSize();
                Info << "[PRINT] Adapter::adjustSolverTimeStepAndReadData() - preCICE max dt = " << maxPreciceDt << endl; // <-- ADDED

                // Guard against zero or negative dt coming from preCICE (causes advance() exception)
                double minDt = 1e-16;
                if (maxPreciceDt <= minDt) {
                    adapterInfo("preCICE returned non-positive max timestep; falling back to solver dt", "warning");
                    maxPreciceDt = std::max(timestepSolverDetermined, minDt);
                    Info << "[PRINT] Adapter::adjustSolverTimeStepAndReadData() - Using fallback dt = " << maxPreciceDt << endl; // <-- ADDED
                }
            } catch (const std::exception& e_std) {
                adapterInfo(std::string("Caught std::exception getting max timestep: ") + e_std.what(), "error");
                Info << "[PRINT] Adapter::adjustSolverTimeStepAndReadData() - Exception getting max dt: " << e_std.what() << endl; // <-- ADDED
                errorsInConfigure = true; return;
            } catch (...) {
                adapterInfo("Caught unknown exception getting max timestep", "error");
                Info << "[PRINT] Adapter::adjustSolverTimeStepAndReadData() - Unknown exception getting max dt" << endl; // <-- ADDED
                errorsInConfigure = true; return;
            }
        } else {
        // If coupling is not ongoing, use the solver's determined step
        timestepSolver_ = timestepSolverDetermined;
        Info << "[PRINT] Adapter::adjustSolverTimeStepAndReadData() - Coupling not ongoing, setting timestepSolver_ = " << timestepSolver_ << endl; // <-- ADDED
        const_cast<Time&>(runTime_).setDeltaT(timestepSolver_, false);
        Info << "[PRINT] Adapter::adjustSolverTimeStepAndReadData() - Set OF dt, skipping readData, END." << endl; // <-- ADDED
         return; // No need to read data if coupling ended
    }
    // --- End get max timestep ---


    // --- Determine final timestep (timestepSolver_) ---
    double tolerance = 1e-14;
    if (!adjustableTimestep_ && useStoredTimestep_ && timestepStored_ > 0.0)
    {
        // Fixed subcycling: nominal deltaT capped by remaining preCICE window time.
        timestepSolver_ = std::min(timestepStored_, maxPreciceDt);
    }
    else if (maxPreciceDt - timestepSolverDetermined > tolerance) {
         adapterInfo(/* ... subcycling message ... */ "info");
         timestepSolver_ = timestepSolverDetermined;
    }
    else if (timestepSolverDetermined - maxPreciceDt > tolerance) {
         adapterInfo(/* ... dt adjustment warning ... */ "warning");
         timestepSolver_ = maxPreciceDt;
    }
    else {
         DEBUG(adapterInfo(/* ... dt same message ... */));
         timestepSolver_ = maxPreciceDt;
    }
    Info << "[PRINT] Adapter::adjustSolverTimeStepAndReadData() - Final calculated timestepSolver_ = " << timestepSolver_ << endl; // <-- ADDED

    // Clamp negative dt
     if (timestepSolver_ < 0) {
         adapterInfo("Calculated negative timestep ("+std::to_string(timestepSolver_)+"), clamping to zero.", "error");
         timestepSolver_ = 0.0;
         errorsInConfigure = true; // Treat as error
     }
    // --- End determine final timestep ---


    // --- Update OpenFOAM timestep ---
    Info << "[PRINT] Adapter::adjustSolverTimeStepAndReadData() - Setting OpenFOAM deltaT to: " << timestepSolver_ << endl; // <-- ADDED
    const_cast<Time&>(runTime_).setDeltaT(timestepSolver_, false);
    // --- End update OpenFOAM timestep ---

    // --- Read Coupling Data ---
    // Serial-implicit endpoint coupling: sample at the remaining window boundary
    // (getMaxTimeStepSize), not the local OpenFOAM substep size.
    DEBUG(adapterInfo("Reading coupling data associated to the calculated time-step size..."));
    if (couplingOngoing && timestepSolver_ > 0) { // Check again coupling status and positive dt
        const double relativeReadTime =
            (maxPreciceDt > 0.0) ? maxPreciceDt : timestepSolver_;
        Info << "[PRINT] Adapter::adjustSolverTimeStepAndReadData() - Calling readCouplingData(relativeReadTime="
             << relativeReadTime << ")" << endl;
         readCouplingData(relativeReadTime);
        Info << "[PRINT] Adapter::adjustSolverTimeStepAndReadData() - Finished readCouplingData()" << endl; // <-- ADDED
    } else {
        Info << "[PRINT] Adapter::adjustSolverTimeStepAndReadData() - Skipping readCouplingData (couplingOngoing="<<couplingOngoing<<", dt="<<timestepSolver_<<")" << endl; // <-- ADDED
        DEBUG(adapterInfo("Skipping readCouplingData as coupling ended or dt is zero.", "debug"));
    }
    // --- End read Coupling Data ---

    Info << "[PRINT] Adapter::adjustSolverTimeStepAndReadData() - END" << endl; // <-- ADDED
    return;
}

bool preciceAdapter::Adapter::isCouplingOngoing()
{
    std::cout << "[DEBUG] isCouplingOngoing() called." << std::endl;
    bool isCouplingOngoing = false;

    // If the coupling ends before the solver ends,
    // the solver would try to access this method again,
    // giving a segmentation fault if precice_
    // was not available.
    if (NULL != precice_)
    {
        std::cout << "[DEBUG] precice_ is not NULL. Calling precice_->isCouplingOngoing()." << std::endl;
        isCouplingOngoing = precice_->isCouplingOngoing();
        std::cout << "[DEBUG] Received isCouplingOngoing = " << isCouplingOngoing << std::endl;
    }
    else
    {
        std::cout << "[DEBUG] precice_ is NULL. Returning false." << std::endl;
    }

    std::cout << "[DEBUG] isCouplingOngoing() returning " << isCouplingOngoing << std::endl;
    return isCouplingOngoing;
}

bool preciceAdapter::Adapter::isCouplingTimeWindowComplete()
{
    std::cout << "[DEBUG] isCouplingTimeWindowComplete() called." << std::endl;
    bool result = precice_->isTimeWindowComplete();
    std::cout << "[DEBUG] isCouplingTimeWindowComplete() returning " << result << std::endl;
    return result;
}

bool preciceAdapter::Adapter::requiresReadingCheckpoint()
{
    std::cout << "[DEBUG] requiresReadingCheckpoint() called." << std::endl;
    bool result = precice_->requiresReadingCheckpoint();
    std::cout << "[DEBUG] requiresReadingCheckpoint() returning " << result << std::endl;
    return result;
}

bool preciceAdapter::Adapter::requiresWritingCheckpoint()
{
    std::cout << "[DEBUG] requiresWritingCheckpoint() called." << std::endl;
    bool result = precice_->requiresWritingCheckpoint();
    std::cout << "[DEBUG] requiresWritingCheckpoint() returning " << result << std::endl;
    return result;
}


void preciceAdapter::Adapter::storeCheckpointTime()
{
    std::cout << "[DEBUG] storeCheckpointTime() called." << std::endl;
    couplingIterationTimeIndex_ = runTime_.timeIndex();
    couplingIterationTimeValue_ = runTime_.value();
    DEBUG(adapterInfo("Stored time value t = " + std::to_string(runTime_.value())));
    std::cout << "[DEBUG] storeCheckpointTime() stored time value t = " << runTime_.value() << std::endl;
    return;
}

void preciceAdapter::Adapter::reloadCheckpointTime()
{
    std::cout << "[DEBUG] reloadCheckpointTime() called." << std::endl;
    const_cast<Time&>(runTime_).setTime(couplingIterationTimeValue_, couplingIterationTimeIndex_);
    // TODO also reset the current iteration?!
    DEBUG(adapterInfo("Reloaded time value t = " + std::to_string(runTime_.value())));
    std::cout << "[DEBUG] reloadCheckpointTime() reloaded time value t = " << runTime_.value() << std::endl;
    return;
}

void preciceAdapter::Adapter::storeMeshPoints()
{
    std::cout << "[DEBUG] storeMeshPoints() called." << std::endl;
    DEBUG(adapterInfo("Storing mesh points..."));
    // TODO: In foam-extend, we would need "allPoints()". Check if this gives the same data.
    meshPoints_ = mesh_.points();
    oldMeshPoints_ = mesh_.oldPoints();

    /*
    // TODO  This is only required for subcycling. It should not be called when not subcycling!!
    // Add a bool 'subcycling' which can be evaluated every timestep.
    if ( !oldVolsStored && mesh_.foundObject<volScalarField::Internal>("V00") ) // For Ddt schemes which use one previous timestep
    {
        setupMeshVolCheckpointing();
        oldVolsStored = true;
    }
    // Update any volume fields from the buffer to the checkpointed values (if already exists.)
    */

    DEBUG(adapterInfo("Stored mesh points."));
    std::cout << "[DEBUG] Mesh points stored." << std::endl;
    if (mesh_.moving())
    {
        std::cout << "[DEBUG] Mesh is moving." << std::endl;
        if (!meshCheckPointed)
        {
            std::cout << "[DEBUG] Mesh not checkpointed yet. Setting up mesh checkpointing." << std::endl;
            // Set up the checkpoint for the mesh flux: meshPhi
            setupMeshCheckpointing();
            meshCheckPointed = true;
        }
        writeMeshCheckpoint();
        std::cout << "[DEBUG] Mesh checkpoint written." << std::endl;
        writeVolCheckpoint(); // Does not write anything unless subcycling.
        std::cout << "[DEBUG] Volume checkpoint written (if applicable)." << std::endl;
    }
}

void preciceAdapter::Adapter::reloadMeshPoints()
{
    std::cout << "[DEBUG] reloadMeshPoints() called." << std::endl;
    if (!mesh_.moving())
    {
        DEBUG(adapterInfo("Mesh points not moved as the mesh is not moving"));
        std::cout << "[DEBUG] Mesh is not moving. Exiting reloadMeshPoints()." << std::endl;
        return;
    }

    // In Foam::polyMesh::movePoints.
    // TODO: The function movePoints overwrites the pointer to the old mesh.
    // Therefore, if you revert the mesh, the oldpointer will be set to the points, which are the new values.
    DEBUG(adapterInfo("Moving mesh points to their previous locations..."));
    std::cout << "[DEBUG] Moving mesh points to previous locations." << std::endl;

    // TODO
    // Switch oldpoints on for pure physics. (is this required?). Switch off for better mesh deformation capabilities?
    // const_cast<pointField&>(mesh_.points()) = oldMeshPoints_;
    const_cast<fvMesh&>(mesh_).movePoints(meshPoints_);
    std::cout << "[DEBUG] Mesh points moved." << std::endl;

    DEBUG(adapterInfo("Moved mesh points to their previous locations."));
    // TODO The if statement can be removed in this case, but it is still included for clarity
    if (meshCheckPointed)
    {
        std::cout << "[DEBUG] Reading mesh checkpoint as meshCheckPointed is true." << std::endl;
        readMeshCheckpoint();
    }

    /*  // TODO This part should only be used when sybcycling. See the description in 'storeMeshPoints()'
        // The if statement can be removed in this case, but it is still included for clarity
    if ( oldVolsStored )
    {
        readVolCheckpoint();
    }
    */
}

void preciceAdapter::Adapter::setupMeshCheckpointing()
{
    std::cout << "[DEBUG] setupMeshCheckpointing() called." << std::endl;
    // The other mesh <type>Fields:
    //      C
    //      Cf
    //      Sf
    //      magSf
    //      delta
    // are updated by the function fvMesh::movePoints. Only the meshPhi needs checkpointing.
    DEBUG(adapterInfo("Creating a list of the mesh checkpointed fields..."));
    std::cout << "[DEBUG] Creating a list of the mesh checkpointed fields." << std::endl;

    // Add meshPhi to the checkpointed fields
    addMeshCheckpointField(
        const_cast<surfaceScalarField&>(
            mesh_.phi()));
#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Added " + mesh_.phi().name() + " in the list of checkpointed fields.");
#endif
}

void preciceAdapter::Adapter::setupMeshVolCheckpointing()
{
    std::cout << "[DEBUG] setupMeshVolCheckpointing() called." << std::endl;
    DEBUG(adapterInfo("Creating a list of the mesh volume checkpointed fields..."));
    std::cout << "[DEBUG] Creating a list of the mesh volume checkpointed fields." << std::endl;
    // Add the V0 and the V00 to the list of checkpointed fields.
    // For V0
    addVolCheckpointField(
        const_cast<volScalarField::Internal&>(
            mesh_.V0()));
#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Added " + mesh_.V0().name() + " in the list of checkpointed fields.");
#endif
    // For V00
    addVolCheckpointField(
        const_cast<volScalarField::Internal&>(
            mesh_.V00()));
#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Added " + mesh_.V00().name() + " in the list of checkpointed fields.");
#endif

    // Also add the buffer fields.
    // TODO For V0
    /* addVolCheckpointFieldBuffer
    (
        const_cast<volScalarField::Internal&>
        (
            mesh_.V0()
        )
    ); */
#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Added " + mesh_.V0().name() + " in the list of buffer checkpointed fields.");
#endif
    // TODO For V00
    /* addVolCheckpointFieldBuffer
    (
        const_cast<volScalarField::Internal&>
        (
            mesh_.V00()
        )
    );*/
#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Added " + mesh_.V00().name() + " in the list of buffer checkpointed fields.");
#endif
}


void preciceAdapter::Adapter::setupCheckpointing()
{
    std::cout << "[DEBUG] setupCheckpointing() called." << std::endl;
    SETUP_TIMER();

    // Add fields in the checkpointing list - sorted for parallel consistency
    DEBUG(adapterInfo("Adding in checkpointed fields..."));
    std::cout << "[DEBUG] Adding checkpointed fields." << std::endl;

#undef doLocalCode
#define doLocalCode(GeomField)                                           \
    /* Checkpoint registered GeomField objects */                        \
    for (const word& obj : mesh_.sortedNames<GeomField>())               \
    {                                                                    \
        addCheckpointField(mesh_.thisDb().getObjectPtr<GeomField>(obj)); \
        DEBUG(adapterInfo("Checkpoint " + obj + " : " #GeomField));      \
        std::cout << "[DEBUG] Checkpoint " << obj << " added for " << #GeomField << std::endl; \
    }

    doLocalCode(volScalarField);
    doLocalCode(volVectorField);
    doLocalCode(volTensorField);
    doLocalCode(volSymmTensorField);

    doLocalCode(surfaceScalarField);
    doLocalCode(surfaceVectorField);
    doLocalCode(surfaceTensorField);

    doLocalCode(pointScalarField);
    doLocalCode(pointVectorField);
    doLocalCode(pointTensorField);

    // NOTE: Add here other object types to checkpoint, if needed.

#undef doLocalCode

    ACCUMULATE_TIMER(timeInCheckpointingSetup_);
    std::cout << "[DEBUG] setupCheckpointing() completed." << std::endl;
}


// All mesh checkpointed fields

void preciceAdapter::Adapter::addMeshCheckpointField(surfaceScalarField& field)
{
    std::cout << "[DEBUG] addMeshCheckpointField(surfaceScalarField) called for field: " << field.name() << std::endl;
    {
        meshSurfaceScalarFields_.push_back(&field);
        meshSurfaceScalarFieldCopies_.push_back(new surfaceScalarField(field));
    }
}

void preciceAdapter::Adapter::addMeshCheckpointField(surfaceVectorField& field)
{
    std::cout << "[DEBUG] addMeshCheckpointField(surfaceVectorField) called for field: " << field.name() << std::endl;
    {
        meshSurfaceVectorFields_.push_back(&field);
        meshSurfaceVectorFieldCopies_.push_back(new surfaceVectorField(field));
    }
}

void preciceAdapter::Adapter::addMeshCheckpointField(volVectorField& field)
{
    std::cout << "[DEBUG] addMeshCheckpointField(volVectorField) called." << std::endl;
    {
        meshVolVectorFields_.push_back(&field);
        meshVolVectorFieldCopies_.push_back(new volVectorField(field));
    }
}

// TODO Internal field for the V0 (volume old) and V00 (volume old-old) fields
void preciceAdapter::Adapter::addVolCheckpointField(volScalarField::Internal& field)
{
    std::cout << "[DEBUG] addVolCheckpointField(volScalarField::Internal) called for field: " << field.name() << std::endl;
    {
        volScalarInternalFields_.push_back(&field);
        volScalarInternalFieldCopies_.push_back(new volScalarField::Internal(field));
    }
}


void preciceAdapter::Adapter::addCheckpointField(volScalarField* field)
{
    if (field)
    {
        std::cout << "[DEBUG] addCheckpointField(volScalarField) called for field: " << field->name() << std::endl;
        volScalarFields_.push_back(field);
        volScalarFieldCopies_.push_back(new volScalarField(*field));
    }
}

void preciceAdapter::Adapter::addCheckpointField(volVectorField* field)
{
    if (field)
    {
        std::cout << "[DEBUG] addCheckpointField(volVectorField) called for field: " << field->name() << std::endl;
        volVectorFields_.push_back(field);
        volVectorFieldCopies_.push_back(new volVectorField(*field));
    }
}

void preciceAdapter::Adapter::addCheckpointField(surfaceScalarField* field)
{
    if (field)
    {
        std::cout << "[DEBUG] addCheckpointField(surfaceScalarField) called for field: " << field->name() << std::endl;
        surfaceScalarFields_.push_back(field);
        surfaceScalarFieldCopies_.push_back(new surfaceScalarField(*field));
    }
}

void preciceAdapter::Adapter::addCheckpointField(surfaceVectorField* field)
{
    if (field)
    {
        std::cout << "[DEBUG] addCheckpointField(surfaceVectorField) called for field: " << field->name() << std::endl;
        surfaceVectorFields_.push_back(field);
        surfaceVectorFieldCopies_.push_back(new surfaceVectorField(*field));
    }
}

void preciceAdapter::Adapter::addCheckpointField(pointScalarField* field)
{
    if (field)
    {
        std::cout << "[DEBUG] addCheckpointField(pointScalarField) called for field: " << field->name() << std::endl;
        pointScalarFields_.push_back(field);
        pointScalarFieldCopies_.push_back(new pointScalarField(*field));
    }
}

void preciceAdapter::Adapter::addCheckpointField(pointVectorField* field)
{
    if (field)
    {
        std::cout << "[DEBUG] addCheckpointField(pointVectorField) called for field: " << field->name() << std::endl;
        pointVectorFields_.push_back(field);
        pointVectorFieldCopies_.push_back(new pointVectorField(*field));
        // TODO: Old time
        // pointVectorFieldsOld_.push_back(const_cast<pointVectorField&>(field->oldTime()));
        // pointVectorFieldCopiesOld_.push_back(new pointVectorField(field->oldTime()));
    }
}

void preciceAdapter::Adapter::addCheckpointField(volTensorField* field)
{
    if (field)
    {
        std::cout << "[DEBUG] addCheckpointField(volTensorField) called for field: " << field->name() << std::endl;
        volTensorFields_.push_back(field);
        volTensorFieldCopies_.push_back(new volTensorField(*field));
    }
}

void preciceAdapter::Adapter::addCheckpointField(surfaceTensorField* field)
{
    if (field)
    {
        std::cout << "[DEBUG] addCheckpointField(surfaceTensorField) called for field: " << field->name() << std::endl;
        surfaceTensorFields_.push_back(field);
        surfaceTensorFieldCopies_.push_back(new surfaceTensorField(*field));
    }
}

void preciceAdapter::Adapter::addCheckpointField(pointTensorField* field)
{
    if (field)
    {
        std::cout << "[DEBUG] addCheckpointField(pointTensorField) called for field: " << field->name() << std::endl;
        pointTensorFields_.push_back(field);
        pointTensorFieldCopies_.push_back(new pointTensorField(*field));
    }
}

void preciceAdapter::Adapter::addCheckpointField(volSymmTensorField* field)
{
    if (field)
    {
        std::cout << "[DEBUG] addCheckpointField(volSymmTensorField) called for field: " << field->name() << std::endl;
        volSymmTensorFields_.push_back(field);
        volSymmTensorFieldCopies_.push_back(new volSymmTensorField(*field));
    }
}


// NOTE: Add here methods to add other object types to checkpoint, if needed.

void preciceAdapter::Adapter::readCheckpoint()
{
    std::cout << "[DEBUG] readCheckpoint() called." << std::endl;
    SETUP_TIMER();

    // TODO: To increase efficiency: only the oldTime() fields of the quantities which are used in the time
    //  derivative are necessary. (In general this is only the velocity). Also old information of the mesh
    //  is required.
    //  Therefore, loading the oldTime() and oldTime().oldTime() fields for the other fields can be excluded
    //  for efficiency.
    DEBUG(adapterInfo("Reading a checkpoint..."));
    std::cout << "[DEBUG] Reading checkpoint..." << std::endl;

    // Reload the runTime
    reloadCheckpointTime();

    // Reload all the fields of type volScalarField
    for (uint i = 0; i < volScalarFields_.size(); i++)
    {
        std::cout << "[DEBUG] Reloading volScalarField index " << i << " for field: " << volScalarFields_.at(i)->name() << std::endl;
        *(volScalarFields_.at(i)) == *(volScalarFieldCopies_.at(i));
        int nOldTimes(volScalarFields_.at(i)->nOldTimes());
        if (nOldTimes >= 1)
        {
            volScalarFields_.at(i)->oldTime() == volScalarFieldCopies_.at(i)->oldTime();
        }
        if (nOldTimes == 2)
        {
            volScalarFields_.at(i)->oldTime().oldTime() == volScalarFieldCopies_.at(i)->oldTime().oldTime();
        }
    }

    // Reload all the fields of type volVectorField
    for (uint i = 0; i < volVectorFields_.size(); i++)
    {
        std::cout << "[DEBUG] Reloading volVectorField index " << i << " for field: " << volVectorFields_.at(i)->name() << std::endl;
        *(volVectorFields_.at(i)) == *(volVectorFieldCopies_.at(i));

        int nOldTimes(volVectorFields_.at(i)->nOldTimes());
        if (nOldTimes >= 1)
        {
            volVectorFields_.at(i)->oldTime() == volVectorFieldCopies_.at(i)->oldTime();
        }
        if (nOldTimes == 2)
        {
            volVectorFields_.at(i)->oldTime().oldTime() == volVectorFieldCopies_.at(i)->oldTime().oldTime();
        }
    }

    // Reload all the fields of type surfaceScalarField
    for (uint i = 0; i < surfaceScalarFields_.size(); i++)
    {
        std::cout << "[DEBUG] Reloading surfaceScalarField index " << i << " for field: " << surfaceScalarFields_.at(i)->name() << std::endl;
        *(surfaceScalarFields_.at(i)) == *(surfaceScalarFieldCopies_.at(i));

        int nOldTimes(surfaceScalarFields_.at(i)->nOldTimes());
        if (nOldTimes >= 1)
        {
            surfaceScalarFields_.at(i)->oldTime() == surfaceScalarFieldCopies_.at(i)->oldTime();
        }
        if (nOldTimes == 2)
        {
            surfaceScalarFields_.at(i)->oldTime().oldTime() == surfaceScalarFieldCopies_.at(i)->oldTime().oldTime();
        }
    }

    // Reload all the fields of type surfaceVectorField
    for (uint i = 0; i < surfaceVectorFields_.size(); i++)
    {
        std::cout << "[DEBUG] Reloading surfaceVectorField index " << i << " for field: " << surfaceVectorFields_.at(i)->name() << std::endl;
        *(surfaceVectorFields_.at(i)) == *(surfaceVectorFieldCopies_.at(i));

        int nOldTimes(surfaceVectorFields_.at(i)->nOldTimes());
        if (nOldTimes >= 1)
        {
            surfaceVectorFields_.at(i)->oldTime() == surfaceVectorFieldCopies_.at(i)->oldTime();
        }
        if (nOldTimes == 2)
        {
            surfaceVectorFields_.at(i)->oldTime().oldTime() == surfaceVectorFieldCopies_.at(i)->oldTime().oldTime();
        }
    }

    // Reload all the fields of type pointScalarField
    for (uint i = 0; i < pointScalarFields_.size(); i++)
    {
        std::cout << "[DEBUG] Reloading pointScalarField index " << i << " for field: " << pointScalarFields_.at(i)->name() << std::endl;
        *(pointScalarFields_.at(i)) == *(pointScalarFieldCopies_.at(i));

        int nOldTimes(pointScalarFields_.at(i)->nOldTimes());
        if (nOldTimes >= 1)
        {
            pointScalarFields_.at(i)->oldTime() == pointScalarFieldCopies_.at(i)->oldTime();
        }
        if (nOldTimes == 2)
        {
            pointScalarFields_.at(i)->oldTime().oldTime() == pointScalarFieldCopies_.at(i)->oldTime().oldTime();
        }
    }

    // Reload all the fields of type pointVectorField
    for (uint i = 0; i < pointVectorFields_.size(); i++)
    {
        std::cout << "[DEBUG] Reloading pointVectorField index " << i << " for field: " << pointVectorFields_.at(i)->name() << std::endl;
        *(pointVectorFields_.at(i)) == *(pointVectorFieldCopies_.at(i));

        int nOldTimes(pointVectorFields_.at(i)->nOldTimes());
        if (nOldTimes >= 1)
        {
            pointVectorFields_.at(i)->oldTime() == pointVectorFieldCopies_.at(i)->oldTime();
        }
        if (nOldTimes == 2)
        {
            pointVectorFields_.at(i)->oldTime().oldTime() == pointVectorFieldCopies_.at(i)->oldTime().oldTime();
        }
    }

    // Reload all the fields of type volTensorField
    for (uint i = 0; i < volTensorFields_.size(); i++)
    {
        std::cout << "[DEBUG] Reloading volTensorField index " << i << " for field: " << volTensorFields_.at(i)->name() << std::endl;
        *(volTensorFields_.at(i)) == *(volTensorFieldCopies_.at(i));

        int nOldTimes(volTensorFields_.at(i)->nOldTimes());
        if (nOldTimes >= 1)
        {
            volTensorFields_.at(i)->oldTime() == volTensorFieldCopies_.at(i)->oldTime();
        }
        if (nOldTimes == 2)
        {
            volTensorFields_.at(i)->oldTime().oldTime() == volTensorFieldCopies_.at(i)->oldTime().oldTime();
        }
    }

    // Reload all the fields of type surfaceTensorField
    for (uint i = 0; i < surfaceTensorFields_.size(); i++)
    {
        std::cout << "[DEBUG] Reloading surfaceTensorField index " << i << " for field: " << surfaceTensorFields_.at(i)->name() << std::endl;
        *(surfaceTensorFields_.at(i)) == *(surfaceTensorFieldCopies_.at(i));

        int nOldTimes(surfaceTensorFields_.at(i)->nOldTimes());
        if (nOldTimes >= 1)
        {
            surfaceTensorFields_.at(i)->oldTime() == surfaceTensorFieldCopies_.at(i)->oldTime();
        }
        if (nOldTimes == 2)
        {
            surfaceTensorFields_.at(i)->oldTime().oldTime() == surfaceTensorFieldCopies_.at(i)->oldTime().oldTime();
        }
    }

    // Reload all the fields of type pointTensorField
    for (uint i = 0; i < pointTensorFields_.size(); i++)
    {
        std::cout << "[DEBUG] Reloading pointTensorField index " << i << " for field: " << pointTensorFields_.at(i)->name() << std::endl;
        *(pointTensorFields_.at(i)) == *(pointTensorFieldCopies_.at(i));

        int nOldTimes(pointTensorFields_.at(i)->nOldTimes());
        if (nOldTimes >= 1)
        {
            pointTensorFields_.at(i)->oldTime() == pointTensorFieldCopies_.at(i)->oldTime();
        }
        if (nOldTimes == 2)
        {
            pointTensorFields_.at(i)->oldTime().oldTime() == pointTensorFieldCopies_.at(i)->oldTime().oldTime();
        }
    }

    // Reload all the fields of type volSymmTensorField
    for (uint i = 0; i < volSymmTensorFields_.size(); i++)
    {
        std::cout << "[DEBUG] Reloading volSymmTensorField index " << i << " for field: " << volSymmTensorFields_.at(i)->name() << std::endl;
        *(volSymmTensorFields_.at(i)) == *(volSymmTensorFieldCopies_.at(i));

        int nOldTimes(volSymmTensorFields_.at(i)->nOldTimes());
        if (nOldTimes >= 1)
        {
            volSymmTensorFields_.at(i)->oldTime() == volSymmTensorFieldCopies_.at(i)->oldTime();
        }
        if (nOldTimes == 2)
        {
            volSymmTensorFields_.at(i)->oldTime().oldTime() == volSymmTensorFieldCopies_.at(i)->oldTime().oldTime();
        }
    }

#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Checkpoint was read. Time = " + std::to_string(runTime_.value()));
#endif

    std::cout << "[DEBUG] readCheckpoint() completed. Current time value: " << runTime_.value() << std::endl;
    ACCUMULATE_TIMER(timeInCheckpointingRead_);

    return;
}
void preciceAdapter::Adapter::writeCheckpoint()
{
    std::cout << "[DEBUG] writeCheckpoint() called." << std::endl;
    SETUP_TIMER();

    DEBUG(adapterInfo("Writing a checkpoint..."));
    std::cout << "[DEBUG] Writing a checkpoint..." << std::endl;

    // Store the runTime
    std::cout << "[DEBUG] Storing checkpoint time." << std::endl;
    storeCheckpointTime();

    // Store all the fields of type volScalarField
    for (uint i = 0; i < volScalarFields_.size(); i++)
    {
        std::cout << "[DEBUG] Storing volScalarField index " << i << " for field: " << volScalarFields_.at(i)->name() << std::endl;
        *(volScalarFieldCopies_.at(i)) == *(volScalarFields_.at(i));
    }

    // Store all the fields of type volVectorField
    for (uint i = 0; i < volVectorFields_.size(); i++)
    {
        std::cout << "[DEBUG] Storing volVectorField index " << i << " for field: " << volVectorFields_.at(i)->name() << std::endl;
        *(volVectorFieldCopies_.at(i)) == *(volVectorFields_.at(i));
    }

    // Store all the fields of type volTensorField
    for (uint i = 0; i < volTensorFields_.size(); i++)
    {
        std::cout << "[DEBUG] Storing volTensorField index " << i << " for field: " << volTensorFields_.at(i)->name() << std::endl;
        *(volTensorFieldCopies_.at(i)) == *(volTensorFields_.at(i));
    }

    // Store all the fields of type volSymmTensorField
    for (uint i = 0; i < volSymmTensorFields_.size(); i++)
    {
        std::cout << "[DEBUG] Storing volSymmTensorField index " << i << " for field: " << volSymmTensorFields_.at(i)->name() << std::endl;
        *(volSymmTensorFieldCopies_.at(i)) == *(volSymmTensorFields_.at(i));
    }

    // Store all the fields of type surfaceScalarField
    for (uint i = 0; i < surfaceScalarFields_.size(); i++)
    {
        std::cout << "[DEBUG] Storing surfaceScalarField index " << i << " for field: " << surfaceScalarFields_.at(i)->name() << std::endl;
        *(surfaceScalarFieldCopies_.at(i)) == *(surfaceScalarFields_.at(i));
    }

    // Store all the fields of type surfaceVectorField
    for (uint i = 0; i < surfaceVectorFields_.size(); i++)
    {
        std::cout << "[DEBUG] Storing surfaceVectorField index " << i << " for field: " << surfaceVectorFields_.at(i)->name() << std::endl;
        *(surfaceVectorFieldCopies_.at(i)) == *(surfaceVectorFields_.at(i));
    }

    // Store all the fields of type surfaceTensorField
    for (uint i = 0; i < surfaceTensorFields_.size(); i++)
    {
        std::cout << "[DEBUG] Storing surfaceTensorField index " << i << " for field: " << surfaceTensorFields_.at(i)->name() << std::endl;
        *(surfaceTensorFieldCopies_.at(i)) == *(surfaceTensorFields_.at(i));
    }

    // Store all the fields of type pointScalarField
    for (uint i = 0; i < pointScalarFields_.size(); i++)
    {
        std::cout << "[DEBUG] Storing pointScalarField index " << i << " for field: " << pointScalarFields_.at(i)->name() << std::endl;
        *(pointScalarFieldCopies_.at(i)) == *(pointScalarFields_.at(i));
    }

    // Store all the fields of type pointVectorField
    for (uint i = 0; i < pointVectorFields_.size(); i++)
    {
        std::cout << "[DEBUG] Storing pointVectorField index " << i << " for field: " << pointVectorFields_.at(i)->name() << std::endl;
        *(pointVectorFieldCopies_.at(i)) == *(pointVectorFields_.at(i));
    }

    // Store all the fields of type pointTensorField
    for (uint i = 0; i < pointTensorFields_.size(); i++)
    {
        std::cout << "[DEBUG] Storing pointTensorField index " << i << " for field: " << pointTensorFields_.at(i)->name() << std::endl;
        *(pointTensorFieldCopies_.at(i)) == *(pointTensorFields_.at(i));
    }

#ifdef ADAPTER_DEBUG_MODE
    adapterInfo("Checkpoint for time t = " + std::to_string(runTime_.value()) + " was stored.");
#endif

    ACCUMULATE_TIMER(timeInCheckpointingWrite_);
    std::cout << "[DEBUG] writeCheckpoint() completed." << std::endl;
    return;
}

void preciceAdapter::Adapter::readMeshCheckpoint()
{
    std::cout << "[DEBUG] readMeshCheckpoint() called." << std::endl;
    DEBUG(adapterInfo("Reading a mesh checkpoint..."));
    std::cout << "[DEBUG] Reading a mesh checkpoint..." << std::endl;

    // Reload all the fields of type mesh surfaceScalarField
    for (uint i = 0; i < meshSurfaceScalarFields_.size(); i++)
    {
        std::cout << "[DEBUG] Reading meshSurfaceScalarField index " << i << " for field: " << meshSurfaceScalarFields_.at(i)->name() << std::endl;
        *(meshSurfaceScalarFields_.at(i)) == *(meshSurfaceScalarFieldCopies_.at(i));

        int nOldTimes = meshSurfaceScalarFields_.at(i)->nOldTimes();
        if (nOldTimes >= 1)
        {
            std::cout << "[DEBUG] Reading oldTime for meshSurfaceScalarField index " << i << std::endl;
            meshSurfaceScalarFields_.at(i)->oldTime() == meshSurfaceScalarFieldCopies_.at(i)->oldTime();
        }
        if (nOldTimes == 2)
        {
            std::cout << "[DEBUG] Reading second oldTime for meshSurfaceScalarField index " << i << std::endl;
            meshSurfaceScalarFields_.at(i)->oldTime().oldTime() == meshSurfaceScalarFieldCopies_.at(i)->oldTime().oldTime();
        }
    }

    // Reload all the fields of type mesh surfaceVectorField
    for (uint i = 0; i < meshSurfaceVectorFields_.size(); i++)
    {
        std::cout << "[DEBUG] Reading meshSurfaceVectorField index " << i << " for field: " << meshSurfaceVectorFields_.at(i)->name() << std::endl;
        *(meshSurfaceVectorFields_.at(i)) == *(meshSurfaceVectorFieldCopies_.at(i));

        int nOldTimes = meshSurfaceVectorFields_.at(i)->nOldTimes();
        if (nOldTimes >= 1)
        {
            std::cout << "[DEBUG] Reading oldTime for meshSurfaceVectorField index " << i << std::endl;
            meshSurfaceVectorFields_.at(i)->oldTime() == meshSurfaceVectorFieldCopies_.at(i)->oldTime();
        }
        if (nOldTimes == 2)
        {
            std::cout << "[DEBUG] Reading second oldTime for meshSurfaceVectorField index " << i << std::endl;
            meshSurfaceVectorFields_.at(i)->oldTime().oldTime() == meshSurfaceVectorFieldCopies_.at(i)->oldTime().oldTime();
        }
    }

    // Reload all the fields of type mesh volVectorField
    for (uint i = 0; i < meshVolVectorFields_.size(); i++)
    {
        std::cout << "[DEBUG] Reading meshVolVectorField index " << i << " for field." << std::endl;
        *(meshVolVectorFields_.at(i)) == *(meshVolVectorFieldCopies_.at(i));

        int nOldTimes = meshVolVectorFields_.at(i)->nOldTimes();
        if (nOldTimes >= 1)
        {
            std::cout << "[DEBUG] Reading oldTime for meshVolVectorField index " << i << std::endl;
            meshVolVectorFields_.at(i)->oldTime() == meshVolVectorFieldCopies_.at(i)->oldTime();
        }
        if (nOldTimes == 2)
        {
            std::cout << "[DEBUG] Reading second oldTime for meshVolVectorField index " << i << std::endl;
            meshVolVectorFields_.at(i)->oldTime().oldTime() == meshVolVectorFieldCopies_.at(i)->oldTime().oldTime();
        }
    }

#ifdef ADAPTER_DEBUG_MODE
    adapterInfo("Mesh checkpoint was read. Time = " + std::to_string(runTime_.value()));
#endif

    std::cout << "[DEBUG] readMeshCheckpoint() completed." << std::endl;
    return;
}

void preciceAdapter::Adapter::writeMeshCheckpoint()
{
    std::cout << "[DEBUG] writeMeshCheckpoint() called." << std::endl;
    DEBUG(adapterInfo("Writing a mesh checkpoint..."));
    std::cout << "[DEBUG] Writing a mesh checkpoint..." << std::endl;

    // Store all the fields of type mesh surfaceScalarField
    for (uint i = 0; i < meshSurfaceScalarFields_.size(); i++)
    {
        std::cout << "[DEBUG] Writing meshSurfaceScalarField index " << i << " for field: " << meshSurfaceScalarFields_.at(i)->name() << std::endl;
        *(meshSurfaceScalarFieldCopies_.at(i)) == *(meshSurfaceScalarFields_.at(i));
    }

    // Store all the fields of type mesh surfaceVectorField
    for (uint i = 0; i < meshSurfaceVectorFields_.size(); i++)
    {
        std::cout << "[DEBUG] Writing meshSurfaceVectorField index " << i << " for field: " << meshSurfaceVectorFields_.at(i)->name() << std::endl;
        *(meshSurfaceVectorFieldCopies_.at(i)) == *(meshSurfaceVectorFields_.at(i));
    }

    // Store all the fields of type mesh volVectorField
    for (uint i = 0; i < meshVolVectorFields_.size(); i++)
    {
        std::cout << "[DEBUG] Writing meshVolVectorField index " << i << std::endl;
        *(meshVolVectorFieldCopies_.at(i)) == *(meshVolVectorFields_.at(i));
    }

#ifdef ADAPTER_DEBUG_MODE
    adapterInfo("Mesh checkpoint for time t = " + std::to_string(runTime_.value()) + " was stored.");
#endif

    std::cout << "[DEBUG] writeMeshCheckpoint() completed." << std::endl;
    return;
}

void preciceAdapter::Adapter::readVolCheckpoint()
{
    std::cout << "[DEBUG] readVolCheckpoint() called." << std::endl;
    DEBUG(adapterInfo("Reading the mesh volumes checkpoint..."));
    std::cout << "[DEBUG] Reading the mesh volumes checkpoint..." << std::endl;

    // Reload all the fields of type mesh volScalarField::Internal
    for (uint i = 0; i < volScalarInternalFields_.size(); i++)
    {
        std::cout << "[DEBUG] Reading volScalarInternalField index " << i << " for field: " << volScalarInternalFields_.at(i)->name() << std::endl;
        *(volScalarInternalFields_.at(i)) = *(volScalarInternalFieldCopies_.at(i));
    }

#ifdef ADAPTER_DEBUG_MODE
    adapterInfo("Mesh volumes were read. Time = " + std::to_string(runTime_.value()));
#endif

    std::cout << "[DEBUG] readVolCheckpoint() completed." << std::endl;
    return;
}

void preciceAdapter::Adapter::writeVolCheckpoint()
{
    std::cout << "[DEBUG] writeVolCheckpoint() called." << std::endl;
    DEBUG(adapterInfo("Writing a mesh volumes checkpoint..."));
    std::cout << "[DEBUG] Writing a mesh volumes checkpoint..." << std::endl;

    // Store all the fields of type mesh volScalarField::Internal
    for (uint i = 0; i < volScalarInternalFields_.size(); i++)
    {
        std::cout << "[DEBUG] Writing volScalarInternalField index " << i << " for field: " << volScalarInternalFields_.at(i)->name() << std::endl;
        *(volScalarInternalFieldCopies_.at(i)) = *(volScalarInternalFields_.at(i));
    }

#ifdef ADAPTER_DEBUG_MODE
    adapterInfo("Mesh volumes checkpoint for time t = " + std::to_string(runTime_.value()) + " was stored.");
#endif

    std::cout << "[DEBUG] writeVolCheckpoint() completed." << std::endl;
    return;
}

void preciceAdapter::Adapter::end()
{
    std::cout << "[DEBUG] end() called." << std::endl;
    // If preCICE is still active and the coupling has not yet finished, keep exchanging
    // data until preCICE reports that all participants are done. This avoids leaving the
    // partner (BioDynaMo) stuck waiting for additional communication.
    if (precice_ && preciceInitialized_ && precice_->isCouplingOngoing()) {
        adapterInfo(
            "The solver requested shutdown while the coupling was still ongoing. "
            "Completing the remaining preCICE exchanges before terminating.",
            "warning");
        std::cout << "[DEBUG] Coupling still ongoing during end(). Completing remaining exchanges."
                  << std::endl;

        const unsigned int kMaxEmergencyIterations = 100000;
        unsigned int iterations = 0;

        while (precice_->isCouplingOngoing()) {
            // Reuse existing utilities to stay in sync with the participant order.
            adjustTimeStep();
            writeCouplingData();
            advance();

            if (requiresReadingCheckpoint()) {
                readCheckpoint();
            }
            if (requiresWritingCheckpoint()) {
                writeCheckpoint();
            }

            if (++iterations >= kMaxEmergencyIterations) {
                adapterInfo(
                    "Reached emergency limit while completing coupling during end(); "
                    "forcing finalize to avoid an infinite loop.",
                    "warning");
                std::cout << "[DEBUG] Emergency limit reached while finalizing coupling."
                          << std::endl;
                break;
            }
        }

        // Make sure we notify preCICE that we are done.
        finalize();
        std::cout << "[DEBUG] Completed coupling cleanup in end() after "
                  << iterations << " iterations." << std::endl;
    } else if (precice_ && !preciceInitialized_) {
        std::cout << "[DEBUG] end() called with preCICE created but not initialized."
                  << std::endl;
    } else if (!precice_) {
        std::cout << "[DEBUG] end() called with no active preCICE participant." << std::endl;
    }

    std::cout << "[DEBUG] end() completed." << std::endl;
    return;
}

void preciceAdapter::Adapter::teardown()
{
    std::cout << "[DEBUG] teardown() called." << std::endl;
    // If the solver interface was not deleted before, delete it now.
    // Normally it should be deleted when isCouplingOngoing() becomes false.
    if (NULL != precice_)
    {
        DEBUG(adapterInfo("Destroying the preCICE solver interface..."));
        std::cout << "[DEBUG] Destroying preCICE solver interface." << std::endl;
        delete precice_;
        precice_ = NULL;
    }

    // Delete the preCICE solver interfaces
    if (interfaces_.size() > 0)
    {
        DEBUG(adapterInfo("Deleting the interfaces..."));
        std::cout << "[DEBUG] Deleting interfaces." << std::endl;
        for (uint i = 0; i < interfaces_.size(); i++)
        {
            std::cout << "[DEBUG] Deleting interface index " << i << std::endl;
            delete interfaces_.at(i);
        }
        interfaces_.clear();
    }

    // Delete the copied fields for checkpointing
    if (checkpointing_)
    {
        DEBUG(adapterInfo("Deleting the checkpoints... "));
        std::cout << "[DEBUG] Deleting checkpoint copies." << std::endl;

        // Fields
        // volScalarFields
        for (uint i = 0; i < volScalarFieldCopies_.size(); i++)
        {
            std::cout << "[DEBUG] Deleting volScalarFieldCopy index " << i << std::endl;
            delete volScalarFieldCopies_.at(i);
        }
        volScalarFieldCopies_.clear();
        // volVector
        for (uint i = 0; i < volVectorFieldCopies_.size(); i++)
        {
            std::cout << "[DEBUG] Deleting volVectorFieldCopy index " << i << std::endl;
            delete volVectorFieldCopies_.at(i);
        }
        volVectorFieldCopies_.clear();
        // surfaceScalar
        for (uint i = 0; i < surfaceScalarFieldCopies_.size(); i++)
        {
            std::cout << "[DEBUG] Deleting surfaceScalarFieldCopy index " << i << std::endl;
            delete surfaceScalarFieldCopies_.at(i);
        }
        surfaceScalarFieldCopies_.clear();
        // surfaceVector
        for (uint i = 0; i < surfaceVectorFieldCopies_.size(); i++)
        {
            std::cout << "[DEBUG] Deleting surfaceVectorFieldCopy index " << i << std::endl;
            delete surfaceVectorFieldCopies_.at(i);
        }
        surfaceVectorFieldCopies_.clear();
        // pointScalar
        for (uint i = 0; i < pointScalarFieldCopies_.size(); i++)
        {
            std::cout << "[DEBUG] Deleting pointScalarFieldCopy index " << i << std::endl;
            delete pointScalarFieldCopies_.at(i);
        }
        pointScalarFieldCopies_.clear();
        // pointVector
        for (uint i = 0; i < pointVectorFieldCopies_.size(); i++)
        {
            std::cout << "[DEBUG] Deleting pointVectorFieldCopy index " << i << std::endl;
            delete pointVectorFieldCopies_.at(i);
        }
        pointVectorFieldCopies_.clear();

        // Mesh fields
        // meshSurfaceScalar
        for (uint i = 0; i < meshSurfaceScalarFieldCopies_.size(); i++)
        {
            std::cout << "[DEBUG] Deleting meshSurfaceScalarFieldCopy index " << i << std::endl;
            delete meshSurfaceScalarFieldCopies_.at(i);
        }
        meshSurfaceScalarFieldCopies_.clear();

        // meshSurfaceVector
        for (uint i = 0; i < meshSurfaceVectorFieldCopies_.size(); i++)
        {
            std::cout << "[DEBUG] Deleting meshSurfaceVectorFieldCopy index " << i << std::endl;
            delete meshSurfaceVectorFieldCopies_.at(i);
        }
        meshSurfaceVectorFieldCopies_.clear();

        // meshVolVector
        for (uint i = 0; i < meshVolVectorFieldCopies_.size(); i++)
        {
            std::cout << "[DEBUG] Deleting meshVolVectorFieldCopy index " << i << std::endl;
            delete meshVolVectorFieldCopies_.at(i);
        }
        meshVolVectorFieldCopies_.clear();

        // volScalarInternal
        for (uint i = 0; i < volScalarInternalFieldCopies_.size(); i++)
        {
            std::cout << "[DEBUG] Deleting volScalarInternalFieldCopy index " << i << std::endl;
            delete volScalarInternalFieldCopies_.at(i);
        }
        volScalarInternalFieldCopies_.clear();

        // volTensorField
        for (uint i = 0; i < volTensorFieldCopies_.size(); i++)
        {
            std::cout << "[DEBUG] Deleting volTensorFieldCopy index " << i << std::endl;
            delete volTensorFieldCopies_.at(i);
        }
        volTensorFieldCopies_.clear();

        // surfaceTensorField
        for (uint i = 0; i < surfaceTensorFieldCopies_.size(); i++)
        {
            std::cout << "[DEBUG] Deleting surfaceTensorFieldCopy index " << i << std::endl;
            delete surfaceTensorFieldCopies_.at(i);
        }
        surfaceTensorFieldCopies_.clear();

        // pointTensorField
        for (uint i = 0; i < pointTensorFieldCopies_.size(); i++)
        {
            std::cout << "[DEBUG] Deleting pointTensorFieldCopy index " << i << std::endl;
            delete pointTensorFieldCopies_.at(i);
        }
        pointTensorFieldCopies_.clear();

        // volSymmTensorField
        for (uint i = 0; i < volSymmTensorFieldCopies_.size(); i++)
        {
            std::cout << "[DEBUG] Deleting volSymmTensorFieldCopy index " << i << std::endl;
            delete volSymmTensorFieldCopies_.at(i);
        }
        volSymmTensorFieldCopies_.clear();

        checkpointing_ = false;
        std::cout << "[DEBUG] Checkpoints deleted." << std::endl;
    }

    // Delete the FA module
    if (NULL != FA_)
    {
        DEBUG(adapterInfo("Destroying the FA module..."));
        std::cout << "[DEBUG] Destroying FA module." << std::endl;
        delete FA_;
        FA_ = NULL;
    }

    // NOTE: Delete your new module here

    std::cout << "[DEBUG] teardown() completed." << std::endl;
    return;
}

preciceAdapter::Adapter::~Adapter()
{
    std::cout << "[DEBUG] ~Adapter() destructor called." << std::endl;
    teardown();

    TIMING_MODE(
        // Continuing the output started in the destructor of preciceAdapterFunctionObject
        Info << "Time exclusively in the adapter: " << (timeInConfigRead_ + timeInMeshSetup_ + timeInCheckpointingSetup_ + timeInWrite_ + timeInRead_ + timeInCheckpointingWrite_ + timeInCheckpointingRead_).str() << nl;
        Info << "  (S) reading preciceDict:       " << timeInConfigRead_.str() << nl;
        Info << "  (S) constructing preCICE:      " << timeInPreciceConstruct_.str() << nl;
        Info << "  (S) setting up the interfaces: " << timeInMeshSetup_.str() << nl;
        Info << "  (S) setting up checkpointing:  " << timeInCheckpointingSetup_.str() << nl;
        Info << "  (I) writing data:              " << timeInWrite_.str() << nl;
        Info << "  (I) reading data:              " << timeInRead_.str() << nl;
        Info << "  (I) writing checkpoints:       " << timeInCheckpointingWrite_.str() << nl;
        Info << "  (I) reading checkpoints:       " << timeInCheckpointingRead_.str() << nl;
        Info << "  (I) writing OpenFOAM results:  " << timeInWriteResults_.str() << " (at the end of converged time windows)" << nl << nl;
        Info << "Time exclusively in preCICE:     " << (timeInInitialize_ + timeInAdvance_ + timeInFinalize_).str() << nl;
        Info << "  (S) initialize():              " << timeInInitialize_.str() << nl;
        Info << "  (I) advance():                 " << timeInAdvance_.str() << nl;
        Info << "  (I) finalize():                " << timeInFinalize_.str() << nl;
        Info << "  These times include time waiting for other participants." << nl;
        Info << "  See also precice-profiling on the website https://precice.org/tooling-performance-analysis.html." << nl;
        Info << "-------------------------------------------------------------------------------------" << nl;)

    std::cout << "[DEBUG] ~Adapter() destructor completed." << std::endl;
    return;
}
