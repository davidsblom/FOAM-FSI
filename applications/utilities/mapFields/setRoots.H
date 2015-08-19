    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("source dir");

    argList::validOptions.insert("sourceTime", "scalar");
    argList::validOptions.insert("parallelSource", "");
    argList::validOptions.insert("parallelTarget", "");
    argList::validOptions.insert("consistent", "");
    argList::validOptions.insert("copyphi", "");

    argList args(argc, argv);

    if (!args.check())
    {
         FatalError.exit();
    }

    fileName rootDirTarget(args.rootPath());
    fileName caseDirTarget(args.globalCaseName());

    fileName casePath(args.additionalArgs()[0]);
    fileName rootDirSource = casePath.path();
    fileName caseDirSource = casePath.name();

    Info<< "Source: " << rootDirSource << " " << caseDirSource << nl
        << "Target: " << rootDirTarget << " " << caseDirTarget << endl;

    bool parallelSource = args.optionFound("parallelSource");
    bool parallelTarget = args.optionFound("parallelTarget");
    bool consistent     = args.optionFound("consistent");
    bool copyphi        = args.optionFound("copyphi");
