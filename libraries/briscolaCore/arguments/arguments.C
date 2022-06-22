#include "arguments.H"
#include "fileOperation.H"
#include "fileOperationInitialise.H"

namespace Foam
{

namespace briscola
{

SLList<Foam::string>    arguments::validArgs;
HashTable<Foam::string> arguments::validOptions;
HashTable<Foam::string> arguments::optionUsage;
Foam::string::size_type arguments::usageMin = 20;
Foam::string::size_type arguments::usageMax = 80;

void arguments::addBoolOption
(
    const word& opt,
    const string& usage
)
{
    addOption(opt, "", usage);
}

void arguments::addOption
(
    const word& opt,
    const string& param,
    const string& usage
)
{
    validOptions.set(opt, param);
    if (!usage.empty())
    {
        optionUsage.set(opt, usage);
    }
}

void arguments::addUsage
(
    const word& opt,
    const string& usage
)
{
    if (usage.empty())
    {
        optionUsage.erase(opt);
    }
    else
    {
        optionUsage.set(opt, usage);
    }
}

void arguments::removeOption(const word& opt)
{
    validOptions.erase(opt);
    optionUsage.erase(opt);
}

void arguments::printOptionUsage
(
    const label location,
    const string& str
)
{
    const string::size_type textWidth = usageMax - usageMin;
    const string::size_type strLen = str.size();

    if (strLen)
    {
        // Minimum of 2 spaces between option and usage:
        if (string::size_type(location) + 2 <= usageMin)
        {
            for (string::size_type i = location; i < usageMin; ++i)
            {
                Info<<' ';
            }
        }
        else
        {
            // or start a new line
            Info<< nl;
            for (string::size_type i = 0; i < usageMin; ++i)
            {
                Info<<' ';
            }
        }

        // Text wrap
        string::size_type pos = 0;
        while (pos != string::npos && pos + textWidth < strLen)
        {
            // Potential end point and next point
            string::size_type curr = pos + textWidth - 1;
            string::size_type next = string::npos;

            if (isspace(str[curr]))
            {
                // We were lucky: ended on a space
                next = str.find_first_not_of(" \t\n", curr);
            }
            else if (isspace(str[curr+1]))
            {
                // The next one is a space - so we are okay
                curr++;  // otherwise the length is wrong
                next = str.find_first_not_of(" \t\n", curr);
            }
            else
            {
                // Search for end of a previous word break
                string::size_type prev = str.find_last_of(" \t\n", curr);

                // Reposition to the end of previous word if possible
                if (prev != string::npos && prev > pos)
                {
                    curr = prev;
                }
            }

            if (next == string::npos)
            {
                next = curr + 1;
            }

            // Indent following lines (not the first one)
            if (pos)
            {
                for (string::size_type i = 0; i < usageMin; ++i)
                {
                    Info<<' ';
                }
            }

            Info<< str.substr(pos, (curr - pos)).c_str() << nl;
            pos = next;
        }

        // Output the remainder of the string
        if (pos != string::npos)
        {
            // Indent following lines (not the first one)
            if (pos)
            {
                for (string::size_type i = 0; i < usageMin; ++i)
                {
                    Info<<' ';
                }
            }

            Info<< str.substr(pos).c_str() << nl;
        }
    }
    else
    {
        Info<< nl;
    }
}

// Convert argv -> args_
bool arguments::regroupArgv(int& argc, char**& argv)
{
    int nArgs = 0;
    int listDepth = 0;
    string tmpString;

    // Note: we also re-write directly into args_
    // and use a second pass to sort out args/options
    for (int argI = 0; argI < argc; ++argI)
    {
        if (strcmp(argv[argI], "(") == 0)
        {
            ++listDepth;
            tmpString += "(";
        }
        else if (strcmp(argv[argI], ")") == 0)
        {
            if (listDepth)
            {
                --listDepth;
                tmpString += ")";
                if (listDepth == 0)
                {
                    args_[nArgs++] = tmpString;
                    tmpString.clear();
                }
            }
            else
            {
                args_[nArgs++] = argv[argI];
            }
        }
        else if (listDepth)
        {
            // Quote each string element
            tmpString += "\"";
            tmpString += argv[argI];
            tmpString += "\"";
        }
        else
        {
            args_[nArgs++] = argv[argI];
        }
    }

    if (tmpString.size())
    {
        args_[nArgs++] = tmpString;
    }

    args_.setSize(nArgs);

    return nArgs < argc;
}

void arguments::getRootCase()
{
    fileName casePath;

    casePath = cwd();

    rootPath_ = casePath.path();
    case_ = casePath.name();

    if (rootPath_.isAbsolute())
    {
        setEnv("BRISCOLA_CASE", rootPath_/case_, true);
        setEnv("BRISCOLA_CASENAME", case_, true);
    }
    else
    {
        casePath = cwd()/rootPath_/case_;
        casePath.clean();

        setEnv("BRISCOLA_CASE", casePath, true);
        setEnv("BRISCOLA_CASENAME", casePath.name(), true);
    }
}

arguments::arguments
(
    int& argc,
    char**& argv,
    bool checkArgs,
    bool checkOpts
)
:
    args_(argc),
    options_(argc)
{
    word handlerType = fileOperation::defaultFileHandler;

    bool needsThread = fileOperations::fileOperationInitialise::New
    (
        handlerType,
        argc,
        argv
    )().needsThreading();

    for (int argI = 0; argI < argc; ++argI)
    {
        if (argv[argI][0] == '-')
        {
            const char *optionName = &argv[argI][1];

            if (word(optionName) == "parallel")
            {
                parRunControl_.runPar(argc, argv, needsThread);
                break;
            }
        }
    }

    regroupArgv(argc, argv);

    args_[0] = fileName(argv[0]);
    executable_ = fileName(argv[0]).name();

    int nArgs = 1;
    argumentsStr_ = args_[0];

    for (int argI = 1; argI < args_.size(); ++argI)
    {
        argumentsStr_ += ' ';
        argumentsStr_ += args_[argI];

        if (args_[argI][0] == '-')
        {
            const char *optionName = &args_[argI][1];

            if
            (
                validOptions.found(optionName)
             && !validOptions[optionName].empty()
            )
            {
                ++argI;
                if (argI >= args_.size())
                {
                    FatalError
                        <<"Option '-" << optionName
                        << "' requires an argument" << endl;
                    printUsage();
                    FatalError.exit();
                }

                argumentsStr_ += ' ';
                argumentsStr_ += args_[argI];
                options_.insert(optionName, args_[argI]);
            }
            else
            {
                options_.insert(optionName, "");
            }
        }
        else
        {
            if (nArgs != argI)
            {
                args_[nArgs] = args_[argI];
            }
            ++nArgs;
        }
    }

    args_.setSize(nArgs);

    parse(checkArgs, checkOpts);
}

void arguments::parse
(
    bool checkArgs,
    bool checkOpts
)
{
    if
    (
        options_.found("help")
    )
    {
        printUsage();
        ::exit(0);
    }

    if (!check(checkArgs, checkOpts))
    {
        FatalError.exit();
    }

    // Set default fileHandler

    autoPtr<fileOperation> handler
    (
        fileOperation::New
        (
            fileOperation::defaultFileHandler,
            false
        )
    );

    fileHandler(handler);

    getRootCase();

    if (parRunControl_.parRun())
    {
        if (Pstream::master())
        {
            for
            (
                int slave = Pstream::firstSlave();
                slave <= Pstream::lastSlave();
                slave++
            )
            {
                OPstream toSlave(Pstream::commsTypes::scheduled, slave);
                toSlave << args_ << options_;
            }
        }
        else
        {
            IPstream fromMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );
            fromMaster >> args_ >> options_;

            getRootCase();
        }
    }
}

arguments::~arguments()
{
    autoPtr<fileOperation> dummy(nullptr);
    fileHandler(dummy);
}

bool arguments::setOption(const word& opt, const string& param)
{
    bool changed = false;

    // Only allow valid options

    if (validOptions.found(opt))
    {
        // Protect the parallel option

        if (opt == "parallel")
        {
            FatalError
                <<"used arguments::setOption on a protected option: '"
                << opt << "'" << endl;
            FatalError.exit();
        }

        if (validOptions[opt].empty())
        {
            // Bool option
            if (!param.empty())
            {
                // Disallow change of type
                FatalError
                    <<"used arguments::setOption to change bool to non-bool: '"
                    << opt << "'" << endl;
                FatalError.exit();
            }
            else
            {
                // Did not previously exist
                changed = !options_.found(opt);
            }
        }
        else
        {
            // Non-bool option
            if (param.empty())
            {
                // Disallow change of type
                FatalError
                    <<"used arguments::setOption to change non-bool to bool: '"
                    << opt << "'" << endl;
                FatalError.exit();
            }
            else
            {
                // Existing value needs changing, or did not previously exist
                changed = options_.found(opt) ? options_[opt] != param : true;
            }
        }
    }
    else
    {
        FatalError
            <<"used arguments::setOption on an invalid option: '"
            << opt << "'" << nl << "allowed are the following:"
            << validOptions << endl;
        FatalError.exit();
    }

    // Set/change the option as required
    if (changed)
    {
        options_.set(opt, param);
    }

    return changed;
}


bool arguments::unsetOption(const word& opt)
{
    // Only allow valid options
    if (validOptions.found(opt))
    {
        // Protect the parallel option
        if (opt == "parallel")
        {
            FatalError
                <<"used arguments::unsetOption on a protected option: '"
                << opt << "'" << endl;
            FatalError.exit();
        }

        // Remove the option, return true if state changed
        return options_.erase(opt);
    }
    else
    {
        FatalError
            <<"used arguments::unsetOption on an invalid option: '"
            << opt << "'" << nl << "allowed are the following:"
            << validOptions << endl;
        FatalError.exit();
    }

    return false;
}

void arguments::printUsage() const
{
    Info<< "\nUsage: " << executable_ << " [OPTIONS]";

    forAllConstIter(SLList<string>, validArgs, iter)
    {
        Info<< " <" << iter().c_str() << '>';
    }

    Info<< "\noptions:\n";

    wordList opts = validOptions.sortedToc();
    forAll(opts, optI)
    {
        const word& optionName = opts[optI];

        HashTable<string>::const_iterator iter = validOptions.find(optionName);
        Info<< "  -" << optionName;
        label len = optionName.size() + 3;  // Length includes leading '  -'

        if (iter().size())
        {
            // Length includes space and between option/param and '<>'
            len += iter().size() + 3;
            Info<< " <" << iter().c_str() << '>';
        }

        HashTable<string>::const_iterator usageIter =
            optionUsage.find(optionName);

        if (usageIter != optionUsage.end())
        {
            printOptionUsage
            (
                len,
                usageIter()
            );
        }
        else
        {
            Info<< nl;
        }
    }

    Info<< "  -help";
    printOptionUsage
    (
        7,
        "print the usage"
    );
}

bool arguments::check(bool checkArgs, bool checkOpts) const
{
    bool ok = true;

    if (Pstream::master())
    {
        if (checkArgs && args_.size() - 1 != validArgs.size())
        {
            FatalError
                << "Wrong number of arguments, expected " << validArgs.size()
                << " found " << args_.size() - 1 << endl;
            ok = false;
        }

        if (checkOpts)
        {
            forAllConstIter(HashTable<string>, options_, iter)
            {
                if (!validOptions.found(iter.key()))
                {
                    FatalError
                        << "Invalid option: -" << iter.key() << endl;
                    ok = false;
                }
            }
        }

        if (!ok)
        {
            printUsage();
        }
    }

    return ok;
}

bool arguments::checkRootCase() const
{
    if (!fileHandler().isDir(rootPath()))
    {
        FatalError
            << executable_
            << ": cannot open root directory " << rootPath()
            << endl;

        return false;
    }

    fileName pathDir(fileHandler().filePath(path()));

    if (pathDir.empty() && Pstream::master())
    {
        // Allow slaves on non-existing processor directories, created later
        // (e.g. redistributePar)
        FatalError
            << executable_
            << ": cannot open case directory " << path()
            << endl;

        return false;
    }

    return true;
}

}

}
