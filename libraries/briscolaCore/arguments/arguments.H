// Lightweight version of OpenFOAM's argList class. Processor-dependent case
// paths were removed and some parts were stripped

#ifndef arguments_H
#define arguments_H

#include "parRun.H"
#include "IStringStream.H"
#include "OSspecific.H"

#include "sigFpe.H"

namespace Foam
{

namespace briscola
{

class arguments
{
    // Parallel mode
    ParRunControl parRunControl_;

    // List of arguments
    stringList args_;

    // Option table
    HashTable<string> options_;

    // Executable name
    word executable_;

    // Root path of the case
    fileName rootPath_;

    // Name of the case
    fileName case_;

    // String of arguments
    string argumentsStr_;

    // Floating point exception trapping
    sigFpe sigFpe_;

    // Helper function for printUsage
    static void printOptionUsage
    (
        const label location,
        const string& str
    );

    // Get root path as cwd
    void getRootCase();

    // Transcribe argv into internal args_
    bool regroupArgv(int& argc, char**& argv);

public:

    // A list of valid (mandatory) arguments
    static SLList<string> validArgs;

    // A list of valid options
    static HashTable<string> validOptions;

    // A list of valid parallel options
    static HashTable<string> validParOptions;

    // Short usage information for validOptions
    static HashTable<string> optionUsage;

    // Min offset for displaying usage (default: 20)
    static string::size_type usageMin;

    // Max screen width for displaying usage (default: 80)
    static string::size_type usageMax;

    // Construct from argc and argv
    arguments
    (
        int& argc,
        char**& argv,
        bool checkArgs = true,
        bool checkOpts = true
    );

    // Destructor
    virtual ~arguments();

    // Parse
    void parse
    (
        bool checkArgs,
        bool checkOpts
    );

    // Name of executable without the path
    inline const word& executable() const
    {
        return executable_;
    }

    // Return root path
    inline const fileName& rootPath() const
    {
        return rootPath_;
    }

    // Return case name
    inline const fileName& caseName() const
    {
        return case_;
    }

    // Return parRunControl
    inline const ParRunControl& parRunControl() const
    {
        return parRunControl_;
    }

    // Return the path to the caseName
    inline fileName path() const
    {
        return rootPath()/caseName();
    }

    // Return arguments
    inline const stringList& args() const
    {
        return args_;
    }

    // Return non-const access to arguments
    inline stringList& args()
    {
        return args_;
    }

    // Return the argument corresponding to index.
    inline const string& arg(const label index) const
    {
        return args_[index];
    }

    // Return the number of arguments
    inline label size() const
    {
        return args_.size();
    }

    // Read a value from the argument at index.
    template<class T>
    inline T argRead(const label index) const;

    // Return options
    inline const Foam::HashTable<string>& options() const
    {
        return options_;
    }

    // Return non-const access to options
    inline Foam::HashTable<string>& options()
    {
        return options_;
    }

    // Return the argument string associated with the named option
    inline const string& option(const word& opt) const
    {
        return options_[opt];
    }

    // Return true if the named option is found
    inline bool optionFound(const word& opt) const
    {
        return options_.found(opt);
    }

    // Return an IStringStream from the named option
    inline IStringStream optionLookup(const word& opt) const
    {
        return IStringStream(options_[opt]);
    }

    // Read a value from the named option
    template<class T>
    inline T optionRead(const word& opt) const;

    // Read a value from the named option if present
    template<class T>
    inline bool optionReadIfPresent(const word& opt, T&) const;

    // Read a value from the named option if present
    template<class T>
    inline bool optionReadIfPresent
    (
        const word& opt,
        T&,
        const T& deflt
    ) const;

    // Read a value from the named option if present
    template<class T>
    inline T optionLookupOrDefault
    (
        const word& opt,
        const T& deflt
    ) const;

    // Read a List of values from the named option
    template<class T>
    List<T> optionReadList(const word& opt) const
    {
        return readList<T>(optionLookup(opt)());
    }

    // Return the argument corresponding to index
    inline const string& operator[](const label index) const;

    // Return the argument string associated with the named option
    //  \sa option()
    inline const string& operator[](const word& opt) const;

    // Add to a bool option to validOptions with usage information
    static void addBoolOption
    (
        const word& opt,
        const string& usage = ""
    );

    // Add to an option to validOptions with usage information
    static void addOption
    (
        const word& opt,
        const string& param = "",
        const string& usage = ""
    );

    // Add option usage information to optionUsage
    static void addUsage
    (
        const word& opt,
        const string& usage
    );

    // Remove option from validOptions and from optionUsage
    static void removeOption(const word& opt);

    // Set option directly
    bool setOption(const word& opt, const string& param = "");

    // Unset option directly
    bool unsetOption(const word& opt);

    // Print usage
    void printUsage() const;

    // Check argument list
    bool check(bool checkArgs=true, bool checkOpts=true) const;

    // Check root path and case path
    bool checkRootCase() const;
};

}

}

#include "argumentsI.H"

#endif
