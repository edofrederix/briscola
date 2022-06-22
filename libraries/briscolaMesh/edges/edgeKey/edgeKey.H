#ifndef edgeKey_H
#define edgeKey_H

#include "Pair.H"
#include "label.H"

namespace Foam
{

namespace briscola
{

class edgeKey;

bool operator==(const edgeKey&, const edgeKey&);
bool operator!=(const edgeKey&, const edgeKey&);

Istream& operator>>(Istream&, edgeKey&);
Ostream& operator<<(Ostream&, const edgeKey&);

class edgeKey
:
    public Pair<label>
{
public:

        class hash
        :
            public Hash<edgeKey>
        {
        public:

            // Construct null
            hash();

            // Generate a hash from a phase pair key
            label operator()(const edgeKey& key) const;
        };

public:

    // Construct null
    edgeKey();

    // Construct from vertex labels
    edgeKey(const label v1, const label v2);

    // Destructor
    virtual ~edgeKey();

    // Test if keys are equal
    friend bool operator==(const edgeKey& a, const edgeKey& b);

    // Test if keys are unequal
    friend bool operator!=(const edgeKey& a, const edgeKey& b);

    // Read from stdin
    friend Istream& operator>>(Istream& is, edgeKey& key);

    // Write to stdout
    friend Ostream& operator<<(Ostream& os, const edgeKey& key);
};

inline edgeKey reverse(const edgeKey& e)
{
    return edgeKey(e.second(), e.first());
}

}

}

#endif
