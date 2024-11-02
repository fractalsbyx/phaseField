#ifndef PHASE_H
#define PHASE_H


#include <map>
#include <string>

struct PhaseCompInfo{
    double _M;
    double M;
};

struct Phase {
public:
    std::string name;
    std::map<std::string, PhaseCompInfo> comps;
    double _sigma;
    double _mu;
    double sigma;
    double mu;
};

#endif