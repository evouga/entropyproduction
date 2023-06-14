#include <deque>
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <random>
#include <limits>
#include <bitset>
#include <unordered_map>
#include <numeric>
#include <map>
#include <cassert>

struct Edge
{
    int state;
    double rate;
    bool emit;
};

struct Emission
{
    int state;
    double time;
};

void simulate(const std::vector<std::vector<Edge> >& edges, int nsteps, std::vector<Emission> &results)
{
    int nstates = edges.size();
    double time = 0;

    std::random_device rd;
    std::uniform_real_distribution<> dist(0.0, 1.0);
    std::mt19937 rng(rd());

    std::uniform_int_distribution<> statedist(0, nstates-1);
    int curstate = statedist(rng);

    while(results.size() < nsteps)
    {
        double totrate = 0;
        for (auto e : edges[curstate])
            totrate += e.rate;
        int nout = edges[curstate].size();
        if (nout == 0 || totrate == 0.0)
        {
            std::cerr << "Stuck in unescapable state: " << curstate << std::endl;
            exit(-1);
        }
        std::vector<double> prefixsums(nout + 1);
        double sum = 0;
        for (int i = 0; i < nout; i++)
        {
            prefixsums[i + 1] = prefixsums[i] + edges[curstate][i].rate;
        }
        double sample = dist(rng) * totrate;
        int next = 0;
        while ((next + 1 < nout) && (prefixsums[next + 1] < sample))
            next++;
        time -= 1.0 / totrate * std::log(dist(rng));
        if (edges[curstate][next].emit)
            results.push_back({ curstate, time });
        curstate = edges[curstate][next].state;
    }
}

void computeLagTimeBreakpoints(const std::vector<Emission>& emissions, int nTimeBuckets, std::vector<double>& breakpoints)
{
    breakpoints.resize(nTimeBuckets);
    int nemissions = emissions.size();
    std::vector<double> lagtimes(nemissions - 1);
    for (int i = 1; i < nemissions; i++)
    {
        lagtimes[i - 1] = emissions[i].time - emissions[i - 1].time;
    }
    std::sort(lagtimes.begin(), lagtimes.end());
    int nlagtimes = lagtimes.size();
    for (int i = 0; i < nTimeBuckets; i++)
    {
        int idx = ((nlagtimes - 1) * (i + 1)) / nTimeBuckets;
        double breakpoint = lagtimes[idx];
        breakpoints[i] = breakpoint;
    }
}

void discretizeStates(const std::vector<Emission>& emissions, std::map<int, int> &statemap, const std::vector<double> &lagTimeBreakpoints, std::vector<int>& discreteStates)
{
    int nemissions = emissions.size();
    int nTimeBuckets = lagTimeBreakpoints.size();
    discreteStates.resize(nemissions-1);

    for (int i = 1; i < nemissions; i++)
    {
        int stateidx = statemap[emissions[i].state];
        double lagtime = emissions[i].time - emissions[i - 1].time;
        int lagtimebucket = 0;
        while (lagtimebucket + 1 < nTimeBuckets && lagtime > lagTimeBreakpoints[lagtimebucket])
            lagtimebucket++;
        discreteStates[i - 1] = stateidx * nTimeBuckets + lagtimebucket;
    }
}

void reverse(const std::vector<Emission>& forward, std::vector<Emission>& backward)
{
    int nresults = forward.size();
    backward.resize(nresults);
    for (int i = 0; i < nresults; i++)
    {
        backward[i].state = forward[nresults - i - 1].state;
        backward[i].time = forward[nresults - 1].time - forward[nresults - i - 1].time;
    }
}

int main()
{
    int nstates, nedges;
    std::cin >> nstates >> nedges;
    std::vector<std::vector<Edge> > edges(nstates);
    std::map<int, int> emittermap;
    int nemitters = 0;
    for (int i = 0; i < nedges; i++)
    {
        Edge e;
        int src;
        std::cin >> src >> e.state >> e.rate >> e.emit;
        edges[src].push_back(e);
        if (e.emit)
        {
            emittermap[src] = nemitters++;
        }
    }

    int64_t nsteps;
    int nbuckets;
    std::cin >> nsteps >> nbuckets;
    std::vector<Emission> results;
    simulate(edges, nsteps, results);
    std::vector<int> discreteResults;
    std::vector<double> lagTimeBreakpoints;
    computeLagTimeBreakpoints(results, nbuckets, lagTimeBreakpoints);
    discretizeStates(results, emittermap, lagTimeBreakpoints, discreteResults);
    std::vector<Emission> reversedResults;
    reverse(results, reversedResults);
    std::vector<int> reversedDiscreteResults;
    discretizeStates(reversedResults, emittermap, lagTimeBreakpoints, reversedDiscreteResults);
    int nbins = nemitters * nbuckets;
    std::vector<int> fwdhisto(nbins * nbins * nbins);
    std::vector<int> backhisto(nbins * nbins * nbins);
    int nresults = discreteResults.size();
    for (int i = 2; i < nresults; i++)
    {
        fwdhisto[discreteResults[i - 2] * nbins * nbins + discreteResults[i - 1] * nbins + discreteResults[i]]++;
        backhisto[reversedDiscreteResults[i - 2] * nbins * nbins + reversedDiscreteResults[i - 1] * nbins + reversedDiscreteResults[i]]++;
    }
    for (int i = 0; i < nbins; i++)
    {
        for (int j = 0; j < nbins; j++)
        {
            for (int k = 0; k < nbins; k++)
            {
                int idx = i * nbins * nbins + j * nbins + k;
                std::cout << "(" << i << ", " << j << ", " << k << "): " << fwdhisto[idx] << " " << backhisto[idx] << std::endl;
            }
        }
    }
}