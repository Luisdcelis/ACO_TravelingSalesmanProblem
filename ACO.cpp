#include <iostream>
#include <vector>
#include <fstream>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <utility>

#include "cronometro.hpp"

using namespace std;

typedef vector<vector<double> > graph;
typedef vector<int> tour;


double d_euclid(pair<double, double> A, pair<double, double> B)
{
    return sqrt((A.first - B.first) * (A.first - B.first) + (A.second - B.second) * (A.second - B.second));
}

graph getGraph(string file)
{
    ifstream f;
    f.open(file, ios::in);

    vector<pair<double, double> > coord;

    string s;
    for (size_t i = 0; i < 6; i++)
        getline(f, s);
    double cx, cy;
    while (s != "EOF")
    {
        f >> s;
        f >> cx;
        f >> cy;
        cout << "(" << cx << "," << cy << ")" << endl;
        coord.push_back(make_pair(cx, cy));
    }
    coord.pop_back();
    f.close();
    
    graph res;
    for(int i = 0; i < coord.size(); i++)
    {
        vector<double> aux;
        double dis;
        for(int j = 0; j < coord.size(); j++)
        {
            if(i == j)
                dis = RAND_MAX;
            else
                dis = d_euclid(coord[i], coord[j]);
            
            aux.push_back(dis);
        }
        res.push_back(aux);
    }
    return res;

}

tour getTour(string file)
{
    ifstream f;
    f.open(file, ios::in);

    tour path;

    string s;
    for (size_t i = 0; i < 4; i++)
        getline(f, s);
    int c;
    while (c != -1)
    {

        f >> c;
        path.push_back(c);
    }
    path.pop_back();            
    f.close();
    
    return path;

}

double costTravelSolution(graph cities, tour path)
{
    double cost = 0;
    for(int i = 0; i < path.size()-1; i++)
    {
        cost += cities[path[i]-1][path[i+1]-1];
    }
    cost += cities[path[path.size()-1]-1][path[0]-1];
    return cost;
}

double costTravel(graph cities, tour path)
{
    double cost = 0;
    for(int i = 0; i < path.size()-1; i++)
    {
        cost += cities[path[i]][path[i+1]];
    }
    cost += cities[path[path.size()-1]][path[0]];
    return cost;
}

// El camino de la solucion no vuelve a la primera ciudad
// pero esto se toma en cuenta en la funcion coste
tour CreateHeuristicSolution(int nCities)
{
    tour path;

    for(int i = 0; i < nCities; i++)
    {
        path.push_back(i);
    }

    return path;
}

graph initializePheromones(int nCities, double Pheromone_i)
{
    graph Pheromone;
    
    for(int i = 0; i < nCities; i++)
    {
        vector<double> aux;
        for(int j = 0; j < nCities; j++)
        {
            aux.push_back(Pheromone_i);
        }
        Pheromone.push_back(aux);
    }
    return Pheromone;
}

graph calculateEta(graph cities)
{
    for(int i = 0; i < cities.size(); i++)
    {
        for(int j = 0; j < cities.size(); j++)
        {
            cities[i][j] = 1.0/cities[i][j];
        }
    }
    return cities;
}

vector<double> calculateProb(graph Pheromones, graph eta, double alpha,
                            double beta, vector<int> UnvisitedNodes, int actNode)
{
    vector<double> probs;
    double denominator = 0;

    for(int i = 0; i < UnvisitedNodes.size(); i++)
    {
        denominator += pow(Pheromones[actNode][UnvisitedNodes[i]], alpha)*pow(eta[actNode][UnvisitedNodes[i]], beta);
    }

    double numerator;
    for(int i = 0; i < UnvisitedNodes.size(); i++)
    {
        numerator = pow(Pheromones[actNode][UnvisitedNodes[i]], alpha)*pow(eta[actNode][UnvisitedNodes[i]], beta);
        probs.push_back(numerator/denominator);
    }

    return probs;
     
}

tour ConstructSolution(graph cities, graph Pheromones,double alpha,
                       double beta, double q0)
{
    tour res;
    res.push_back(0);

    graph eta = calculateEta(cities);
    vector<int> UnvisitedNodes;

    for(int i = 1; i < cities.size(); i++)
    {
        UnvisitedNodes.push_back(i);
    }

    int actNode = 0;
    while(!UnvisitedNodes.empty())
    {
        
        vector<double> probs = calculateProb(Pheromones, eta, alpha, beta, UnvisitedNodes, actNode);
        
        
        bool selected = false;
        int epoch = 0;
        while(!selected)
        {
        double ra = ((double)rand() / (RAND_MAX + 1.0));
            if(ra <= q0)
            {
                double minimo = RAND_MAX;
                int pos;
                for(int i = 0; i <  UnvisitedNodes.size(); i++)
                {
                    if(cities[actNode][UnvisitedNodes[i]] < minimo)
                    {
                        minimo = cities[actNode][UnvisitedNodes[i]];
                        pos = UnvisitedNodes[i];
                    }
                }
                
                actNode =  pos;
                
                res.push_back(actNode);
                auto it = find(UnvisitedNodes.begin(), UnvisitedNodes.end(), actNode);
            
                UnvisitedNodes.erase(it);
                selected = true;
            }else{
                while(!selected && epoch < probs.size())
                {
                    double r = ((double)rand() / (RAND_MAX + 1.0));

                    if(r <= probs[epoch])
                    {
                        auto it = UnvisitedNodes.begin() + epoch;
                        actNode = UnvisitedNodes[epoch];
                        res.push_back(actNode);
                        UnvisitedNodes.erase(it);
                        selected = true;
                    }
                    epoch++;
                }
            }
        } 
    }
    return res;
}

void LocalUpdateAndDecayPheromone(graph& Pheromones, tour Si, graph Pheromones_0, double sigma)
{
    for(int i = 0; i < Si.size()-1; i++)
    {
        Pheromones[Si[i]][Si[i+1]] =(1-sigma)*Pheromones[Si[i]][Si[i+1]] + sigma * Pheromones_0[Si[i]][Si[i+1]];
    }
}

void GlobalUpdateAndDecayPheromone(graph& Pheromones, graph cities, tour pBest, double decay_factor)
{
    //primero solo  evaporamos
    for(int i = 0; i < Pheromones.size(); i++)
    {
        for(int j = 0; j < Pheromones.size(); j++)
        {
            Pheromones[i][j] = (1-decay_factor)*Pheromones[i][j];
        }
    }

    for(int i = 0; i < pBest.size()-1; i++)
    {
        Pheromones[pBest[i]][pBest[i+1]] +=  decay_factor * 1/(cities[pBest[i]][pBest[i+1]]);
    }
}




tour AntColonyOptimization(graph cities, int max_it, int max_rep, 
                           int m, double decay_factor, double alpha,
                           double beta, double q0, double sigma)
{
    tour pBest = CreateHeuristicSolution(cities.size());
    double pBestCost = costTravel(cities, pBest);
    double Pheromone_i = 1.0/(cities.size()*pBestCost);
    graph Pheromones = initializePheromones(cities.size(), Pheromone_i);
    graph Pheromones_0 = Pheromones;

    double pBestCost_ant;
    int cont = 0;
    
    int StopCondition = 0;
    while(StopCondition < max_it && cont < max_rep)
    {
        pBestCost_ant = pBestCost;
        for(int i = 0; i < m; i++)
        {
            tour Si = ConstructSolution(cities, Pheromones, alpha, beta, q0);
            double SiCost = costTravel(cities, Si);

            if(SiCost < pBestCost)
            {
                pBestCost = SiCost;
                pBest = Si;
                cout << "he actualizado pBest con: " << pBestCost << endl;
            }
            LocalUpdateAndDecayPheromone(Pheromones, Si, Pheromones_0, sigma);
            
            cout << "Coste hormiga actual: " << SiCost << endl;
        }
        GlobalUpdateAndDecayPheromone(Pheromones, cities, pBest, decay_factor);
        StopCondition++;

        if(pBestCost == pBestCost_ant)
            cont++;
        else   
            cont = 0;
    }
    return pBest;
}


tour AntColonyOptimization_Elitist(graph cities, int max_it, int max_rep, 
                           int m, double decay_factor, double alpha,
                           double beta, double q0, double sigma)
{
    tour pBest = CreateHeuristicSolution(cities.size());
    double pBestCost = costTravel(cities, pBest);
    double Pheromone_i = 1.0/(cities.size()*pBestCost);
    graph Pheromones = initializePheromones(cities.size(), Pheromone_i);
    graph Pheromones_0 = Pheromones;

    double pBestCost_ant;
    int cont = 0;
    
    
    int StopCondition = 0;
    while(StopCondition < max_it && cont < max_rep)
    {
        pBestCost_ant = pBestCost;
        for(int i = 0; i < m; i++)
        {
            tour Si = ConstructSolution(cities, Pheromones, alpha, beta, q0);
            double SiCost = costTravel(cities, Si);

            if(SiCost < pBestCost)
            {
                pBestCost = SiCost;
                pBest = Si;
                cout << "he actualizado pBest con: " << pBestCost << endl;
            }
            LocalUpdateAndDecayPheromone(Pheromones, Si, Pheromones_0, sigma);

            //Elitist Update
            GlobalUpdateAndDecayPheromone(Pheromones, cities, pBest, decay_factor);
            
            //cout << "Coste Si: " << SiCost << endl;
        }
        GlobalUpdateAndDecayPheromone(Pheromones, cities, pBest, decay_factor);
        StopCondition++;
        
        if(pBestCost == pBestCost_ant)
            cont++;
        else
            cont = 0;
    }
    return pBest;
}




int main(int argc, char const *argv[])
{   
    time_t t;
    srand(time(&t));
    graph cities;
    tour path, res;
    cities = getGraph("TSPLIB/berlin52.tsp");

    cronometro c;
    int max_it = 100;
    int max_rep = 5;
    int m = 30;
    c.activar();
    res = AntColonyOptimization(cities, max_it, max_rep, m, 0.1, 1, 2.5, 0.9, 0.1);
    c.parar();


    cout << endl;
    double mio = costTravel(cities, res);
    cout << "El algoritmo ha tardado " << c.tiempo() << " segundos " << "para " << cities.size() << " ciudades " << endl;

    //Cost of optimal solution:
    path = getTour("TSPLIB/berlin52.opt.tour");
    double opt = costTravelSolution(cities, path)
    cout << "mejor encontrado:" << mio << endl;
    cout << "optimo: " << opt << endl;
    cout << "diferencia: " << mio-opt << endl;

    return 0;
}


