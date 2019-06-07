#ifndef CRONOMETRO_HPP_ 
#define CRONOMETRO_HPP_

#include <ctime>
using namespace std;

const double cps = CLOCKS_PER_SEC; //lo necesitamos para obtener el tiempo
                                   // segudnos ya que clock() lo devuelve en pulsos/segundo

class cronometro 
{
public:
    cronometro(){activo = false;}
    void activar() {activo = true; t0 = clock();}
    void parar() {if(activo){t1 = clock(); activo = false;}}
    double tiempo() const {return ((activo ? clock() : t1) - t0) / cps;}

private:
    bool activo;
    clock_t t0, t1;
};



#endif