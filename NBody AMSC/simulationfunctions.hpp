//Committed by Elena Nuttini

// File adibito alla gestione delle funzioni usate globalmente per fare la simulazione. 
// possiamo già identificare un metodo utile per far proseguire la simulazione dove specifica ciò che accade per ogni timestep

#ifndef PARTICLE_UTIL_HPP
#define PARTICLE_UTIL_HPP

#include "Particle.hpp"
#include "Constants.hpp"
#include <vector>
#include <random>


class simulationfunctions{
public:
  
  
   // Metodo statico per generare una collezione di oggetti Particle
   static std::vector<Particle<dim>> generateParticles() {
       std::vector<Particle<dim>> particles;
       //PER DEBUGGING
       //particles.push_back(Particle<dim>(1, Arrows<dim>({1.0, 2.0, 3.0}), Arrows<dim>({0.0, 0.0, 0.0}), Arrows<dim>({0.0,0.0,0.0}), Arrows<dim>({0.0,0.0,0.0}), 5.0));
       //particles.push_back(Particle<dim>(2, Arrows<dim>({3.0, 2.0, 1.0}), Arrows<dim>({0.0, 0.0, 0.0}), Arrows<dim>({0.0,0.0,0.0}), Arrows<dim>({0.0,0.0,0.0}), 7.0));
      
       // Metodo che genera un numero prefissato di particelle di dimensione DIM con valori casuali
       std::random_device rd;
       std::mt19937 gen(rd());
       std::uniform_real_distribution<double> distribution(-1000.0, 1000.0); // Cambia il range se necessario
       std::uniform_real_distribution<double> massDistribution(0.1, 10.0); // Massa compresa tra 0.1 e 10.0 (valori arbitrari)


       for (int i = 0; i < 5; ++i) {
           // Genera valori casuali per posizione, velocità, massa, ecc.
           Arrows<dim> randomPosition;
           Arrows<dim> randomVelocity;
           //Arrows<dim> randomAcceleration;
           //Arrows<dim> randomCoefficients;


           for (unsigned int j = 0; j < dim; ++j) {
               randomPosition[j] = distribution(gen);
               randomVelocity[j] = distribution(gen);
               // Genera valori casuali per altre proprietà se necessario
           }


           double randomMass = massDistribution(gen);


           // Crea e aggiungi una particella con valori casuali alla collezione
           particles.push_back(Particle<dim>(i + 1, randomPosition, randomVelocity, Arrows<dim>({0.0, 0.0, 0.0}), Arrows<dim>({0.0, 0.0, 0.0}), randomMass));
       }




       return particles;
   }


   //Metodo che esegue la simulazione. Tale metodo provvederà a chiamare il generatore di Particles. PER ORA SOLVER BASE!!!
   static void doSim(){
      
       std::vector<Particle<dim>> particles = std::vector<Particle<dim>>();
       particles = generateParticles();
       for(Particle particle : particles){
           particle.setToZero();
           particle.print();
           std::cout<<" "<<std::endl;
       }
       stepSim(particles);
       updateSim(particles);
  
   }


   static void stepSim(std::vector<Particle<dim>>& particles){ 
       std::cout<<"!!!! INIZIO DELLA SIMULAZIONE !!!!"<<std::endl;
       std::cout<<" "<<std::endl;
      
       for(unsigned int i = 0; i<particles.size(); ++i){
           Arrows<dim> temp = Arrows<dim>();
           for(unsigned int j=0; j<particles.size(); ++j){
               if(i==j){continue;}
               temp += particles[i].calcCoefficients(particles[j]);
           }
           particles[i].coefficientsSetter(temp);
           particles[i].calcAccelleration();
       }
   }   
      
   static void updateSim(std::vector<Particle<dim>>& particles){
       double dt = 10;
       for(Particle particle : particles){
           particle.updateVelocity(dt);
           particle.updatePosition(dt);
           particle.print();
           std::cout<<" "<<std::endl;
       }
   }
  
 
   void solveCollisions(std::vector<Particle<dim>> particles){
      
   }


  
};


#endif // PARTICLE_UTIL_HPP


