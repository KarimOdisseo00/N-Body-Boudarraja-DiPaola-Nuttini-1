#ifndef PARTICLE_UTIL_HPP_PARALLEL
#define PARTICLE_UTIL_HPP_PARALLEL

#include "Particle.hpp"
#include "Constants.hpp"
#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>


class simulationFunctionsParallel{
public:
    
    
    // Metodo statico per generare una collezione di oggetti Particle
    static std::vector<Particle<dim>> generateParticles() {
        std::vector<Particle<dim>> particles;
        particles.reserve(numberOfParticles);
        //double M=1.989e30;
        //PER DEBUGGING
      /* particles.push_back(Particle<dim>(1, Arrows<dim>({0.0, 0.0, 0.0}), Arrows<dim>({0.0, 0.0, 0.0}), Arrows<dim>({0.0,0.0,0.0}), Arrows<dim>({0.0,0.0,0.0}), Arrows<dim>({0.0,0.0,0.0}),M)); // Sun

        double r = 5.79e10; // distance from the sun
double v = sqrt(G * M / r); // orbital velocity
particles.push_back(Particle<dim>(2, Arrows<dim>({r, 0.0, 0.0}), Arrows<dim>({0.0, v, 0.0}), Arrows<dim>({0.0,0.0,0.0}), Arrows<dim>({0.0,0.0,0.0}), Arrows<dim>({0.0,0.0,0.0}),0.0553));

// Venus
r = 1.082e11;
v = sqrt(G * M / r);
particles.push_back(Particle<dim>(3, Arrows<dim>({r, 0.0, 0.0}), Arrows<dim>({0.0, v, 0.0}), Arrows<dim>({0.0,0.0,0.0}), Arrows<dim>({0.0,0.0,0.0}), Arrows<dim>({0.0,0.0,0.0}),0.815));

// Earth
r = 1.496e11;
v = sqrt(G * M / r);
particles.push_back(Particle<dim>(4, Arrows<dim>({r, 0.0, 0.0}), Arrows<dim>({0.0, v, 0.0}), Arrows<dim>({0.0,0.0,0.0}), Arrows<dim>({0.0,0.0,0.0}), Arrows<dim>({0.0,0.0,0.0}),1.0));
        //particles.push_back(Particle<dim>(3, Arrows<dim>({4.0, 5.0, 6.0}), Arrows<dim>({0.0, 0.0, 0.0}), Arrows<dim>({0.0,0.0,0.0}), Arrows<dim>({0.0,0.0,0.0}), 10.0));
        //particles.push_back(Particle<dim>(4, Arrows<dim>({6.0, 5.0, 4.0}), Arrows<dim>({0.0, 0.0, 0.0}), Arrows<dim>({0.0,0.0,0.0}), Arrows<dim>({0.0,0.0,0.0}), 1.0));*/
        // Metodo che genera un numero prefissato di particelle di dimensione DIM con valori casuali*/
        // Sun
double M = 1.989e30; // mass of the sun
particles.push_back(Particle<dim>(1, Arrows<dim>({0.0, 0.0, 0.0}), Arrows<dim>({0.0, 0.0, 0.0}), Arrows<dim>({0.0,0.0,0.0}), Arrows<dim>({0.0,0.0,0.0}), Arrows<dim>({0.0,0.0,0.0}), M));

// Earth
double r = 1.496e11; // distance from the sun
double v = sqrt(G * M / r); // orbital velocity
particles.push_back(Particle<dim>(2, Arrows<dim>({r, 0.0, 0.0}), Arrows<dim>({0.0, v, 0.0}), Arrows<dim>({0.0,0.0,0.0}), Arrows<dim>({0.0,0.0,0.0}), Arrows<dim>({0.0,0.0,0.0}), 5.972e24));
       /*double v = sqrt(G * M / 1.496e11);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> distribution(-1000, 1000); // Cambia il range se necessario
        std::uniform_real_distribution<double> distributionVel(-10, 10);
        std::uniform_real_distribution<double> massDistribution(10.0, 10000); // Massa compresa tra 0.1 e 10.0 (valori arbitrari)
        
        for (unsigned int i = 0; i < numberOfParticles; ++i) {
            // Genera valori casuali per posizione, velocità, massa, ecc.
            Arrows<dim> randomPosition;
            Arrows<dim> randomVelocity;
            //Arrows<dim> randomAcceleration;
            //Arrows<dim> randomCoefficients;
            
            for (unsigned int j = 0; j < dim; ++j) {
                randomPosition[j] = distribution(gen);
                randomVelocity[j] = distributionVel(gen);
                // Genera valori casuali per altre proprietà se necessario
            }

            double randomMass = massDistribution(gen);

            // Crea e aggiungi una particella con valori casuali alla collezione
            particles.push_back(Particle<dim>(i + 1, randomPosition, randomVelocity, Arrows<dim>(), Arrows<dim>(), Arrows<dim>(),randomMass));


    }*/
        


        return particles;
    }
   
     //Salvataggio delle posizioni su una vector. utile per generare il file txt che dovrà essere letto per l'animazione grafica
    static void memPos(const std::vector<Particle<dim>> particles, std::vector<double>& positions){
        for(Particle particle : particles){
            for(unsigned int i=0; i<dim; ++i){
                positions.emplace_back(particle.getPositionCoordinate(i));
            }
        }
    }

     //INSERIMENTO BLOCCHI PARALLELI:
     static void doParallelSim(){
        
        double local_totalTime = totalTime;
        std::vector<double> positions = std::vector<double>();
        positions.reserve(numberOfParticles);

        
        std::vector<Particle<dim>> particles = std::vector<Particle<dim>>();
        particles.reserve(numberOfParticles);
        particles = generateParticles();
        
        for(Particle particle : particles){
            particle.setToZero();
            //particle.print();
            //std::cout<<" "<<std::endl;
        }
        
      
       int local_cycles = cycles;
      
       for(unsigned int i=0 ; i<cycles;i++){
           //std::cout<<"!!!! INIZIO DELLA SIMULAZIONE !!!!"<<std::endl;
           //std::cout<<" "<<std::endl;
           //std::cout<<"CICLO "<< i << ":" << std::endl;
           //std::cout<<" "<<std::endl;
           stepParallelSim(particles,positions);
           //updateSim(particles);      //sezione computata in quanto inserita nel processo di parallelizzazione 
                                        //stepParallelSim
       }
        writePositionToTXT(positions);

       
   }
   
static void stepParallelSim(std::vector<Particle<dim>>& particles, std::vector<double>& positions)
{
    double local_dt = dt;
    double half_dt = dt / 2;
    memPos(particles, positions); // Fill the positions at each cycle

    std::vector<Arrows<dim>> tempCoefficients(particles.size());

    std::vector<bool> collisions(particles.size(), false);

    #pragma omp parallel num_threads(8) shared(particles, local_dt, half_dt, tempCoefficients, collisions)
    {
        #pragma omp for schedule(static)
        for (unsigned int i = 0; i < particles.size(); ++i) {
            Arrows<dim> temp1 = Arrows<dim>();

            for (unsigned int j = 0; j < particles.size(); ++j) {
                if (i == j) { continue; }
                temp1 += particles[i].calcCoefficients(particles[j]);
                if (particles[i].collision(particles[j])) { collisions[i] = true; }
            }

            tempCoefficients[i] = temp1;
        }

        #pragma omp for schedule(static)
        for (unsigned int i = 0; i < particles.size(); ++i) {
            particles[i].coefficientsSetter(tempCoefficients[i]);
            if (collisions[i]) {
                particles[i].calcForcesAfterCollision();
            }
            else {
                particles[i].calcForces();
            }
            particles[i].calcAccelleration();
        }

        #pragma omp for schedule(static)
        for (unsigned int i = 0; i < particles.size(); ++i) {
            particles[i].updatePosition(local_dt);
        }

        #pragma omp for schedule(static)
        for (unsigned int i = 0; i < particles.size(); ++i) {
            particles[i].updateVelocity(half_dt);
        }

        #pragma omp for schedule(static)
        for (unsigned int i = 0; i < particles.size(); ++i) {
            particles[i].setToZero();
        }

        #pragma omp for schedule(static)
        for (unsigned int i = 0; i < particles.size(); ++i) {
            particles[i].updateVelocity(half_dt);
        }
    }
}
   

    
    //Determinazione del numero di cicli. Deve essere fornito un lasso di tempo totale e un dt
    unsigned int numberOfCycles(double totalTime, double dt){
        double temp = totalTime/dt;
        return static_cast<unsigned int>(std::floor(temp));
    }
    
    //Scrivi le coordinate delle particelle su file TXT. Necessario per implementare la graficazione della simulazione
    static void writePositionToTXT(const std::vector<double> positions) {
        std::ofstream outputFile("positions.txt");
        if (!outputFile.is_open()) {
            std::cerr << "Impossibile aprire il file!" << std::endl;
            return;
        }

        int coordinates = dim;
        int rows = positions.size() / coordinates;
        int dif = 0;
        int cycle = 0;
        int counter = cycle + 1;

        for (int i = 0; i < rows; ++i) {
            outputFile << cycle << ","; // Scrivi il numero del ciclo
            if(i >= numberOfParticles)
            {
                outputFile << i - dif << ","; //Scrivi l'ID della particella dopo il primo ciclo
            } 
            else 
            {
                outputFile << i << ","; // Scrivi l'ID della particella per il primo ciclo
            }

            for (int j = 0; j < coordinates; ++j) {
                outputFile << positions[i * coordinates + j]; // Scrivi le coordinate
                if (j != coordinates - 1) {
                    outputFile << ",";
                }
            }
            outputFile << std::endl; // Vai a capo dopo ogni riga
            

            if(i == numberOfParticles*counter - 1)
            {
                dif = counter*numberOfParticles;
                ++counter;
                ++cycle;
            }
            
        }

        outputFile.close();
    }

};

#endif // PARTICLE_UTIL_HPP
