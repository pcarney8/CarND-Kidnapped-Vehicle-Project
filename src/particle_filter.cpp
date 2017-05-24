/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"
using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	default_random_engine gen;
	num_particles = 10;
	cout << "init particles" << endl;
	// Add random Gaussian noise to each particle.
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	for(int i = 0; i < num_particles; i++){
		Particle particle;
		particle.id = i;
		particle.weight = 1;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);

		particles.push_back(particle);
	}
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    //for readability
    double std_x = std_pos[0];
    double std_y = std_pos[1];
    double std_theta = std_pos[2];
    cout << "prediction" << endl;
    //Zero mean Gaussians
    normal_distribution<double> dist_x(0, std_x);
    normal_distribution<double> dist_y(0, std_y);
    normal_distribution<double> dist_theta(0, std_theta);
	default_random_engine gen;

    for(int i = 0; i < num_particles; i++){
        double x_0 = particles[i].x;
        double y_0 = particles[i].y;
        double theta_0 = particles[i].theta;
        double theta_add = theta_0 + yaw_rate*delta_t;

        //state prediction
        particles[i].x = x_0 + (velocity/yaw_rate)*(sin(theta_add) - sin(theta_0));
        particles[i].y = y_0 + (velocity/yaw_rate)*(cos(theta_0) - cos(theta_add));
        particles[i].theta = theta_add;

        //Add noise
        particles[i].x += dist_x(gen);
        particles[i].y += dist_y(gen);
        particles[i].theta += dist_theta(gen);
    }


}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html
    double pi = atan(1)*4;
    cout << "update weights" << endl;
    for(auto &particle : particles){
        for(auto &observation : observations){
            //transform observations to map coordinates
            observation.x = particle.x + observation.x*cos(particle.theta) - observation.y*sin(particle.theta);
            observation.y = particle.y + observation.x*sin(particle.theta) + observation.y*cos(particle.theta);
            //todo: potentially do the dataAssociation

	   //todo: there's a bug somewhere here becuase weights are all 0
            //calculate the new particle weight and multiply it back into the particle.weight
            double std_x = std_landmark[0];
            double std_y = std_landmark[1];
            double exp_x = pow(observation.x - particle.x, 2)/(2. * std_x * std_x);
            double exp_y = pow(observation.y - particle.y,  2)/(2. * std_y * std_y);
            double w = exp(-1*(exp_x + exp_y))/(2. * pi * std_x * std_y);
	    cout << "calculated w: " << w << ", ";
            particle.weight = particle.weight * w;
        }
    }

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
//    default_random_engine gen;
//    uniform_real_distribution<> dis(0, 1);
    int index = rand() % num_particles;
    double beta = 0.0;
    cout << "resample " << endl;
    Particle max_weight_particle = *max_element(begin(particles), end(particles),
        [] (const Particle& p1, const Particle& p2) {
        return p1.weight < p2.weight;
    });

    double max_weight = max_weight_particle.weight;
    vector<Particle> new_particles;
    cout << "debug, max_weight: " << max_weight << endl;
//    cout << "debug, dis(gen): " << dis(gen) << endl;

    for(int i = 0; i < num_particles; i++){
//        beta += dis(gen) * 2.0 * max_weight;
        beta += rand() * 2.0 * max_weight;

        //wheel
        while(beta > particles[index].weight){
            beta -= particles[index].weight;
            index = (index + 1) % num_particles;
        }

        //add to new particles vector
        new_particles.push_back(particles[index]);
    }
    cout << "finished resample wheel" << endl;
    //assign back to particles with the new particles
    particles = new_particles;
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
