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
	random_device device;
	mt19937 gen(device());
	num_particles = 2;
//	cout << "init particles" << endl;
	// Add random Gaussian noise to each particle.
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	for(int i = 0; i < num_particles; i++){
		Particle particle;
		particle.id = i;
		particle.weight = 1.0;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);

		particles.push_back(particle);
	}
	cout << "particle: " << particles[0].x << ", " << particles[0].y << ", " << particles[0].theta << endl;
	is_initialized = true;
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
  //  cout << "prediction" << endl;
    //Zero mean Gaussians
    random_device device;
    mt19937 gen(device());
    normal_distribution<double> dist_x(0.0, std_x);
    normal_distribution<double> dist_y(0.0, std_y);
    normal_distribution<double> dist_theta(0.0, std_theta);

    for(int i = 0; i < num_particles; i++){
        double x_0 = particles[i].x;
        double y_0 = particles[i].y;
        double theta_0 = particles[i].theta;
        double theta_add = theta_0 + yaw_rate*delta_t;

	if (fabs(yaw_rate) == 0){
		particles[i].x = x_0 + velocity*delta_t*cos(theta_0);
	        particles[i].y = y_0 + velocity*delta_t*sin(theta_0);
       		particles[i].theta = theta_0;
	} else {
	        particles[i].x = x_0 + (velocity/yaw_rate)*(sin(theta_add) - sin(theta_0));
	        particles[i].y = y_0 + (velocity/yaw_rate)*(cos(theta_0) - cos(theta_add));
       		particles[i].theta = theta_add;
	}
        //Add noise
        particles[i].x += dist_x(gen);
        particles[i].y += dist_y(gen);
        particles[i].theta += dist_theta(gen);
    }

	cout << "particle: " << particles[0].x << ", " << particles[0].y << ", " << particles[0].theta << endl;
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
   // cout << "update weights" << endl;
   // cout << "particle weight" << particles[0].weight << endl;
    for(auto& particle : particles){
	particle.weight = 1.0;
        for(auto& observation : observations){
            //transform observations to map coordinates
            double new_observation_x = particle.x + observation.x*cos(particle.theta) - observation.y*sin(particle.theta);
            double new_observation_y = particle.y + observation.x*sin(particle.theta) + observation.y*cos(particle.theta);
            
            //todo: potentially do the dataAssociation

	   //todo: there's a bug somewhere here becuase weights are all 0
            //calculate the new particle weight and multiply it back into the particle.weight
            double std_x = std_landmark[0];
            double std_y = std_landmark[1];
            double exp_x = pow(new_observation_x - particle.x, 2.0)/(2.0 * std_x * std_x);
            double exp_y = pow(new_observation_y - particle.y,  2.0)/(2.0 * std_y * std_y);
            double w = (1.0/(2.0 * M_PI * std_x * std_y))*exp(-1.0*(exp_x + exp_y));
	    cout << "calculated w: " << w << endl;
	    cout << "particle w: " << particle.weight << endl;
            particle.weight = particle.weight * w;
        }
	weights.push_back(particle.weight);
    }

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    random_device device;
    mt19937 gen(device());
    discrete_distribution<> dist(weights.begin(), weights.end());

    vector<Particle> new_particles;
   // cout << "debug, max_weight: " << max_weight << endl;

    for(int i = 0; i < num_particles; i++){
        //add to new particles vector
        new_particles.push_back(particles[dist(gen)]);
    }
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
