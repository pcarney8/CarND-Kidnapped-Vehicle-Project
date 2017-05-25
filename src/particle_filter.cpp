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
	num_particles = 20;
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
//	cout << "particle: " << particles[0].x << ", " << particles[0].y << ", " << particles[0].theta << endl;
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
//    cout << "prediction" << endl;
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

//	cout << "particle: " << particles[0].x << ", " << particles[0].y << ", " << particles[0].theta << endl;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
    LandmarkObs temp_landmark;
    LandmarkObs transformed_observation;
    vector<LandmarkObs> transformed_observations;
    vector<LandmarkObs> landmarks_in_range;
    double distance;
    double min_dist;

    weights.clear();
//	cout << "update weights" << endl;
    
    for(auto& particle : particles){
	particle.weight = 1.0;
        landmarks_in_range.clear();
        //filter out the closest landmarks
//	cout << "about to open map" << endl;
        for (auto landmark : map_landmarks.landmark_list) {
            distance = dist(landmark.x_f, landmark.y_f, particle.x, particle.y);
            if (distance <= sensor_range) {
                temp_landmark.x = landmark.x_f;
                temp_landmark.y = landmark.y_f;
                temp_landmark.id = landmark.id_i;
                landmarks_in_range.push_back(temp_landmark);
            }
        }

        transformed_observations.clear();
        //transform observations
//	cout << "about to transform obs" << endl;
        for (size_t i=0; i<observations.size(); ++i) {
            transformed_observation.x = particle.x + (observations[i].x * cos(particle.theta)) - (observations[i].y * sin(particle.theta));
            transformed_observation.y = particle.y + (observations[i].x * sin(particle.theta)) + (observations[i].y * cos(particle.theta));
            transformed_observations.push_back(transformed_observation);
        }

        distance = 0.0;
        //DataAssociation find closest landmark
//	cout << "about to data associate" << endl; 
        for (auto& o : transformed_observations) {
            min_dist = 999999;

            for (size_t j=0; j<landmarks_in_range.size(); ++j) {
                distance = dist(o.x, o.y, landmarks_in_range[j].x, landmarks_in_range[j].y);
                if (distance < min_dist) {
                    min_dist = distance;
                    o.id = j;
                }
            }
        }

        //calculate the weight
        for(auto& observation : transformed_observations){
            double mu_x = landmarks_in_range[observation.id].x;
            double mu_y = landmarks_in_range[observation.id].y;

            double std_x = std_landmark[0];
            double std_y = std_landmark[1];
            double exp_x = pow(observation.x - mu_x, 2.0)/(2.0 * std_x * std_x);
            double exp_y = pow(observation.y - mu_y,  2.0)/(2.0 * std_y * std_y);
            double w = (1.0/(2.0 * M_PI * std_x * std_y))*exp(-1.0*(exp_x + exp_y));
	    //cout << "calculated w: " << w << endl;
	    //cout << "particle w: " << particle.weight << endl;
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
