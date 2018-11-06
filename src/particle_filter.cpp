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
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"
#include "map.h"
using namespace std;

#define NUM_PARTICLE 100

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).    
    this->num_particles = NUM_PARTICLE;
    default_random_engine gen;
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    // initialize all weights
    for (unsigned int i = 0; i < NUM_PARTICLE; i++) { 
        weights.push_back(1.0);
        Particle p; 
        p.id = i;
        p.x = dist_x(gen);
        p.y = dist_y(gen);
        p.theta = dist_theta(gen);
        particles.push_back(p);
    }
    this->is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    default_random_engine gen;
    for (int i = 0; i < this->num_particles; i++) {
        double x_updated = particles[i].x;
        double y_updated = particles[i].y;
        double theta_updated = particles[i].theta;

        if (yaw_rate == 0) {
           // yaw rate is zero
           x_updated += velocity * delta_t * cos(particles[i].theta);
           y_updated += velocity * delta_t * sin(particles[i].theta);
        }
        else {
           x_updated += velocity / yaw_rate * 
                            (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
  
           y_updated += velocity / yaw_rate * 
                            ((-1) * cos(particles[i].theta + yaw_rate * delta_t) + cos(particles[i].theta));
           theta_updated += yaw_rate * delta_t;
        }
      
    //    theta_updated = constrainRadian(theta_updated);        
        normal_distribution<double> dist_x(x_updated, std_pos[0]);
        normal_distribution<double> dist_y(y_updated, std_pos[1]);
        normal_distribution<double> dist_theta(theta_updated, std_pos[2]);
        // update particle     
        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta  = dist_theta(gen);
    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    for (unsigned int i = 0; i < observations.size(); i++) {
        int min_id = -1;
        double curr_min_dist = -1;
        std::cout << "predicted size: " << predicted.size() << std::endl;
        for (unsigned int j = 0; j < predicted.size(); j++) {
            double distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
            if ((curr_min_dist < 0) || (curr_min_dist > distance)) {
                min_id = predicted[j].id;
                curr_min_dist = distance;
            }
       }
       observations[i].id = min_id;
    }
}

inline double prob_weight(double x, double y, double mu_x, double mu_y, double std_x, double std_y) {
    return 0.5 * exp((-1) * ((pow((x - mu_x), 2) / (2 *(std_x * std_x)))
                    + (pow((y - mu_y), 2) / (2 * (std_y * std_y))))) / (M_PI * std_x * std_y);
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
    


    std::vector<Map::single_landmark_s> landmarks = map_landmarks.landmark_list;
    for (int i = 0; i < this->num_particles; i++){
        std::vector<LandmarkObs> predicted;
        std::vector<LandmarkObs> trans_observations;
        // set predict landmark distance 
        for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
            double distance = dist(particles[i].x, particles[i].y, landmarks[j].x_f, landmarks[j].y_f);
            if (distance <= sensor_range)
            {
                LandmarkObs pred_j;
                pred_j.id = landmarks[j].id_i;
                pred_j.x = landmarks[j].x_f;
                pred_j.y = landmarks[j].y_f;
                predicted.push_back(pred_j);
            }
        }
        // transform observation to particle spaces
        for (unsigned int k = 0; k < observations.size(); k++) {
            LandmarkObs transf_obs;
            double xp = particles[i].x;
            double yp = particles[i].y;
            double theta = particles[i].theta;
            double xc = observations[k].x;
            double yc = observations[k].y;
            xp += cos(theta) * xc - sin(theta) * yc;
            yp += sin(theta) * xc + cos(theta) * yc;
            transf_obs.x = xp;
            transf_obs.y = yp;
            trans_observations.push_back(transf_obs);
        }
        dataAssociation(predicted, trans_observations);
        
        // now draw the weight
        particles[i].weight = 1.0;
        std::cout << "Landmarks size: " << landmarks.size() << std::endl;
        for(unsigned int m = 0; m < trans_observations.size(); m++) {
             int idx = trans_observations[m].id - 1;
             std::cout << "Id: " << trans_observations[m].id << std::endl;
             std::cout << "Idx: " << idx << std::endl;
             particles[i].weight *= prob_weight(trans_observations[m].x, trans_observations[m].y, 
                               landmarks[idx].x_f, landmarks[idx].y_f, std_landmark[0], std_landmark[1]);
        }
        std::cout << "Weights :" << particles[i].weight << std::endl;
    }

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    std::vector<double> new_weights;
    for (Particle p: this->particles) {
        new_weights.push_back(p.weight);
    }

    std::discrete_distribution<> d(new_weights.begin(), new_weights.end());
    default_random_engine gen;
    std::vector<Particle> new_particles;
    std::vector<double> updated_weights;
    for (int i = 0; i < this->num_particles; i++) {
        int idx = d(gen);
        new_particles.push_back(this->particles[idx]);        
        updated_weights.push_back(this->particles[idx].weight);
    }
    this->particles = new_particles;
    this->weights = new_weights;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
