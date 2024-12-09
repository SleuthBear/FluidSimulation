//
// Created by Dylan Beaumont on 7/12/2024.
//

#include "SimulationParticle2D.h"

#include <algorithm>
#include <iostream>
#include <omp.h>

SimulationParticle2D::SimulationParticle2D(int nParticles, float radius) {
    // heap allocated vector2 arrays
    pos = new Vector2[nParticles];
    predicted_pos = new Vector2[nParticles];

    vel = new Vector2[nParticles];

    spatial_lookup = new spatial_point[nParticles];

    densities = new float[nParticles];
    n = nParticles;
    smoothing_radius = radius;
    particle_radius = 10;

    int screen_width = GetScreenWidth();
    int screen_height = GetScreenHeight();

    h_cells = (int) screen_width / radius + 1;
    v_cells = (int) screen_height / radius + 1;
    start_indices = new int[h_cells*v_cells*2];

    // populate the point arrays
    for(int i=0; i<n; i++) {
        pos[i].x = screen_width/2 - screen_width/6 + rand() % screen_width/3;
        pos[i].y = screen_height/2 + -screen_height/6 + rand() % screen_height/3;
        vel[i] = {0,0};
    }
}

void SimulationParticle2D::Draw() {
    for(int i=0; i<n; i++) {
        DrawCircle(pos[i].x, pos[i].y, smoothing_radius, BLUE);
        DrawCircle(pos[i].x, pos[i].y, 2, WHITE);

    }
}

void SimulationParticle2D::Update(float delta) {

    // Gravity + density
    #pragma omp parallel for num_threads(8)
    for(int i=0; i<n; i++) {
        vel[i].y += GRAVITY*delta;
        predicted_pos[i] = pos[i] + vel[i] * (1/60);
    }

    UpdateSpatialLookup();

    #pragma omp parallel for num_threads(8)
    for(int i=0; i<n; i++) {
        std::vector<int> relevant = FindRelevantPoints(pos[i]);
        densities[i] = CalculateDensity(pos[i], relevant);
    }

    #pragma omp parallel for num_threads(8)
    for(int i=0; i<n; i++) {
        std::vector<int> relevant = FindRelevantPoints(pos[i]);
        Vector2 pressure_force = CalculatePressureForce(i, relevant);
        Vector2 pressure_acceleration = pressure_force / densities[i];
        vel[i] += pressure_acceleration;
    }

    #pragma omp parallel for num_threads(8)
    for(int i=0; i<n; i++) {
        pos[i] += vel[i] * delta;
        vel[i] *= (1-FRICTION_FACTOR);
        HandleCollision(i, delta);
    }
}

void SimulationParticle2D::HandleCollision(int i, float delta) {
    if(pos[i].y > GetScreenHeight()-particle_radius) {
        pos[i].y = GetScreenHeight()-particle_radius;
        vel[i].y = -vel[i].y;
    } else if (pos[i].y < particle_radius) {
        pos[i].y = particle_radius;
        vel[i].y = -vel[i].y;
    }
    if(pos[i].x > GetScreenWidth()-particle_radius) {
        pos[i].x = GetScreenWidth()-particle_radius;
        vel[i].x = -vel[i].x;
    } else if (pos[i].x < particle_radius) {
        pos[i].x = particle_radius;
        vel[i].x = -vel[i].x;
    }
}

float SimulationParticle2D::SmoothingKernel(float dst, float radius) {
    if(dst >= radius) return 0.f;
    float volume = PI*std::pow(radius, 4)/6;
    return (radius-dst) * (radius-dst) / volume;
}

float SimulationParticle2D::SmoothingKernelDerivative(float dst, float radius) {
    if(dst >= radius) return 0.f;
    float scale = 12/(PI*std::pow(radius,4));
    return (dst-radius) * scale;
}

float SimulationParticle2D::CalculateDensity(Vector2 point, std::vector<int> relevant) {
    float density = 0.f;
    for(int r=0; r<relevant.size(); r++) {
        int i = relevant.at(r);
        float dst = Vector2Distance(point, pos[i]);
        density += MASS*SmoothingKernel(dst, smoothing_radius);
    }
    return density;
}

Vector2 SimulationParticle2D::CalculatePressureForce(int point, std::vector<int> relevant) {
    Vector2 pressure_force = {0, 0};
    for(int r=0; r<relevant.size(); r++) {
        int i = relevant.at(r);
        if(i == point) continue;
        float dst = Vector2Distance(pos[i], pos[point]);
        Vector2 dir = (pos[i] - pos[point]) / dst;
        float slope = SmoothingKernelDerivative(dst, smoothing_radius);
        float shared_pressure = CalculateSharedPressure(densities[i], densities[point]);
        pressure_force -= dir*shared_pressure*slope*MASS / densities[i];
    }
    return pressure_force;
}

float SimulationParticle2D::CalculateSharedPressure(float density_1, float density_2) {
    float pressure_1 = ConvertDensityToPressure(density_1);
    float pressure_2 = ConvertDensityToPressure(density_2);
    return (pressure_1 + pressure_2) / 2.f;
}

float SimulationParticle2D::ConvertDensityToPressure(float density) {
    float density_error = TARGET_DENSITY - density;
    float pressure = -density_error * PRESSURE_MULTIPLIER;
    return pressure;
}

int SimulationParticle2D::PositionToCell(Vector2 point) {
    // we will create an array of cells of size radius*radius. Then find the index of the cell, and the surround cell
    // by using **MATHS**(tm)

    // e.g. if radius = 10, and y = 15, then it is in the second row (index 1), so y / radius = 1.5, casting to int -> 1.
    int row = (int) point.y / smoothing_radius;
    int column = (int) point.x / smoothing_radius;
    return row*h_cells + column;
}

void SimulationParticle2D::UpdateSpatialLookup() {
    #pragma omp parallel for num_threads(8)
    for(int i=0; i<n; i++) {
        spatial_lookup[i].cell_index = PositionToCell(predicted_pos[i]);
        spatial_lookup[i].particle_index = i;
    }
    for(int i=0; i<h_cells*v_cells; i++) {
        start_indices[i] = 9999999;
    }
    // sort by cell index
    std::sort(spatial_lookup, spatial_lookup+n, CompareSpatial);

    // get start index for each cell key
    #pragma omp parallel for num_threads(8)
    for(int i = 0; i<n; i++) {
        int key = spatial_lookup[i].cell_index;
        int prev_key = i==0 ? 999999 : spatial_lookup[i-1].cell_index;
        if(key != prev_key) {
            start_indices[key] = i;
        }
    }
}

int SimulationParticle2D::CompareSpatial(spatial_point a, spatial_point b) {
    return a.cell_index < b.cell_index;
}

std::vector<int> SimulationParticle2D::FindRelevantPoints(Vector2 point) {
    int cell = PositionToCell(point);
    std::vector<int> relevant_points;
    int cells[9];
    cells[0] = std::max(cell-1-h_cells, 0); // up and left
    cells[1] = std::max(cell-h_cells, 0); // up
    cells[2] = std::max(cell-h_cells+1, 0); // up and right
    cells[3] = std::max(cell-1, 0); // left
    cells[4] = cell; // centre;
    cells[5] = std::min(cell+1, h_cells*v_cells-1); // right
    cells[6] = std::min(cell-1+h_cells, h_cells*v_cells-1); // down left
    cells[7] = std::min(cell+h_cells, h_cells*v_cells-1); // down
    cells[8] = std::min(cell+1+h_cells, h_cells*v_cells-1); // down right
    int start_index = 0;
    for(int c=0; c<9; c++) {
        // DrawRectangleLines((cells[c] % h_cells )*smoothing_radius, (cells[c]/h_cells)*smoothing_radius, smoothing_radius, smoothing_radius, GREEN);
        start_index = start_indices[cells[c]];
        for(int i=start_index; i<n; i++) {
            if(spatial_lookup[i].cell_index != cells[c]) break;
            if(Vector2Distance(pos[spatial_lookup[i].particle_index], point) < smoothing_radius) {
                relevant_points.push_back(spatial_lookup[i].particle_index);
            }
        }
    }
    return relevant_points;
}

void SimulationParticle2D::DrawPoint(int point, Color color) {
    DrawCircle(pos[point].x, pos[point].y, 4, color);
}

void SimulationParticle2D::PushAway(Vector2 point, float radius) {
    for(int i=0; i<n; i++) {
        float dst = Vector2Distance(pos[i], point);
        if(dst < radius) {
            vel[i] += (pos[i]-point) * 500 / dst;
        }
    }
}

