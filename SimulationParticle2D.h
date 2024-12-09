//
// Created by Dylan Beaumont on 7/12/2024.
//

#ifndef SIMULATIONPARTICLE2D_H
#define SIMULATIONPARTICLE2D_H
#include <vector>

#include "raylib.h"
#include "raymath.h"

inline float GRAVITY = 200.f; //490.f;
inline float TARGET_DENSITY = 1.f;
inline float PRESSURE_MULTIPLIER = 30.f;
inline float FRICTION_FACTOR = 0.05;
inline float MASS = 1.f;

struct spatial_point {
    int cell_index;
    int particle_index;
};

class SimulationParticle2D {
    Vector2* pos;
    Vector2* predicted_pos;
    Vector2* vel;
    float* densities;
    int n;
    float smoothing_radius;
    int particle_radius;
    int h_cells;
    int v_cells;
    spatial_point *spatial_lookup;
    int* start_indices;


public:
    SimulationParticle2D(int n, float radius);
    void Draw();
    void Update(float delta);
    void HandleCollision(int i, float delta);
    static float SmoothingKernel(float value, float dst);
    static float SmoothingKernelDerivative(float radius, float dst);
    float CalculateDensity(Vector2 point, std::vector<int> relevant);
    void UpdateDensities();
    float CalculateProperty(Vector2 point);
    Vector2 CalculatePressureForce(int point, std::vector<int> relevant);
    static float ConvertDensityToPressure(float density);
    void DrawPoint(int point, Color color);

    int PositionToCell(Vector2 point);

    void UpdateSpatialLookup();

    static int CompareSpatial(spatial_point a, spatial_point b);

    std::vector<int> FindRelevantPoints(Vector2 point);

    static float CalculateSharedPressure(float density_1, float density_2);

    void PushAway(Vector2 point, float radius);
};



#endif //SIMULATIONPARTICLE2D_H
