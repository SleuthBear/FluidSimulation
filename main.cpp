#include <iostream>
#include "raylib.h"
#include "raymath.h"
#include "SimulationParticle2D.h"
void GetPoints(SimulationParticle2D sim);
int main()
{
    const int screenWidth = 1300;
    const int screenHeight = 720;
    InitWindow(screenWidth, screenHeight, "Simulation");
    SimulationParticle2D sim = SimulationParticle2D( 3000, 10.f);
    while(!WindowShouldClose()) {

        SetTargetFPS(60);
        float delta = GetFrameTime();

        BeginDrawing();

        ClearBackground(BLACK);
        if(IsMouseButtonDown(MOUSE_BUTTON_LEFT)){
            sim.PushAway(GetMousePosition(), 100);
        }
        sim.Update(delta);
        sim.Draw();
        GetPoints(sim);

        DrawText(TextFormat("CURRENT FPS: %i", (int)(1.0f/delta)), GetScreenWidth() - 220, 40, 20, GREEN);

        EndDrawing();
    }
    CloseWindow();
}

void GetPoints(SimulationParticle2D sim) {
    std::vector<int> points = sim.FindRelevantPoints(GetMousePosition());
    for(int point:points) {
        sim.DrawPoint(point, GREEN);
    }
}