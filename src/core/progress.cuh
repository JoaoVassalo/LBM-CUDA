#pragma once

#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>

#include "../config/stencilConfig.cuh"

inline void printProgress(
    size_t completedIterations,
    size_t totalIterations,
    std::chrono::steady_clock::time_point startTime)
{
    cudaDeviceSynchronize();

    const auto now = std::chrono::steady_clock::now();
    const double elapsedSeconds = std::chrono::duration<double>(now - startTime).count();
    const double completedUpdates = static_cast<double>(completedIterations) * Geometry::NX * Geometry::NY;
    const double mlups = elapsedSeconds > 0.0 ? completedUpdates / elapsedSeconds / 1'000'000.0 : 0.0;
    const double percent = totalIterations > 0
                               ? 100.0 * static_cast<double>(completedIterations) / totalIterations
                               : 100.0;
    const double remainingSeconds = completedIterations > 0
                                        ? elapsedSeconds * (totalIterations - completedIterations) / completedIterations
                                        : 0.0;

    const auto eta = static_cast<size_t>(std::ceil(remainingSeconds));
    const size_t hours = eta / 3600;
    const size_t minutes = (eta % 3600) / 60;
    const size_t seconds = eta % 60;

    std::cout << '\r'
              << std::fixed << std::setprecision(2)
              << "Iteration " << completedIterations << '/' << totalIterations
              << " | " << percent << "%"
              << " | " << mlups << " MLUPS"
              << " | ETA " << std::setfill('0') << std::setw(2) << hours
              << ':' << std::setw(2) << minutes
              << ':' << std::setw(2) << seconds
              << std::setfill(' ') << std::flush;
}
