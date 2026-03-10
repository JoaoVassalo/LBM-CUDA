#pragma once

#include "../presets/stencil.h"

__device__ inline int pop_id(int g_id, int i){
    return (g_id*Q + i);
}