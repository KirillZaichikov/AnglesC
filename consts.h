#pragma once

// this is max tax dimension that can be calculated. if task with greater dimension appears - this value should be adjusted
//     but consider that it affects device memory usage.
//     If cuda fails with "out of bound error" CUDA_STACK_SIZE may be adjusted
#define MAX_DIMENSION 7