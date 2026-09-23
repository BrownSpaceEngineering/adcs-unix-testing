#ifndef QUEST
#define QUEST
#include <stdbool.h>

#define QUEST_MAX_MSMTS 8

// body/ref are msmt_ct stacked 3-vectors (need not be unit length; they are normalized
// internally). Result is the body->ref quaternion [w, x, y, z]. Returns false (and writes the
// identity) if the inputs are degenerate (fewer than 2 usable, non-parallel vector pairs).
bool quest(const float* body, const float* ref, int msmt_ct, float* result);
// Same, with per-measurement weights (normalized internally; NULL means equal weights)
bool quest_weighted(const float* body, const float* ref, const float* weights, int msmt_ct,
                    float* result);
#endif
