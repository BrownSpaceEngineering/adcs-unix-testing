#ifndef SUNVEC
#define SUNVEC
// Unit Earth->Sun vector in ECI (port of sunVectorECI.m, which takes a JD instead of unix time)
void sun_vec(int unix_time, float* sun);
#endif
