#ifndef VOLUME_INTERSECTION_H_INCLUDED
#define VOLUME_INTERSECTION_H_INCLUDED

#include "mat_vec_fns_II.h"

// scalar volumeIntersection(arr3 (&tetA)[4], arr3 (&tetB)[4]);
template <class t_scalar, class brr3> scalar volumeIntersection(brr3 (&tetA)[4], brr3 (&tetB)[4]);
template <class t_scalar, class brr3> void volumeAndAreaIntersection(brr3 (&tetA)[4], brr3 (&tetB)[4], t_scalar &vol, t_scalar &area);

#endif
