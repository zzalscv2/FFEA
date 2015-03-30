#ifndef CONNECTIVITYTYPES_H_INCLUDED
#define CONNECTIVITYTYPES_H_INCLUDED

struct connectivity_entry {
    /* Number of connections this object has with other objects */
    int num_connections;

    /* Array of size num_connections containing all the indices of connecting objects */
    int *connected_to;
};

#endif
